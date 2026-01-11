import shutil
import subprocess
import sys
import time
import urllib.request
from contextlib import contextmanager
from pathlib import Path
from urllib.error import HTTPError, URLError

ROOT = Path(__file__).resolve().parents[1]
VIREMA_SCRIPT = ROOT / "src" / "ViReMa.py"
COMPILER_SCRIPT = ROOT / "src" / "Compiler_Module.py"
VISUALIZE_SCRIPT = ROOT / "src" / "visualize.py"

DATASET = "PR8_Giulia"  # "Test_Data", "PR8", or "PR8_Giulia"
BATCH_MODE = True  # Only supported for PR8_Giulia
SRR_LIST_FILE = ROOT / "data" / "PR8_Giulia" / "SRR_Acc_List.txt"
DOWNLOAD_DATA = True
DOWNLOAD_METHOD = "auto"  # "auto", "ena", or "sra-tools"
SRA_THREADS = "4"

PR8_READ = "1"  # used only when DATASET == "PR8" ("1" or "2")
GIULIA_READ = "1"  # used only when DATASET == "PR8_Giulia" ("1" or "2")
USE_PADDED_REFERENCE = True
OVERWRITE = True
OUTPUT_SAM = None  # default: <input_stem>.sam
EXTRA_VIREMA_ARGS = []  # e.g. ["--Seed", "25", "--MicroInDel_Length", "5"]

RUN_COMPILER = True
EXTRA_COMPILER_ARGS = []  # e.g. ["--MicroInDel_Length", "5"]
COMPILER_SUBDIR = "compiled"

RUN_VISUALIZE = True
PLOTS_SUBDIR = "plots"
PLOT_CHROMS = []  # empty means all chroms in the BED file
GENOME_LENGTH = None  # optional int for axis limits
TOP_N = 15
SHOW_PLOTS = False


class _Tee:
    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for stream in self._streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self._streams:
            stream.flush()


@contextmanager
def _redirect_streams(stdout, stderr):
    old_stdout, old_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = stdout, stderr
    try:
        yield
    finally:
        sys.stdout, sys.stderr = old_stdout, old_stderr


def _pr8_reference() -> Path:
    padded = ROOT / "data" / "PR8" / "PR8_Reference_padded.fasta"
    if USE_PADDED_REFERENCE and padded.exists():
        return padded
    return ROOT / "data" / "PR8" / "PR8_Reference.fasta"


def _dataset_paths(dataset: str) -> tuple[Path, Path, Path]:
    if dataset == "Test_Data":
        input_path = ROOT / "data" / "Test_Data" / "FHV_10k.txt"
        reference_path = ROOT / "data" / "Test_Data" / "FHV_Genome.txt"
        output_dir = ROOT / "output" / "Test_Data"
    elif dataset == "PR8":
        input_path = ROOT / "data" / "PR8" / f"SRR16862569_{PR8_READ}.fastq"
        reference_path = _pr8_reference()
        output_dir = ROOT / "output" / "PR8"
    elif dataset == "PR8_Giulia":
        input_path = ROOT / "data" / "PR8_Giulia" / f"SRR14352105_{GIULIA_READ}.fastq"
        reference_path = ROOT / "data" / "PR8_Giulia" / "PR8_AF_Reference_padded.fasta"
        output_dir = ROOT / "output" / "PR8_Giulia"
    else:
        raise ValueError('Unknown DATASET. Use "Test_Data", "PR8", or "PR8_Giulia".')
    return input_path, reference_path, output_dir


def _resolve_output_dir(base_dir: Path) -> Path:
    if OVERWRITE or not base_dir.exists():
        return base_dir
    return Path(f"{base_dir}{int(time.time())}")


def _resolve_child_dir(parent_dir: Path, name: str) -> Path:
    candidate = parent_dir / name
    if not candidate.exists():
        return candidate
    return parent_dir / f"{name}_{int(time.time())}"


def _default_sam_name(input_path: Path) -> str:
    name = input_path.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            return name[: -len(suffix)] + ".sam"
    return f"{input_path.stem}.sam"


def _load_srr_list(list_path: Path) -> list[str]:
    if not list_path.exists():
        raise FileNotFoundError(f"SRR list not found: {list_path}")
    srr_ids = []
    with list_path.open("r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            srr_ids.append(line.split()[0])
    if not srr_ids:
        raise ValueError(f"SRR list is empty: {list_path}")
    return srr_ids


def _download_file(url: str, dest: Path):
    tmp = dest.with_name(dest.name + ".part")
    if tmp.exists():
        tmp.unlink()
    try:
        with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
            shutil.copyfileobj(response, out)
        tmp.replace(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def _fetch_ena_fastq_urls(srr_id: str) -> list[str]:
    api_url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={srr_id}&result=read_run&fields=fastq_ftp&format=tsv"
    )
    try:
        with urllib.request.urlopen(api_url) as response:
            text = response.read().decode("utf-8")
    except (HTTPError, URLError) as exc:
        raise RuntimeError(f"ENA query failed for {srr_id}: {exc}") from exc

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"No ENA fastq URLs found for {srr_id}")

    fastq_field = lines[1].split("\t")[0]
    urls = []
    for entry in fastq_field.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        urls.append(entry if entry.startswith("http") else f"https://{entry}")
    if not urls:
        raise RuntimeError(f"No ENA fastq URLs found for {srr_id}")
    return urls


def _download_via_ena(srr_id: str, dest_dir: Path) -> Path:
    urls = _fetch_ena_fastq_urls(srr_id)
    read1_url = next((url for url in urls if url.endswith("_1.fastq.gz")), None)
    if read1_url is None and len(urls) == 1:
        read1_url = urls[0]
    if read1_url is None:
        raise RuntimeError(f"Could not find read1 URL for {srr_id}")

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / Path(read1_url).name
    if dest.exists():
        return dest

    print(f"Downloading {srr_id} read1 from ENA...")
    _download_file(read1_url, dest)
    return dest


def _download_via_sra_tools(srr_id: str, dest_dir: Path) -> Path:
    prefetch = shutil.which("prefetch")
    fasterq = shutil.which("fasterq-dump")
    if not prefetch or not fasterq:
        raise FileNotFoundError("sra-tools not found in PATH")

    dest_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run([prefetch, srr_id, "--output-directory", str(dest_dir)], check=True)

    sra_path = dest_dir / srr_id / f"{srr_id}.sra"
    source = str(sra_path) if sra_path.exists() else srr_id
    cmd = [fasterq, source, "--outdir", str(dest_dir), "--split-files"]
    if SRA_THREADS:
        cmd.extend(["--threads", str(SRA_THREADS)])
    subprocess.run(cmd, check=True)

    fastq = dest_dir / f"{srr_id}_1.fastq"
    if fastq.exists():
        return fastq
    fastq_single = dest_dir / f"{srr_id}.fastq"
    if fastq_single.exists():
        return fastq_single
    raise FileNotFoundError(f"Expected FASTQ output for {srr_id} in {dest_dir}")


def _ensure_fastq1(srr_id: str, dest_dir: Path) -> Path:
    fastq = dest_dir / f"{srr_id}_1.fastq"
    fastq_gz = dest_dir / f"{srr_id}_1.fastq.gz"
    fastq_single = dest_dir / f"{srr_id}.fastq"
    fastq_single_gz = dest_dir / f"{srr_id}.fastq.gz"
    if fastq.exists():
        return fastq
    if fastq_gz.exists():
        return fastq_gz
    if fastq_single.exists():
        return fastq_single
    if fastq_single_gz.exists():
        return fastq_single_gz
    if not DOWNLOAD_DATA:
        raise FileNotFoundError(f"Missing FASTQ for {srr_id} in {dest_dir}")

    if DOWNLOAD_METHOD in ("auto", "sra-tools"):
        try:
            return _download_via_sra_tools(srr_id, dest_dir)
        except FileNotFoundError:
            if DOWNLOAD_METHOD == "sra-tools":
                raise
        except subprocess.CalledProcessError as exc:
            if DOWNLOAD_METHOD == "sra-tools":
                raise
            print(f"sra-tools failed for {srr_id}: {exc}. Falling back to ENA.")

    if DOWNLOAD_METHOD in ("auto", "ena"):
        return _download_via_ena(srr_id, dest_dir)

    raise ValueError('Unknown DOWNLOAD_METHOD. Use "auto", "ena", or "sra-tools".')


def _run_and_log(cmd, cwd):
    print("Running:", " ".join(cmd))
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        sys.stdout.write(line)
    exit_code = process.wait()
    if exit_code != 0:
        raise subprocess.CalledProcessError(exit_code, cmd)


def _find_bed_file(output_dir: Path) -> Path:
    bed_dir = output_dir / "BED_Files"
    if not bed_dir.exists():
        raise FileNotFoundError(f"BED_Files directory not found: {bed_dir}")
    candidates = sorted(
        bed_dir.glob("*Virus_Recombination_Results.bed"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError(f"No Virus_Recombination_Results.bed found in {bed_dir}")
    return candidates[-1]


def _run_pipeline(
    input_path: Path,
    reference_path: Path,
    output_dir: Path,
    resolve_output_dir: bool,
    label=None,
):
    if resolve_output_dir:
        output_dir = _resolve_output_dir(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "log.txt"

    if not VIREMA_SCRIPT.exists():
        raise FileNotFoundError(f"ViReMa script not found at {VIREMA_SCRIPT}")
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference file not found: {reference_path}")

    output_sam = OUTPUT_SAM or _default_sam_name(input_path)

    virema_cmd = [
        sys.executable,
        str(VIREMA_SCRIPT),
        str(reference_path),
        str(input_path),
        output_sam,
        "--Output_Dir",
        str(output_dir),
        "-Overwrite",
    ]
    virema_cmd.extend(EXTRA_VIREMA_ARGS)

    with log_path.open("w") as log_file:
        stdout_tee = _Tee(sys.stdout, log_file)
        stderr_tee = _Tee(sys.stderr, log_file)
        with _redirect_streams(stdout_tee, stderr_tee):
            if label:
                print(f"=== Processing {label} ===")
            _run_and_log(virema_cmd, cwd=str(ROOT))

            compiler_output_dir = None
            if RUN_COMPILER:
                compiler_input = output_dir / output_sam
                if not COMPILER_SCRIPT.exists():
                    raise FileNotFoundError(f"Compiler script not found at {COMPILER_SCRIPT}")
                if not compiler_input.exists():
                    raise FileNotFoundError(f"Compiler input not found: {compiler_input}")
                compiler_output_dir = _resolve_child_dir(output_dir, COMPILER_SUBDIR)
                compiler_cmd = [
                    sys.executable,
                    str(COMPILER_SCRIPT),
                    str(compiler_input),
                    str(reference_path),
                    "--Output_Dir",
                    str(compiler_output_dir),
                    "-Overwrite",
                    "-BED",
                ]
                compiler_cmd.extend(EXTRA_COMPILER_ARGS)
                _run_and_log(compiler_cmd, cwd=str(ROOT))

            if RUN_VISUALIZE:
                if not VISUALIZE_SCRIPT.exists():
                    raise FileNotFoundError(f"Visualize script not found at {VISUALIZE_SCRIPT}")
                bed_source_dir = compiler_output_dir if RUN_COMPILER else output_dir
                bed_file = _find_bed_file(bed_source_dir)
                plots_dir = output_dir / PLOTS_SUBDIR
                visualize_cmd = [
                    sys.executable,
                    str(VISUALIZE_SCRIPT),
                    str(bed_file),
                    "--output-dir",
                    str(plots_dir),
                    "--top-n",
                    str(TOP_N),
                ]
                if GENOME_LENGTH:
                    visualize_cmd.extend(["--genome-length", str(GENOME_LENGTH)])
                for chrom in PLOT_CHROMS:
                    visualize_cmd.extend(["--chrom", chrom])
                if SHOW_PLOTS:
                    visualize_cmd.append("--show")
                _run_and_log(visualize_cmd, cwd=str(ROOT))


def main():
    if BATCH_MODE:
        if DATASET != "PR8_Giulia":
            raise ValueError("BATCH_MODE is only supported for DATASET == 'PR8_Giulia'.")
        srr_ids = _load_srr_list(SRR_LIST_FILE)
        data_dir = ROOT / "data" / "PR8_Giulia"
        reference_path = ROOT / "data" / "PR8_Giulia" / "PR8_AF_Reference_padded.fasta"
        base_output_dir = ROOT / "output" / "PR8_Giulia"
        for srr_id in srr_ids:
            input_path = _ensure_fastq1(srr_id, data_dir)
            output_dir = base_output_dir / srr_id
            _run_pipeline(
                input_path,
                reference_path,
                output_dir,
                resolve_output_dir=False,
                label=srr_id,
            )
        return

    input_path, reference_path, output_dir = _dataset_paths(DATASET)
    _run_pipeline(input_path, reference_path, output_dir, resolve_output_dir=True)


if __name__ == "__main__":
    main()
