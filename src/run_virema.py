import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VIREMA_SCRIPT = ROOT / "src" / "ViReMa.py"
COMPILER_SCRIPT = ROOT / "src" / "Compiler_Module.py"
VISUALIZE_SCRIPT = ROOT / "src" / "visualize.py"

DATASET = "PR8_Giulia"  # "Test_Data", "PR8", or "PR8_Giulia"
PR8_READ = "1"  # used only when DATASET == "PR8" ("1" or "2")
GIULIA_READ = "1"  # used only when DATASET == "PR8_Giulia" ("1" or "2")
USE_PADDED_REFERENCE = True
OVERWRITE = True
OUTPUT_SAM = None  # default: <input_stem>.sam
EXTRA_VIREMA_ARGS = []  # e.g. ["--Seed", "25", "--MicroInDel_Length", "5"]

RUN_COMPILER = True
EXTRA_COMPILER_ARGS = []  # e.g. ["--MicroInDel_Length", "5"]
COMPILER_SUBDIR = "compiled"

RUN_VISUALIZE = Truew
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


def main():
    input_path, reference_path, output_dir = _dataset_paths(DATASET)
    output_dir = _resolve_output_dir(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "log.txt"

    if not VIREMA_SCRIPT.exists():
        raise FileNotFoundError(f"ViReMa script not found at {VIREMA_SCRIPT}")
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference file not found: {reference_path}")

    output_sam = OUTPUT_SAM or f"{input_path.stem}.sam"

    virema_cmd = [
        sys.executable,
        str(VIREMA_SCRIPT),
        str(reference_path),
        str(input_path),
        output_sam,
        "--Output_Dir",
        str(output_dir),
    ]
    # Force ViReMa to use the exact output_dir we prepared.
    virema_cmd.append("-Overwrite")
    virema_cmd.extend(EXTRA_VIREMA_ARGS)

    with log_path.open("w") as log_file:
        stdout_tee = _Tee(sys.stdout, log_file)
        stderr_tee = _Tee(sys.stderr, log_file)
        with _redirect_streams(stdout_tee, stderr_tee):
            _run_and_log(virema_cmd, cwd=str(ROOT))

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


if __name__ == "__main__":
    main()
