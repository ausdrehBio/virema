#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess
import sys


INPUT_DIR = "A_PuertoRico_8_1934"
BASE_DIR = Path("data") / INPUT_DIR
BIOPROJECT_FILE = BASE_DIR / "Bioproject_ID.txt"
RUNINFO_CSV = BASE_DIR / f"{INPUT_DIR}_runinfo.csv"
SRR_LIST_FILE = BASE_DIR / "SRR_Acc_List.txt"
RUNINFO_TMP = BASE_DIR / ".runinfo_tmp.csv"


def _require_tools(tools):
    for tool in tools:
        if shutil.which(tool) is None:
            print("wrong env: conda activate virema_stable")
            sys.exit(1)


def _load_bioproject_ids(path):
    ids = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.split("#", 1)[0]
            line = "".join(line.split())
            if line:
                ids.append(line)
    return ids


def _run_pipeline(cmds, output_path):
    processes = []
    prev_stdout = None
    for index, cmd in enumerate(cmds):
        is_last = index == len(cmds) - 1
        if is_last:
            with output_path.open("wb") as out:
                proc = subprocess.Popen(cmd, stdin=prev_stdout, stdout=out)
                processes.append(proc)
                if prev_stdout is not None:
                    prev_stdout.close()
                for process in processes:
                    process.wait()
                    if process.returncode != 0:
                        raise subprocess.CalledProcessError(process.returncode, process.args)
            return
        proc = subprocess.Popen(cmd, stdin=prev_stdout, stdout=subprocess.PIPE)
        processes.append(proc)
        if prev_stdout is not None:
            prev_stdout.close()
        prev_stdout = proc.stdout


def _extract_srr_ids(runinfo_path):
    srr_ids = []
    with runinfo_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            first = line.split(",", 1)[0].strip()
            if first.startswith("SRR"):
                srr_ids.append(first)
    return srr_ids


def main():
    if not BASE_DIR.is_dir():
        print(f"ERROR: Base directory not found: {BASE_DIR}")
        return 1
    if not BIOPROJECT_FILE.is_file():
        print(f"ERROR: Missing file: {BIOPROJECT_FILE}")
        return 1

    _require_tools(["esearch", "elink", "efetch", "prefetch", "fasterq-dump"])

    bioproject_ids = _load_bioproject_ids(BIOPROJECT_FILE)
    if not bioproject_ids:
        print(f"ERROR: {BIOPROJECT_FILE} is empty.")
        return 1

    RUNINFO_CSV.write_text("", encoding="utf-8")
    SRR_LIST_FILE.write_text("", encoding="utf-8")

    header_written = False
    total_srr = 0

    try:
        for bioproject_id in bioproject_ids:
            bioproject_dir = BASE_DIR / bioproject_id
            bioproject_dir.mkdir(parents=True, exist_ok=True)

            print(f"Fetching metadata for {bioproject_id}...")
            _run_pipeline(
                [
                    ["esearch", "-db", "bioproject", "-query", bioproject_id],
                    ["elink", "-target", "sra"],
                    ["efetch", "-format", "runinfo"],
                ],
                RUNINFO_TMP,
            )

            srr_ids = _extract_srr_ids(RUNINFO_TMP)
            if not srr_ids:
                print(f"WARNING: No SRR runs found for {bioproject_id}.")
                continue

            with RUNINFO_CSV.open("a", encoding="utf-8") as out:
                with RUNINFO_TMP.open("r", encoding="utf-8") as inp:
                    for index, line in enumerate(inp):
                        if index == 0 and header_written:
                            continue
                        out.write(line)
            header_written = True

            with SRR_LIST_FILE.open("a", encoding="utf-8") as out:
                for srr_id in srr_ids:
                    out.write(f"{srr_id}\n")

            print(f"Found {len(srr_ids)} runs for {bioproject_id}. Starting download...")
            for srr_id in srr_ids:
                print(f"Downloading {srr_id}...")
                prefetch = subprocess.run(
                    ["prefetch", srr_id, "--output-directory", str(bioproject_dir)]
                )
                if prefetch.returncode != 0:
                    print(f"WARNING: Prefetch failed for {srr_id}")
                    continue

                sra_path = bioproject_dir / srr_id / f"{srr_id}.sra"
                if sra_path.is_file():
                    fasterq_cmd = [
                        "fasterq-dump",
                        "--split-files",
                        "--outdir",
                        str(bioproject_dir),
                        str(sra_path),
                    ]
                else:
                    fasterq_cmd = [
                        "fasterq-dump",
                        "--split-files",
                        "--outdir",
                        str(bioproject_dir),
                        srr_id,
                    ]
                fasterq = subprocess.run(fasterq_cmd)
                if fasterq.returncode != 0:
                    print(f"WARNING: fasterq-dump failed for {srr_id}")

            total_srr += len(srr_ids)
    finally:
        if RUNINFO_TMP.exists():
            RUNINFO_TMP.unlink()

    if total_srr == 0:
        print("ERROR: No SRR runs found for any BioProject.")
        return 1

    print("Done.")
    print(f"RunInfo CSV: {RUNINFO_CSV}")
    print(f"SRR list: {SRR_LIST_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
