# virema

Simple pipeline around ViReMa to download SRR data, run mapping, compile results, and generate plots.

## Requirements
- create conda env: conda create -n virema_stable python=3.10 bowtie
- Python 3.10+
- bowtie in PATH (bowtie, bowtie-build, bowtie-inspect)
- Python packages: numpy, pandas, matplotlib
- Optional: sra-tools (prefetch, fasterq-dump) for downloads

## Data layout

Prepare a subfolder under `data/` with:

- `Reference_padded.fasta`
- `SRR_Acc_List.txt`
- (downloaded) `SRR*_1.fastq` or `SRR*_1.fastq.gz`

Only `SRR_Acc_List.txt` (and `reference_accessions.txt` if you use it) are tracked in git; FASTQ files are ignored.

## Create a reference (optional helper)

`src/get_reference.py` downloads accessions listed in `data/<DATA_SUBDIR>/reference_accessions.txt` and writes
`Reference_padded.fasta` in the same folder.

Steps:

1. Set `pathto` in `src/get_reference.py` to your accession list.
2. Run: `python src/get_reference.py`

## Run the pipeline

Edit `src/run_virema.py`:

- `DATA_SUBDIR` - subfolder under `data/` (e.g. `PR8`)
- `BATCH_MODE` - `True` to process all SRR IDs, `False` for a single run
- `SINGLE_SRR` - set an SRR ID if `BATCH_MODE` is `False`
- `DOWNLOAD_METHOD` - `auto`, `ena`, or `sra-tools`

Run:

```
python src/run_virema.py
```

## Output structure

For each SRR in batch mode:

- `output/<DATA_SUBDIR>/<SRR_ID>/log.txt`
- `output/<DATA_SUBDIR>/<SRR_ID>/compiled*/` (compiler outputs, BED files)
- `output/<DATA_SUBDIR>/<SRR_ID>/plots/` (PNG plots)

## Notes

- Only read1 (`*_1.fastq`) is processed.
- Downloads are saved under `data/<DATA_SUBDIR>/`.
