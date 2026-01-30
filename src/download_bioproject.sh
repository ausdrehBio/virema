#!/usr/bin/env zsh

INPUT_DIR="B_Lee_1940"
BASE_DIR="data/${INPUT_DIR}"
BIOPROJECT_FILE="${BASE_DIR}/Bioproject_ID.txt"
RUNINFO_CSV="${BASE_DIR}/${INPUT_DIR}_runinfo.csv"
SRR_LIST_FILE="${BASE_DIR}/SRR_Acc_List.txt"
RUNINFO_TMP="${BASE_DIR}/.runinfo_tmp.csv"

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: Base directory not found: $BASE_DIR"
    exit 1
fi

if [[ ! -f "$BIOPROJECT_FILE" ]]; then
    echo "ERROR: Missing file: $BIOPROJECT_FILE"
    exit 1
fi

for cmd in esearch elink efetch prefetch fasterq-dump; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "wrong env: conda activate virema_stable"
        exit 1
    fi
done

BIOPROJECT_IDS=()
while IFS= read -r raw; do
    line="${raw%%#*}"
    line="${line//[[:space:]]/}"
    if [[ -n "$line" ]]; then
        BIOPROJECT_IDS+=("$line")
    fi
done < "$BIOPROJECT_FILE"

if (( ${#BIOPROJECT_IDS[@]} == 0 )); then
    echo "ERROR: $BIOPROJECT_FILE is empty."
    exit 1
fi

: > "$RUNINFO_CSV"
: > "$SRR_LIST_FILE"

header_written=0
total_srr=0

for BIOPROJECT_ID in "${BIOPROJECT_IDS[@]}"; do
    BIOPROJECT_DIR="${BASE_DIR}/${BIOPROJECT_ID}"
    mkdir -p "$BIOPROJECT_DIR"

    echo "Fetching metadata for $BIOPROJECT_ID..."
    esearch -db bioproject -query "$BIOPROJECT_ID" | \
        elink -target sra | \
        efetch -format runinfo > "$RUNINFO_TMP"

    SRR_LIST=($(cut -d "," -f 1 "$RUNINFO_TMP" | grep "^SRR"))

    if (( ${#SRR_LIST[@]} == 0 )); then
        echo "WARNING: No SRR runs found for $BIOPROJECT_ID."
        continue
    fi

    if (( header_written == 0 )); then
        cat "$RUNINFO_TMP" >> "$RUNINFO_CSV"
        header_written=1
    else
        tail -n +2 "$RUNINFO_TMP" >> "$RUNINFO_CSV"
    fi

    for SRR_ID in "${SRR_LIST[@]}"; do
        print -r -- "$SRR_ID" >> "$SRR_LIST_FILE"
    done

    echo "Found ${#SRR_LIST[@]} runs for $BIOPROJECT_ID. Starting download..."
    for SRR_ID in "${SRR_LIST[@]}"; do
        echo "Downloading $SRR_ID..."
        if prefetch "$SRR_ID" --output-directory "$BIOPROJECT_DIR"; then
            sra_path="${BIOPROJECT_DIR}/${SRR_ID}/${SRR_ID}.sra"
            if [[ -f "$sra_path" ]]; then
                fasterq-dump --split-files --outdir "$BIOPROJECT_DIR" "$sra_path"
            else
                fasterq-dump --split-files --outdir "$BIOPROJECT_DIR" "$SRR_ID"
            fi
        else
            echo "WARNING: Prefetch failed for $SRR_ID"
        fi
    done

    total_srr=$(( total_srr + ${#SRR_LIST[@]} ))
done

rm -f "$RUNINFO_TMP"

if (( total_srr == 0 )); then
    echo "ERROR: No SRR runs found for any BioProject."
    exit 1
fi

echo "Done."
echo "RunInfo CSV: $RUNINFO_CSV"
echo "SRR list: $SRR_LIST_FILE"
