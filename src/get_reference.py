import os
import urllib.request
from pathlib import Path

pathto = "data/PR8_Giulia/reference_accessions.txt"

def load_accession_ids(path):
    ids = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.append(line.split()[0])
    return ids

def prepare_reference():
    ids = load_accession_ids(pathto)
    if not ids:
        raise ValueError(f"No accession IDs found in {pathto}")
    
    output_dir = Path(pathto).parent
    output_file = f"{output_dir}/PR8_AF_Reference_padded.fasta"
    os.makedirs(output_dir, exist_ok=True)
    
    print("Downloading Reference Genome...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(ids)}&rettype=fasta&retmode=text"
    
    try:
        with urllib.request.urlopen(url) as response:
            raw_data = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error downloading: {e}")
        return

    print(f"Processing to {output_file}...")
    with open(output_file, "w") as f:
        entries = raw_data.split(">")
        count = 0
        for entry in entries:
            if not entry.strip(): continue
            
            lines = entry.split("\n")
            # CLEAN HEADER: Keep only the ID (e.g., AF389115.1), drop the rest
            clean_header = lines[0].split()[0]
            
            # PAD SEQUENCE: Rejoin lines and add 160 'A's
            sequence = "".join(lines[1:])
            padded_seq = sequence + ("A" * 160)
            
            f.write(f">{clean_header}\n{padded_seq}\n")
            print(f"  > Saved {clean_header}")
            count += 1
            
    print(f"Done! {count} segments ready.")

if __name__ == "__main__":
    prepare_reference()
