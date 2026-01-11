import urllib.request
import os

def prepare_reference():
    # Your specific Accession List
    ids = [
        "AF389115.1",
        "AF389116.1",
        "AF389117.1",
        "AF389118.1",
        "AF389119.1",
        "AF389120.1",
        "AF389121.1",
        "AF389122.1"
    ]
    
    output_dir = "data/PR8_Giulia"
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