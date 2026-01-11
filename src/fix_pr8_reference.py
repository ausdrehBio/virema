import urllib.request
import sys
import os

def fix_pr8_reference():
    output_file = "data/PR8/PR8_Reference_padded.fasta"
    
    # 1. Define the 8 Segments of PR8
    segments = {
        "NC_002023.1": "Segment 1 (PB2)",
        "NC_002021.1": "Segment 2 (PB1)",
        "NC_002022.1": "Segment 3 (PA)",
        "NC_002017.1": "Segment 4 (HA)",    
        "NC_002016.1": "Segment 5 (NP)",
        "NC_002018.1": "Segment 6 (NA)",
        "NC_002019.1": "Segment 7 (M)",
        "NC_002020.1": "Segment 8 (NS)"
    }
    
    print("Downloading all 8 segments of Influenza PR8...")
    
    # Construct NCBI URL
    ids = ",".join(segments.keys())
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids}&rettype=fasta&retmode=text"
    
    try:
        with urllib.request.urlopen(url) as response:
            raw_fasta = response.read().decode('utf-8')
    except Exception as e:
        print(f"Error downloading: {e}")
        return

    # 2. Process, Clean, and Pad
    print(f"Processing and writing to {output_file}...")
    
    with open(output_file, "w") as f_out:
        # Split by '>' to handle each sequence
        entries = raw_fasta.split(">")
        count = 0
        found_segments = []
        
        for entry in entries:
            if not entry.strip(): continue
            
            # Separate Header from Sequence
            lines = entry.split("\n")
            raw_header = lines[0]
            # Clean Header: Take only the ID (e.g., NC_002023.1)
            clean_header = raw_header.split()[0]
            
            # Reassemble Sequence (remove newlines)
            sequence = "".join(lines[1:])
            
            # Add Padding
            padded_sequence = sequence + ("A" * 160)
            
            # Write to file
            f_out.write(f">{clean_header}\n{padded_sequence}\n")
            
            if clean_header in segments:
                found_segments.append(clean_header)
                count += 1

    # 3. Final Verification
    print("\n--- Verification Report ---")
    missing = [s for s in segments if s not in found_segments]
    
    if len(missing) == 0:
        print(f"[SUCCESS] All 8 Segments found and written to {output_file}")
    else:
        print(f"[ERROR] Missing segments: {missing}")
        print("Please check your internet connection and try again.")

if __name__ == "__main__":
    # Ensure directory exists
    os.makedirs("data/PR8", exist_ok=True)
    fix_pr8_reference()