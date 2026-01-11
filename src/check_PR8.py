import os

# --- Config ---
REF_FILE = "data/PR8/PR8_Reference_padded.fasta"
FASTQ_FILE = "data/PR8/SRR8526463_1.fastq"

def check_files():
    print(f"--- DIAGNOSTIC REPORT ---")
    
    # 1. Check if files exist and have size
    if not os.path.exists(REF_FILE):
        print(f"[FAIL] Reference file missing: {REF_FILE}")
        return
    if not os.path.exists(FASTQ_FILE):
        print(f"[FAIL] FASTQ file missing: {FASTQ_FILE}")
        return

    ref_size = os.path.getsize(REF_FILE)
    fastq_size = os.path.getsize(FASTQ_FILE)
    print(f"Reference Size: {ref_size/1024:.2f} KB")
    print(f"FASTQ Size:     {fastq_size/1024/1024:.2f} MB")
    
    if ref_size < 1000:
        print("[FAIL] Reference file seems too small (empty?).")
        return

    # 2. Load Reference Genome
    print("\nLoading Reference Genome...")
    genome_seqs = {}
    current_header = ""
    seq_buffer = []
    
    with open(REF_FILE, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_header:
                    genome_seqs[current_header] = "".join(seq_buffer)
                current_header = line.strip().split()[0] # Clean header
                seq_buffer = []
            else:
                seq_buffer.append(line.strip())
        if current_header:
            genome_seqs[current_header] = "".join(seq_buffer)
            
    print(f"Loaded {len(genome_seqs)} segments.")
    for name, seq in genome_seqs.items():
        print(f"  > {name}: {len(seq)} bp")

    # 3. Check Reads
    print(f"\nChecking first 50 reads from {FASTQ_FILE}...")
    matches = 0
    checked = 0
    
    with open(FASTQ_FILE, 'r') as f:
        lines = []
        for line in f:
            lines.append(line.strip())
            if len(lines) == 4:
                checked += 1
                seq = lines[1]
                
                # Check normal and reverse complement
                rev_seq = seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
                
                found = False
                for ref_name, ref_seq in genome_seqs.items():
                    if seq[:30] in ref_seq: # Check just the first 30bp (Seed)
                        found = True
                        break
                    if rev_seq[:30] in ref_seq:
                        found = True
                        break
                
                if found:
                    matches += 1
                else:
                    if checked <= 3: # Print failure for first few only
                        print(f"  [FAIL] Read {checked} start '{seq[:15]}...' NOT found in genome.")
                
                lines = []
                if checked >= 50: break
    
    print(f"\n--- RESULT ---")
    print(f"Checked {checked} reads.")
    print(f"Direct Matches Found: {matches}")
    
    if matches == 0:
        print("CONCLUSION: Your reads do NOT match this Reference. The data is likely from a different organism or the file is corrupt.")
    elif matches < 10:
        print("CONCLUSION: Very low match rate. Sequencing quality might be poor.")
    else:
        print("CONCLUSION: Data looks GOOD! The issue is likely in the Bowtie Indexing step.")

if __name__ == "__main__":
    check_files()