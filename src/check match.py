#check match
def check_match(fasta_path, fastq_path):
    print(f"Checking if reads from {fastq_path} exist in {fasta_path}...")
    
    # 1. Load Genome
    genome = ""
    with open(fasta_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                genome += line.strip()
    
    # 2. Check first 10 reads
    with open(fastq_path, 'r') as f:
        count = 0
        for line in f:
            if count >= 40: break # Check first 10 reads (4 lines per read)
            
            # Line 2 is the sequence
            if count % 4 == 1:
                seq = line.strip()
                # Check for match (Normal or Reverse Complement)
                rev_seq = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
                
                if seq in genome:
                    print(f"[SUCCESS] Read found in genome! (Forward match)")
                    return
                elif rev_seq in genome:
                    print(f"[SUCCESS] Read found in genome! (Reverse match)")
                    return
                else:
                    print(f"[FAIL] Read starting with {seq[:10]}... not found.")
            count += 1
            
    print("\nResult: No direct matches found in first 10 reads.")
    print("If you see [SUCCESS], your data is good and ViReMa/Bowtie is the issue.")
    print("If you see [FAIL], your reads do not belong to this reference genome.")

if __name__ == "__main__":
    # Adjust paths if necessary
    check_match("data/PR8/PR8_Reference_padded.fasta", "data/PR8/SRR16862569_1.fastq")