import sys

def clean_and_pad(input_file, output_file, pad_length=160):
    print(f"Processing {input_file}...")
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = ""
        sequence = []
        
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                # Write previous sequence if it exists
                if header:
                    full_seq = "".join(sequence)
                    padded_seq = full_seq + ("A" * pad_length)
                    f_out.write(f"{header}\n{padded_seq}\n")
                    print(f"Fixed: {header}")
                
                # CLEAN THE NEW HEADER
                # Take everything after '>' but STOP at the first space
                # Example: ">NC_002016.1 Influenza..." becomes ">NC_002016.1"
                raw_header = line.split()[0]
                header = raw_header
                sequence = []
            else:
                sequence.append(line)
        
        # Write the final sequence
        if header:
            full_seq = "".join(sequence)
            padded_seq = full_seq + ("A" * pad_length)
            f_out.write(f"{header}\n{padded_seq}\n")
            print(f"Fixed: {header}")


# Usage: python Pad_Genome.py input.fasta output_padded.fasta
input = "data/PR8/GCF_000865725.1_ViralMultiSegProj15521_genomic.fna"
output = "data/PR8/PR8_Reference_padded.fasta"

clean_and_pad(input, output)