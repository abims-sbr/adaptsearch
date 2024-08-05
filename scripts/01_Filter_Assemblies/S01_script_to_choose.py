import sys
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def reformat_headers(input_file, output_file, prefix):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                parts = line.strip().split('_')
                if len(parts) > 1:
                    # Extract the numeric part after 'ou'
                    numeric_part = parts[0][1:].replace('ou', '')
                    rest = '_'.join(parts[1:])
                    new_header = f">{prefix}{numeric_part}_{rest}"
                    outfile.write(new_header + '\n')
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

def fasta_formatter(input_file, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ''
        header = ''
        for line in infile:
            if line.startswith('>'):
                if header:
                    outfile.write(header + '\n')
                    outfile.write(sequence + '\n')
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
        if header:
            outfile.write(header + '\n')
            outfile.write(sequence + '\n')

def name_formatting(name, prefix):
    if not os.path.isfile(name):
        print(f"Error: File {name} does not exist.")
        return ""

    with open(name, "r") as f:
        f1 = f.readline()  # Only need to check first line to know the assembler which has been used

    name_find_orf_input = ""

    if f1.startswith(">Locus"):
        name_remove_redondancy = os.path.join("outputs", f"02_{name}")
        os.makedirs(os.path.dirname(name_remove_redondancy), exist_ok=True)
        subprocess.run(["python", "S02a_remove_redondancy_from_velvet_oases.py", name, name_remove_redondancy], check=True)
        name_find_orf_input = os.path.join("outputs", f"{prefix}{name}")
        sed_cmd = f"sed -e 's/Locus_/{prefix}/g' -e 's/_Confidence_/_/g' -e 's/_Transcript_/_/g' -e 's/_Length_/_/g' {name_remove_redondancy}"
        sed_output = subprocess.check_output(sed_cmd, shell=True, text=True)
        with open(name_find_orf_input, 'w') as file:
            file.write(sed_output)
    elif f1.startswith(">c"):
        # Format the name of the sequences with good name
        name_format_fasta = os.path.join("outputs", f"03_{name}")
        os.makedirs(os.path.dirname(name_format_fasta), exist_ok=True)
        subprocess.run(["python", "S02b_format_fasta_name_trinity.py", name, name_format_fasta, prefix], check=True)
        # Apply first script to avoid redundant sequences
        name_find_orf_input = os.path.join("outputs", f"04_{name}")
        os.makedirs(os.path.dirname(name_find_orf_input), exist_ok=True)
        subprocess.run(["python", "S03_choose_one_variants_per_locus_trinity.py", name_format_fasta, name_find_orf_input], check=True)

    return name_find_orf_input

def main():
    if len(sys.argv) < 5:
        print("Usage: script.py <input_files> <length_seq_max> <percent_identity> <overlap_length>")
        sys.exit(1)

    os.makedirs("outputs", exist_ok=True)
    length_seq_max = sys.argv[2]
    percent_identity = sys.argv[3]
    overlap_length = sys.argv[4]

    for name in sys.argv[1].split(","):
        if not os.path.isfile(name):
            print(f"Error: Input file {name} does not exist.")
            continue

        prefix = name[0:2]
        name_fasta_formatter = os.path.join("outputs", f"01_{name}")
        fasta_formatter(name, name_fasta_formatter)
        name_find_orf_input = name_formatting(name_fasta_formatter, prefix)

        if name_find_orf_input == "":
            continue

        # Apply the ORF finding script to keep the longest ORF
        name_find_orf = os.path.join("outputs", f"05_{name}")
        os.makedirs(os.path.dirname(name_find_orf), exist_ok=True)
        subprocess.run(["python", "S04_find_orf.py", name_find_orf_input, name_find_orf], check=True)

        # Apply CAP3
        subprocess.run(["cap3", name_find_orf, "-p", percent_identity, "-o", overlap_length], check=True)

        # Merge singlets and contigs
        singlets_output = subprocess.check_output(["zcat", "-f", f"{name_find_orf}.cap.singlets"], text=True)
        singlets_output_file = os.path.join("outputs", f"{prefix}_singlets.fasta")
        with open(singlets_output_file, 'w') as file:
            file.write(singlets_output)

        # Apply filter script
        name_filter = os.path.join("outputs", f"05_{name}")
        subprocess.run(["python", "S05_filter.py", singlets_output_file, length_seq_max], check=True)

        # Reformat headers in the final output file
        final_output_file = os.path.join("outputs", f"{prefix}{name}")
        reformat_headers(name_filter, final_output_file, prefix)

if __name__ == "__main__":
    main()
