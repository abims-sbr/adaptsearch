#!/usr/bin/env python
import os
import subprocess
import sys

from Bio import SeqIO


def fasta_formatter(input_file, output_file):
    """
    Reformats the input FASTA file to ensure that sequences
    are on a single line.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ''
        header = ''
        for line in infile:
            if line.startswith('>'):
                if sequence:
                    outfile.write(sequence + '\n')
                header = line.strip()
                outfile.write(header + '\n')
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            outfile.write(sequence + '\n')


def reformat_headers(input_file, output_file, prefix):
    """
    Reformats the headers of the FASTA records by adding a specified prefix
    and ensures that sequences are on a single line.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ''
        for line in infile:
            if line.startswith('>'):
                if sequence:
                    outfile.write(sequence + '\n')
                # Process header line
                original_id = line[1:].strip()
                header_parts = original_id.split('/')
                numeric_part = header_parts[0].replace('ou', '')
                rest = '/'.join(header_parts[1:]) \
                    if len(header_parts) > 1 else ""
                if rest:
                    new_header = ">{}{}{}".format(prefix, numeric_part, f"/{rest}")
                else:
                    new_header = ">{}{}".format(prefix, numeric_part)
                outfile.write(new_header + '\n')
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            outfile.write(sequence + '\n')


def rename_fasta_headers(input_fasta, output_fasta):
    # Extraire le nom de base du fichier (sans l'extension .fasta)
    base_name_dir = input_fasta.split('.')[0]
    base_name = base_name_dir.split('/')[1]
    # Les deux premières lettres du nom de fichier
    prefix = base_name[3:5]
    # Liste pour stocker les nouvelles séquences
    modified_sequences = []

    # Lire le fichier fasta et modifier les en-têtes
    for index, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
        # Longueur de la séquence
        seq_length = len(record.seq)

        # Créer le nouvel en-tête selon le format demandé
        new_header = f">{prefix}{index}_1/1_1.000_{seq_length}"

        # Modifier l'en-tête du record
        record.id = new_header[1:]  # [1:] pour enlever le ">"
        record.description = ""

        # Ajouter la séquence modifiée à la liste
        modified_sequences.append(record)

    # Ecrire le fichier de sortie avec les nouvelles en-têtes
    SeqIO.write(modified_sequences, output_fasta, "fasta")


def main():
    if len(sys.argv) < 5:
        print(
            "Usage: script.py <input_files> <length_seq_min>",
            "<percent_identity> <overlap_length>")
        sys.exit(1)

    output_dir = "outputs_new"  # Define the output directory
    os.makedirs(output_dir, exist_ok=True)
    percent_identity = sys.argv[3]
    overlap_length = sys.argv[4]

    for name in sys.argv[1].split(","):
        if not os.path.isfile(name):
            print(f"Error: Input file {name} does not exist.")
            continue

        # Apply CAP3
        # Get the base file name
        file_name = os.path.basename(name)
        # Define the output file path in the output directory
        output_file_path = os.path.join(output_dir, file_name)
        # Create a symbolic link for the input file in the output directory
        symlink_path = os.path.join(output_dir, file_name)
        if not os.path.exists(symlink_path):
            os.symlink(os.path.abspath(name), symlink_path)

        # Print and run the CAP3 command
        print(
            f"cap3 {output_file_path} -p {percent_identity}"
            " -o {overlap_length}")
        subprocess.run([
            "cap3", output_file_path, "-p", percent_identity,
            "-o", overlap_length], check=True)

        # Format file to have sequence in one line
        name_fasta_formatter = os.path.join(
            output_dir, f"02_{os.path.basename(name)}")
        fasta_formatter(
            f"{output_file_path}.cap.singlets", name_fasta_formatter)

        # Merge singlets and contigs
        merged_file = os.path.join(output_dir, f"03_{file_name}_merged.fasta")
        # Define paths for CAP3 output files
        cap_singlets_file = os.path.join(
            output_dir, f"{file_name}.cap.singlets")
        cap_contigs_file = os.path.join(output_dir, f"{file_name}.cap.contigs")
        print(f"{cap_singlets_file} and {cap_contigs_file}")

        with open(merged_file, 'w') as outfile:
            # Write the contents of the contigs file first
            if os.path.exists(cap_contigs_file):
                with open(cap_contigs_file, 'r') as contigs:
                    outfile.write(contigs.read())
            # Append the contents of the singlets file
            if os.path.exists(cap_singlets_file):
                with open(cap_singlets_file, 'r') as singlets:
                    outfile.write(singlets.read())

        # Reformat headers
        name_fasta_final = os.path.join(
            output_dir, f"04_{os.path.basename(name)}")
        rename_fasta_headers(merged_file, name_fasta_final)

        # Format final file to have sequence in one line
        prefix = file_name[:2]
        tmp = prefix + os.path.basename(name)
        name_final_file = os.path.join(output_dir, tmp)
        fasta_formatter(name_fasta_final, name_final_file)


if __name__ == "__main__":
    main()
