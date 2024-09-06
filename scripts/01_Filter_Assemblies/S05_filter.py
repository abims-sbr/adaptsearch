import sys
import os
import logging
from Bio import SeqIO

# Configuration du logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

if len(sys.argv) != 3:
    logging.error("Usage: S05_filter.py <input_fasta> <length_seq_min>")
    sys.exit(1)

path_IN = sys.argv[1]
length_seq_min = int(sys.argv[2])
path_OUT = path_IN + "_filtered"

logging.info(f"Filtering sequences in {path_IN} with max length {length_seq_min}")

try:
    with open(path_IN, "r") as f_in:
        sequences = list(SeqIO.parse(f_in, "fasta"))
except FileNotFoundError:
    logging.error(f"Input file {path_IN} not found.")
    sys.exit(1)

if not sequences:
    logging.warning(f"No sequences found in {path_IN}.")
else:
    with open(path_OUT, "w") as f_out:
        count = 0
        for record in sequences:
            #print(f"len(record.seq) {len(record.seq)}")
            if len(record.seq) >= length_seq_min:
                SeqIO.write(record, f_out, "fasta")
                count += 1
        logging.info(f"{count} sequences written to {path_OUT}.")

if os.path.getsize(path_OUT) == 0:
    logging.warning(f"Filtered output file {path_OUT} is empty.")
else:
    logging.info(f"Filtered output file {path_OUT} contains data.")
