import argparse
from Bio import SeqIO
import csv
import re

parser = argparse.ArgumentParser(
    description="Replace labels in CDS FASTA with metadata from CSV"
)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="<input FASTA> Sequences to relabel")
parser.add_argument("csvfile", type=argparse.FileType('r'),
                    help="<input CSV> File with metadata")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output FASTA>  File to write relabeled sequences.")
args = parser.parse_args()

keys = ['accn', 'strain', 'serotype', 'host', 'country', 'coldate']
metadata = {}
reader = csv.DictReader(args.csvfile)
for row in reader:
    metadata.update({row['accn']: row})


records = SeqIO.parse(args.infile, "fasta")
for record in records:
    accn = record.description.split(".")[0].replace('lcl|', '')
    if accn not in metadata:
        args.outfile.write(f">{accn}__\n{record.seq}\n")
    else:
        md = metadata[accn]
        label = '_'.join([md[k] for k in keys])
        args.outfile.write(f">{label}\n{record.seq}\n")        
