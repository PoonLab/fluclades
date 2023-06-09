from Bio import SeqIO
import sys
import argparse

parser = argparse.ArgumentParser(
    description="Remove duplicate sequences and record labels in separate file."
)
parser.add_argument("infile", type=argparse.FileType('r'), 
                    help="<input FASTA> File containing sequences to compress.")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output FASTA> File containing unique sequences.")
parser.add_argument("csvfile", type=argparse.FileType('w'),
                    help="<output CSV> File to write duplicate labels.")
args = parser.parse_args()

records = SeqIO.parse(args.infile, 'fasta')

unique = {}
for i, record in enumerate(records):
    seq = str(record.seq)
    if seq.count('X') / len(seq) > 0.1:
        print(f"Discarding sequence {record.name} with {seq.count('X')} Xs")
        continue
    seq = seq.replace('X', '')
    if seq not in unique:
        unique.update({seq: []})
    unique[seq].append(record.description)

print(f"Reduced {i} to {len(unique)} unique sequences")

args.csvfile.write("label,duplicate\n")

for idx, item in enumerate(unique.items()):
    seq, labels = item
    args.outfile.write(f">{labels[0]}_{len(labels)}\n{seq}\n")
    for label in labels[1:]:
        args.csvfile.write(f'"{labels[0]}","{label}"\n')

args.csvfile.close()
args.outfile.close()


