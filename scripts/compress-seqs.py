from Bio import SeqIO
import sys

if len(sys.argv) != 4:
    print("Usage: python3 unique-seqs.py [input FASTA] [output FASTA] [output CSV]")
    sys.exit()

handle = open(sys.argv[1])
records = SeqIO.parse(handle, 'fasta')

unique = {}
for i, record in enumerate(records):
    seq = str(record.seq)
    if seq not in unique:
        unique.update({seq: []})
    unique[seq].append(record.description)

print(f"Reduced {i} to {len(unique)} unique sequences")

outfasta = open(sys.argv[2], 'w')
outcsv = open(sys.argv[3], 'w')
outcsv.write("label,duplicate\n")

for idx, item in enumerate(unique.items()):
    seq, labels = item
    outfasta.write(f">{labels[0]}_{len(labels)}\n{seq}\n")
    for label in labels[1:]:
        outcsv.write(f'"{labels[0]}","{label}"\n')

outcsv.close()
outfasta.close()


