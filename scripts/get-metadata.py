import csv
import time
import argparse
from Bio import Phylo, Entrez, SeqIO
import sys
import re

pat = re.compile("[A-Z]+_*[0-9]+\.[1-9]")

parser = argparse.ArgumentParser(
    description="Retrieve metadata from Genbank based on accession numbers "
                "in tip labels of input tree.  Assumes tip labels are "
                "underscore delimited, with accession number in first "
                "position.")

parser.add_argument("infile", type=argparse.FileType('r'),
                    help="<input> File in Newick or FASTA format")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output CSV> File to write metadata.")
parser.add_argument("--email", type=str, default="apoon42@uwo.ca",
                    help="<option> Required for Entrez database transactions.")
parser.add_argument("--batchsize", type=int, default=100, 
                    help="<option> Number of records to retrieve in one batch.")

args = parser.parse_args()

# let Genbank know who we are
Entrez.email = args.email

ifn = args.infile.name
if ifn.endswith(".nwk"):
    phy = Phylo.read(args.infile, 'newick')
    labels = [tip.name for tip in phy.get_terminals()]
elif ifn.endswith(".mafft") or ifn.endswith(".fa") or ifn.endswith(".fasta"):
    records = SeqIO.parse(args.infile, 'fasta')
    labels = [rec.description for rec in records]
else:
    print("Error: input file must have .nwk/.fa/.fasta extension")
    sys.exit()

print(f"Processing {len(labels)} labels from input...")

writer = csv.writer(args.outfile)
writer.writerow(['accn', 'strain', 'serotype', 'host', 'country', 'coldate'])

for i in range(0, len(labels), args.batchsize):
    print(i)

    batch = []
    for j in range(i, min(i+args.batchsize, len(labels))):
        matches = pat.findall(labels[j])
        if len(matches) == 0:
            sys.stderr.write(f"WARNING: Failed to parse {labels[j]}, skipping\n")
            continue
        batch.append(matches[0])
    
    #batch = [lab.split("_")[0] for lab in labels[i:(i+args.batchsize)]]
    handle = Entrez.efetch(db='nuccore', id=batch, rettype='gb', 
                           retmode='text', retmax=args.batchsize)
    
    records = SeqIO.parse(handle, 'gb')
    for record in records:
        source = list(filter(lambda f: f.type=='source', record.features))
        quals = source[0].qualifiers
        writer.writerow([
            record.id,  # record.name is LOCUS, not always accession
            quals.get('strain', [''])[0],
            quals.get('serotype', [''])[0],
            quals.get('host', [''])[0],
            quals.get('country', [''])[0],
            quals.get('collection_date', [''])[0]
        ])
    handle.close()
    time.sleep(1)  # wait to avoid spamming server

args.outfile.close()
