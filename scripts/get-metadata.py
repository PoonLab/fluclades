import csv
import time
import argparse
from Bio import Phylo, Entrez, SeqIO

parser = argparse.ArgumentParser(
    description="Retrieve metadata from Genbank based on accession numbers "
                "in tip labels of input tree.  Assumes tip labels are "
                "underscore delimited, with accession number in first "
                "position.")

parser.add_argument("infile", type=argparse.FileType('r'),
                    help="<input NWK> File containing Newick tree string")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="<output CSV> File to write metadata.")
parser.add_argument("--email", type=str, default="apoon42@uwo.ca",
                    help="<option> Required for Entrez database transactions.")
parser.add_argument("--batchsize", type=int, default=100, 
                    help="<option> Number of records to retrieve in one batch.")

args = parser.parse_args()

# let Genbank know who we are
Entrez.email = args.email

phy = Phylo.read(args.infile, 'newick')
tips = phy.get_terminals()
print(f"Processing {len(tips)} tips in tree...")

writer = csv.writer(args.outfile)
writer.writerow(['accn', 'strain', 'serotype', 'host', 'country', 'coldate'])

for i in range(0, len(tips), args.batchsize):
    print(i)
    batch = [tip.name.split("_")[0] for tip in tips[i:(i+args.batchsize)]]
    handle = Entrez.efetch(db='nuccore', id=batch, rettype='gb', 
                           retmode='text', retmax=args.batchsize)
    records = SeqIO.parse(handle, 'gb')
    for record in records:
        source = list(filter(lambda f: f.type=='source', record.features))
        quals = source[0].qualifiers
        writer.writerow([
            record.name,  # LOCUS
            quals.get('strain', [''])[0],
            quals.get('serotype', [''])[0],
            quals.get('host', [''])[0],
            quals.get('country', [''])[0],
            quals.get('collection_date', [''])[0]
        ])
    handle.close()
    time.sleep(1)  # wait to avoid spamming server

args.outfile.close()
