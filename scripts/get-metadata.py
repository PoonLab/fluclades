import csv
import time
from Bio import Phylo, Entrez, SeqIO
Entrez.email = 'apoon42@uwo.ca'

metadata = {}

phy = Phylo.read("data/removeX.ft2.nwk", 'newick')
tips = phy.get_terminals()

writer = csv.writer(open("data/genbank-metadata.csv", 'w'))
writer.writerow(['accn', 'strain', 'serotype', 'host', 'country', 'coldate'])

for i in range(0, len(tips), 100):
    print(i)
    batch = [tip.name.split("_")[0] for tip in tips[i:(i+100)]]
    handle = Entrez.efetch(db='nuccore', id=batch, rettype='gb', 
                           retmode='text', retmax=100)
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
    time.sleep(1)
