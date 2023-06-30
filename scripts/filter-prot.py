from Bio import SeqIO
import re
import sys

if len(sys.argv) != 3:
    print("Usage: python3 filter-prot.py [infile] [outfile]")
    sys.exit()

#pat = re.compile("\[protein=[Nn]euraminidase.*\]|\[gene=[Nn][Aa]\]")
#pat = re.compile("\[[a-z]+=[a-z ]*[Pp][Bb]2[a-z ]*\]|\[protein=polymerase basic protein 2\]")
#pat = re.compile("\[[a-z]+=[a-z ]*[Pp][Bb]1[a-z ]*\]|\[protein=polymerase basic protein 1\]")
#pat = re.compile("\[[a-z]+=[a-z ]*[Pp][Aa][a-z ]*\]|\[protein=polymerase acid[ic]* protein\]|\[protein=polymerase protein A\]")
pat = re.compile("\[protein=[Nn]ucleo[cp].+\]|\[protein=[Nn][Pp]\]|\[gene=[Nn][Pp]\]")

#infile = open("data/gb-prot-na.fa")
infile = open(sys.argv[1])
#outfile = open("data/gb-prot-na.filtered.fa", 'w')
outfile = open(sys.argv[2], 'w')

records = SeqIO.parse(infile, 'fasta')
for record in records:
    matches = pat.findall(record.description)
    if len(matches) == 0:
        print(record.description)
        continue

    if "PA-X" in record.description:
        continue

    SeqIO.write(record, outfile, 'fasta')

outfile.close()
