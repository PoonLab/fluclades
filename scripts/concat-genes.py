import sys
import argparse
import re

"""
location is an inclusive range "i..j" of a feature in the nucleotide sequence

Genbank uses a "<i.." prefix to indicate that the feature start coordinate
precedes the left terminus of the nucleotide sequence (incomplete coverage).

A "..>j" prefix indicates the feature end coordinate comes after the right
terminus of the nucleotide sequence.

A "join(x..y)" operator is used to connect spliced intervals.
"""

pat = re.compile("^.+\[location=[a-z]*([(.,<>0-9)]+)\].+$")
nums = re.compile("[0-9]+")

def convert_fasta (handle):
    result = []
    sequence = ''
    for line in handle:
        if line.startswith('>'):
            if len(sequence) > 0:
                result.append([h,sequence])
                sequence = ''   # reset
            h = line.strip('>\n')
        else:
            sequence += line.strip('\n').upper()
    result.append([h,sequence]) # handle last entry
    return result


def concat(records):
   # extract feature locations from sequence names
    coords = []
    for h, _ in records:
        matches = pat.findall(h)
        if len(matches) == 0:
            print(block)
            sys.exit()
        coords.append([int(i) for i in nums.findall(matches[0])])
    
    # find main (unspliced) sequence first (M1)
    finds = list(filter(lambda x: len(x)==2, coords))
    if len(finds) == 0:
        print("Error: failed to locate M1 sequence")
        print(records)
        sys.exit()
    
    m1 = finds[0]
    which_m1 = coords.index(m1)
    header, seq = records[which_m1]

    # append other sequence
    for i, co in enumerate(coords):
        if i == which_m1:
            continue
        nuclen = co[-1] - (m1[-1]+1)
        aalen = nuclen//3
        aaseq = records[i][1]
        seq += aaseq[-aalen:]
        break
    
    return header, seq



parser = argparse.ArgumentParser(
    description="Concatenate alternative open reading frames and exclude overlaps"
)
parser.add_argument("infile", type=argparse.FileType('r'),
                    help="Path to file with CDS records")
parser.add_argument("outfile", type=argparse.FileType('w'),
                    help="Path to write output FASTA")
args = parser.parse_args()


blocks = args.infile.read().split("\n\n")
for block in blocks:
    try:
        records = convert_fasta(block.split('\n'))
    except UnboundLocalError:
        if block == '':
            continue  # handle extra lines or end of file
        raise
    
    # trivial case, only one feature (ORF) in sequence
    if len(records) == 1:
        args.outfile.write(block+'\n')
        continue

    header, seq = concat(records)
    args.outfile.write(f">{header}\n{seq}\n")
