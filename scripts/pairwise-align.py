import subprocess
import tempfile
import argparse
from seqUtils import convert_fasta
from gotoh2 import Aligner

#from mpi4py import MPI
#comm = MPI.COMM_WORLD
#my_rank = comm.Get_rank()
#nprocs = comm.Get_size()

aligner = Aligner()  # default settings


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref', type=argparse.FileType('r'),
                        help="FASTA file with reference sequence to cut out of queries")
    parser.add_argument('fasta', type=argparse.FileType('r'),
                        help="FASTA file containing query sequences")
    parser.add_argument('outfile', type=str, default=None,
                        help="file prefix to write output")
    parser.add_argument('--cutoff', type=float, default=0.8, 
                        help="")
    return parser.parse_args()


def mafft(query, ref):
    """ Pairwise alignment, removing insertions """
    handle = tempfile.NamedTemporaryFile(delete=False)
    s = '>ref\n{}\n>query\n{}\n'.format(ref, query)
    handle.write(s.encode('utf-8'))
    handle.close()

    # call MAFFT on temporary file
    stdout = subprocess.check_output(['mafft', '--quiet', handle.name])
    stdout = stdout.decode('utf-8')
    result = convert_fasta(stdout.split('\n'))
    aligned_ref = result[0][1]
    aligned_query = result[1][1]

    # exclude any insertions relative to the reference
    trimmed_query = ""
    for i, nt in enumerate(aligned_ref):
        if nt == '-':
            continue
        trimmed_query += aligned_query[i]

    return(trimmed_query)


def gotoh(query, ref):
    aref, aquery, ascore = aligner.align(ref, query)
    trimmed_query = ''.join([nt for i, nt in enumerate(aquery) if aref[i]!='-'])
    return trimmed_query, ascore/len(aref)


args = parse_args()
#outfile = open(f"{args.outfile}.{my_rank}.fa", 'w')
outfile = open(args.outfile, 'w')

# read reference sequence from file
refseq = convert_fasta(args.ref)[0][1]
fasta = convert_fasta(args.fasta)

for count, record in enumerate(fasta):
    #if count % nprocs != my_rank:
    #    continue
    if count > 100:
        break
    h, s = record
    query = s.replace('-', '')

    #trimmed = mafft(query, refseq)
    trimmed, ascore = gotoh(query, refseq)
    print(ascore)
    if ascore > args.cutoff:
        outfile.write('>{}\n{}\n'.format(h, trimmed))
    else:
        print(f"Discarding sequence {h} with alignment score {ascore}")

    if count % 10 == 0:
        print(f"{count}/{len(fasta)}") #({my_rank}/{nprocs})")

outfile.close()
