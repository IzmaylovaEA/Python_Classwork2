

import argparse
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Finding the name of the closest organism')
    parser.add_argument('-i', '--input', help='Input FASTA file', metavar='FILE', required=True)
    parser.add_argument('-t', '--threads', help='Number of threads', metavar='Int', type=int, default=1)
    parser.add_argument('-o', '--output', help='The resulting filename', metavar='FILE', required=True)
    args = parser.parse_args()
    sys.stdout = open(args.output, 'w')
    for record in SeqIO.parse(args.input, 'fasta'):
        seq = str(record.seq)
        result = NCBIWWW.qblast("blastn", "nt", seq)
        records = NCBIXML.read(result)
        min_e_value = 0.04
        seq_name = ''
        for alignment in records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < min_e_value:
                    seq_name = alignment.title.split('|')[4]
                    min_e_value = hsp.expect
        record.id = ''
        record.description = seq_name
        SeqIO.write(record, sys.stdout, 'fasta')
    sys.stdout.close()
