#!/usr/bin/env python3
"""
This script isn't part of Unicycler but is run separately to process UniProt gene sequences.

It operates on UniProt clustered sequences like those you can get here:
http://www.uniprot.org/uniref/?query=uniprot:(gene%3Adnaa+length%3A%5B10+TO+1000%5D+%28taxonomy%3Aarchaea+OR+taxonomy%3Abacteria%29)+identity:0.9
http://www.uniprot.org/uniref/?query=uniprot:(gene%3Arepa+length%3A%5B10+TO+1000%5D+%28taxonomy%3Aarchaea+OR+taxonomy%3Abacteria%29)+identity:0.9
http://www.uniprot.org/uniref/?query=uniprot:(gene%3Aters+taxonomy%3Aviruses)+identity:0.9
http://www.uniprot.org/uniref/?query=uniprot:(gene%3Aterl+taxonomy%3Aviruses)+identity:0.9

Author: Ryan Wick
email: rrwick@gmail.com
"""

import gzip
import argparse
from misc import quit_with_error, get_sequence_file_type, get_compression_type, check_file_exists


def main():
    args = get_arguments()
    check_file_exists(args.gene_fasta)
    genes = load_genes(args.gene_fasta)
    genes = sorted(genes, key=lambda x: x[2], reverse=True)
    for gene in genes:
        if gene[2] >= args.min_n:
            print('>' + gene[0])
            print(gene[1])


def get_arguments():
    parser = argparse.ArgumentParser(description='UniProt gene sequence processor')

    parser.add_argument('gene_fasta', type=str,
                        help='FASTA file (optionally gzipped) with cluster representative gene '
                             'sequences from UniProt')
    parser.add_argument('--min_n', type=int, default=2,
                        help='Clusters with fewer than this number of sequences will not be '
                             'included')
    return parser.parse_args()


def load_genes(fasta_filename):
    """
    This function loads in sequences from a FASTA file and returns a list of Reference objects.
    """
    try:
        if get_sequence_file_type(fasta_filename) != 'FASTA':
            quit_with_error(fasta_filename + ' is not in FASTA format')
    except ValueError:
        quit_with_error(fasta_filename + ' is not in FASTA format')

    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    genes = []
    fasta_file = open_func(fasta_filename, 'rt')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.decode("utf-8").strip()
        if not line:
            continue
        if line.startswith('>'):  # Header line = start of new contig
            if name:
                genes.append((name, sequence))
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    fasta_file.close()
    if name:
        genes.append((name, sequence))

    genes_with_count = []
    for gene in genes:
        count = 0
        name_parts = gene[0].split()
        for part in name_parts:
            if part.startswith('n='):
                count = int(part[2:])
        genes_with_count.append((gene[0], gene[1], count))

    return genes_with_count


if __name__ == '__main__':
    main()
