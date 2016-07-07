"""
Functions relating to BLAST.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
from .misc import quit_with_error, float_to_str


class CannotFindStart(Exception):
    pass


def find_start_gene(sequence, start_genes_fasta, identity_threshold, coverage_threshold, out_dir,
                    makeblastdb_path, tblastn_path, threads, verbosity):
    """
    This function uses tblastn to look for start genes in the sequence. It returns the first gene
    (using the order in the file) which meets the identity and coverage thresholds, as well as
    the position of that gene (including which strand it is on).
    """
    # Create a FASTA file of the replicon sequence.
    blast_dir = os.path.join(out_dir, 'blast_temp')
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    replicon_fasta_filename = os.path.join(blast_dir, 'replicon.fasta')
    replicon_fasta = open(replicon_fasta_filename, 'w')
    replicon_fasta.write('>replicon\n')
    replicon_fasta.write(sequence)
    replicon_fasta.write('\n')
    replicon_fasta.close()

    # Build the BLAST database.
    command = [makeblastdb_path, '-dbtype', 'nucl', '-in', replicon_fasta_filename]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = process.communicate()
    if err:
        quit_with_error('makeblastdb encountered an error:\n' + err.decode())

    # Run the tblastn search.
    command = [tblastn_path, '-db', replicon_fasta_filename, '-query', start_genes_fasta, '-outfmt',
               '6 qseqid sstart send pident qlen qseq qstart', '-num_threads', str(threads)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
        line = process.stdout.readline().strip().decode()
        if line != '':
            parts = line.split('\t')
            qseqid = parts[0]
            sstart = int(parts[1]) - 1
            send = int(parts[2])
            pident = float(parts[3])
            qlen = float(parts[4])
            qseq = parts[5]
            qstart = int(parts[6])
            query_cov = 100.0 * len(qseq) / qlen

            if verbosity > 2:
                hit_str = qseqid + ' ' + float_to_str(query_cov, 2) + '% cov' + \
                    float_to_str(pident, 2) + '% identity'

            good_hit = pident >= identity_threshold and query_cov >= coverage_threshold and \
                qstart == 1

            if not good_hit and verbosity > 2:
                print('  Insufficient BLAST hit:', hit_str, flush=True)

            if good_hit:
                if sstart <= send:
                    start_pos = sstart - 1
                    flip = False
                else:
                    start_pos = sstart
                    flip = True
                process.terminate()
                if verbosity > 2:
                    print('  Successful BLAST hit:', hit_str, flush=True)
                return qseqid, start_pos, flip

    raise CannotFindStart
