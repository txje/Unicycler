"""
Functions relating to BLAST.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
from .misc import float_to_str, load_fasta


class CannotFindStart(Exception):
    pass


def find_start_gene(sequence, start_genes_fasta, identity_threshold, coverage_threshold, blast_dir,
                    makeblastdb_path, tblastn_path, threads, verbosity):
    """
    This function uses tblastn to look for start genes in the sequence. It returns the first gene
    (using the order in the file) which meets the identity and coverage thresholds, as well as
    the position of that gene (including which strand it is on).
    This function assumes that the sequence is circular with no overlap.
    """
    # Prepare the replicon sequence. In order to get a solid, single BLAST hit in cases where the
    # gene overlaps from the end to the start, we have to duplicate some of the replicon sequence
    # for the BLAST database.
    seq_len = len(sequence)
    queries = load_fasta(start_genes_fasta)
    if not queries:
        raise CannotFindStart
    longest_query = max(len(x[1]) for x in queries)
    longest_query *= 3  # amino acids to nucleotides
    dup_length = min(seq_len, longest_query)
    sequence = sequence + sequence[:dup_length]

    # Create a FASTA file of the replicon sequence.
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
        print('\nmakeblastdb encountered an error:\n' + err.decode())
        raise CannotFindStart

    # Run the tblastn search.
    command = [tblastn_path, '-db', replicon_fasta_filename, '-query', start_genes_fasta, '-outfmt',
               '6 qseqid sstart send pident qlen qseq qstart qend', '-num_threads', str(threads)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
        line = process.stdout.readline().strip().decode()
        if line != '':
            parts = line.split('\t')
            if len(parts) > 6:
                qseqid = parts[0]
                sstart = int(parts[1]) - 1
                send = int(parts[2])
                pident = float(parts[3])
                qlen = float(parts[4])
                qseq = parts[5]
                qstart = int(parts[6]) - 1
                qend = int(parts[7])
                query_cov = 100.0 * len(qseq) / qlen

                if pident >= identity_threshold and query_cov >= coverage_threshold and qstart == 0:
                    process.terminate()
                    if sstart <= send:
                        start_pos = sstart
                        flip = False
                    else:
                        start_pos = sstart + 1
                        flip = True
                    if start_pos >= seq_len:
                        start_pos -= seq_len
                    if verbosity > 2:
                        hit_str = qseqid + ', ' + float_to_str(query_cov, 2) + '% cov, ' + \
                            float_to_str(pident, 2) + '% identity, gene start pos = ' + \
                            str(qstart) + ', gene end pos = ' + str(qend) + \
                            ', replicon start pos = ' + str(sstart) + ', replicon end pos = ' + \
                            str(send)
                        print('  Successful BLAST hit:', hit_str, flush=True)
                    return qseqid, start_pos, flip

    blast_error = process.stderr.readline().strip().decode()
    if blast_error and verbosity > 1:
        print('\nBLAST encountered an error:\n' + blast_error)
    raise CannotFindStart
