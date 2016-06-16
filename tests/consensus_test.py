#!/usr/bin/env python
'''
This script is a test of the multiple sequence alignment in Seqan and the corresponding consensus
sequence calling.
You must first run consensus_test_generate_sequences.py to prepare the files.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division

import os
import sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lib'))
from read_ref import load_references, load_long_reads
from alignment import AlignmentScoringScheme
from cpp_function_wrappers import multiple_sequence_alignment


references = load_references('msa_test_full_span.fasta', 0)
full_sequences = [x.sequence for x in references]
full_sequence_qualities = []

references = load_references('msa_test_start_only.fasta', 0)
start_sequences = [x.sequence for x in references]
start_sequence_qualities = []

references = load_references('msa_test_end_only.fasta', 0)
end_sequences = [x.sequence for x in references]
end_sequence_qualities = []

scoring_scheme = AlignmentScoringScheme('3,-6,-5,-2')

print(multiple_sequence_alignment(full_sequences, full_sequence_qualities,
                                  start_sequences, start_sequence_qualities,
                                  end_sequences, end_sequence_qualities,
                                  scoring_scheme, 1000))
