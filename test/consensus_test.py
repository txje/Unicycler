#!/usr/bin/env python
"""
This script is a test of the multiple sequence alignment in Seqan and the corresponding consensus
sequence calling.
You must first run consensus_test_generate_sequences.py to prepare the files.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lib'))
from read_ref import load_references, load_long_reads
from alignment import AlignmentScoringScheme
from cpp_function_wrappers import multiple_sequence_alignment


read_dict, read_names = load_long_reads('msa_test_full_span.fastq', 0)
full_span_sequences = [read_dict[x].sequence for x in read_names]
full_span_qualities = [read_dict[x].qualities for x in read_names]

read_dict, read_names = load_long_reads('msa_test_start_only.fastq', 0)
start_only_sequences = [read_dict[x].sequence for x in read_names]
start_only_qualities = [read_dict[x].qualities for x in read_names]

read_dict, read_names = load_long_reads('msa_test_end_only.fastq', 0)
end_only_sequences = [read_dict[x].sequence for x in read_names]
end_only_qualities = [read_dict[x].qualities for x in read_names]

scoring_scheme = AlignmentScoringScheme('3,-6,-5,-2')

print(multiple_sequence_alignment(full_span_sequences, full_span_qualities,
                                  start_only_sequences, start_only_qualities,
                                  end_only_sequences, end_only_qualities,
                                  scoring_scheme, 1000))
