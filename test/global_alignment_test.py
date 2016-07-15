#!/usr/bin/env python
"""
This script is a test of the banded global alignment in Seqan.
You must first run global_alignment_test_generate_sequences.py to prepare the files.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lib'))
from alignment import AlignmentScoringScheme
from cpp_function_wrappers import fully_global_alignment
from read_ref import load_references

band_size = int(sys.argv[1])

sequence_1 = load_references('global_alignment_test_1.fasta', 0)
sequence_1 = sequence_1[0].sequence

sequence_2 = load_references('global_alignment_test_2.fasta', 0)
sequence_2 = sequence_2[0].sequence

scoring_scheme = AlignmentScoringScheme('3,-6,-5,-2')

print(fully_global_alignment(sequence_1, sequence_2, scoring_scheme, True, band_size))
