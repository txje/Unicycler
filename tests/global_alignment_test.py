
from __future__ import print_function
from __future__ import division

import sys
sys.dont_write_bytecode = True
sys.path.append('/Users/Ryan/Programs/Hybrid_assembler')
from semi_global_aligner import load_references, load_long_reads

band_size = int(sys.argv[1])

sequence_1 = load_references('global_alignment_test_1.fasta', 0)
sequence_1 = sequence_1[0].sequence

sequence_2 = load_references('global_alignment_test_2.fasta', 0)
sequence_2 = sequence_2[0].sequence

from semi_global_aligner import AlignmentScoringScheme
scoring_scheme = AlignmentScoringScheme('3,-6,-5,-2')

from cpp_function_wrappers import fully_global_alignment
print(fully_global_alignment(sequence_1, sequence_2, scoring_scheme, True, band_size))
