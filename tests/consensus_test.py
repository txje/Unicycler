
from __future__ import print_function
from __future__ import division

import sys
sys.dont_write_bytecode = True
sys.path.append('/Users/Ryan/Programs/Hybrid_assembler')
from semi_global_aligner import load_references, load_long_reads


references = load_references('msa_test_full_span.fasta', 0)
full_sequences = [x.sequence for x in references]
full_sequence_qualities = []

references = load_references('msa_test_start_only.fasta', 0)
start_sequences = [x.sequence for x in references]
start_sequence_qualities = []

references = load_references('msa_test_end_only.fasta', 0)
end_sequences = [x.sequence for x in references]
end_sequence_qualities = []

from semi_global_aligner import AlignmentScoringScheme
scoring_scheme = AlignmentScoringScheme('3,-6,-5,-2')

from cpp_function_wrappers import multiple_sequence_alignment
print(multiple_sequence_alignment(full_sequences, full_sequence_qualities,
                                  start_sequences, start_sequence_qualities,
                                  end_sequences, end_sequence_qualities,
                                  scoring_scheme, 1000))
