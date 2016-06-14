'''
This script makes use of several C++ functions which are in cpp_functions.so. They are wrapped in
similarly named Python functions.
'''
from __future__ import print_function
from __future__ import division
import os
from ctypes import CDLL, cast, c_char_p, c_int, c_double, c_void_p, c_bool, POINTER

C_LIB = CDLL(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'cpp_functions.so'))



# This is the big semi-global C++ Seqan alignment function at the heart of the aligner.
C_LIB.semiGlobalAlignment.argtypes = [c_char_p, # Read name
                                      c_char_p, # Read sequence
                                      c_int,    # Verbosity
                                      c_double, # Expected slope
                                      c_void_p, # KmerPositions pointer
                                      c_int,    # Match score
                                      c_int,    # Mismatch score
                                      c_int,    # Gap open score
                                      c_int,    # Gap extension score
                                      c_double, # Low score threshold
                                      c_bool,   # Return bad alignments
                                      c_int]    # K-mer size
C_LIB.semiGlobalAlignment.restype = c_void_p    # String describing alignments

def semi_global_alignment(read_name, read_sequence, verbosity, expected_slop, kmer_positions_ptr,
                          match_score, mismatch_score, gap_open_score, gap_extend_score,
                          low_score_threshold, keep_bad, kmer_size):
    '''
    Python wrapper for semiGlobalAlignment C++ function.
    '''
    ptr = C_LIB.semiGlobalAlignment(read_name, read_sequence, verbosity, expected_slop,
                                    kmer_positions_ptr, match_score, mismatch_score,
                                    gap_open_score, gap_extend_score, low_score_threshold,
                                    keep_bad, kmer_size)
    return c_string_to_python_string(ptr)



# These functions are used to conduct short alignments for the sake of extending the start and end
# of a GraphMap alignment.
C_LIB.startExtensionAlignment.argtypes = [c_char_p, # Read sequence
                                          c_char_p, # Reference sequence
                                          c_int,    # Match score
                                          c_int,    # Mismatch score
                                          c_int,    # Gap open score
                                          c_int]    # Gap extension score
C_LIB.startExtensionAlignment.restype = c_void_p    # String describing alignment

def start_extension_alignment(realigned_read_seq, realigned_ref_seq, scoring_scheme):
    '''
    Python wrapper for startExtensionAlignment C++ function.
    '''
    ptr = C_LIB.startExtensionAlignment(realigned_read_seq, realigned_ref_seq,
                                        scoring_scheme.match, scoring_scheme.mismatch,
                                        scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)

C_LIB.endExtensionAlignment.argtypes = [c_char_p, # Read sequence
                                        c_char_p, # Reference sequence
                                        c_int,    # Match score
                                        c_int,    # Mismatch score
                                        c_int,    # Gap open score
                                        c_int]    # Gap extension score
C_LIB.endExtensionAlignment.restype = c_void_p    # String describing alignment

def end_extension_alignment(realigned_read_seq, realigned_ref_seq, scoring_scheme):
    '''
    Python wrapper for endExtensionAlignment C++ function.
    '''
    ptr = C_LIB.endExtensionAlignment(realigned_read_seq, realigned_ref_seq,
                                      scoring_scheme.match, scoring_scheme.mismatch,
                                      scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)



# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None

def c_string_to_python_string(c_string):
    '''
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    '''
    python_string = cast(c_string, c_char_p).value
    C_LIB.freeCString(c_string)
    return python_string



# These functions make/delete a C++ object that will be used during line-finding.
C_LIB.newKmerPositions.argtypes = []
C_LIB.newKmerPositions.restype = c_void_p

def new_kmer_positions():
    '''
    Python wrapper for newKmerPositions C++ function.
    '''
    return C_LIB.newKmerPositions()

C_LIB.addKmerPositions.argtypes = [c_void_p, # KmerPositions pointer
                                   c_char_p, # Name
                                   c_char_p, # Sequence
                                   c_int]    # K-mer size
C_LIB.addKmerPositions.restype = None

def add_kmer_positions(kmer_positions_ptr, name, sequence, kmer_size):
    '''
    Python wrapper for addKmerPositions C++ function.
    '''
    C_LIB.addKmerPositions(kmer_positions_ptr, name, sequence, kmer_size)

C_LIB.deleteAllKmerPositions.argtypes = [c_void_p]
C_LIB.deleteAllKmerPositions.restype = None

def delete_all_kmer_positions(kmer_positions_ptr):
    '''
    Python wrapper for deleteAllKmerPositions C++ function.
    '''
    C_LIB.deleteAllKmerPositions(kmer_positions_ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentScores.argtypes = [c_int, # Random sequence length
                                                   c_int, # Count
                                                   c_int, # Match score
                                                   c_int, # Mismatch score
                                                   c_int, # Gap open score
                                                   c_int] # Gap extension score
C_LIB.getRandomSequenceAlignmentScores.restype = c_void_p

def get_random_sequence_alignment_mean_and_std_dev(seq_length, count, scoring_scheme):
    '''
    Python wrapper for getRandomSequenceAlignmentScores C++ function.
    '''
    ptr = C_LIB.getRandomSequenceAlignmentScores(seq_length, count,
                                                 scoring_scheme.match, scoring_scheme.mismatch,
                                                 scoring_scheme.gap_open, scoring_scheme.gap_extend)
    return_str = c_string_to_python_string(ptr)
    return_parts = return_str.split(',')
    return float(return_parts[0]), float(return_parts[1])



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.getRandomSequenceAlignmentErrorRates.argtypes = [c_int, # Random sequence length
                                                       c_int, # Count
                                                       c_int, # Match score
                                                       c_int, # Mismatch score
                                                       c_int, # Gap open score
                                                       c_int] # Gap extension score
C_LIB.getRandomSequenceAlignmentErrorRates.restype = c_void_p

def get_random_sequence_alignment_error_rates(seq_length, count, scoring_scheme):
    '''
    Python wrapper for getRandomSequenceAlignmentErrorRate C++ function.
    '''
    ptr = C_LIB.getRandomSequenceAlignmentErrorRates(seq_length, count,
                                                     scoring_scheme.match, scoring_scheme.mismatch,
                                                     scoring_scheme.gap_open,
                                                     scoring_scheme.gap_extend)
    return c_string_to_python_string(ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.simulateDepths.argtypes = [POINTER(c_int), # Alignment lengths
                                 c_int,          # Alignment count
                                 c_int,          # Reference length
                                 c_int,          # Iterations
                                 c_int]          # Threads
C_LIB.simulateDepths.restype = c_void_p

def simulate_depths(read_lengths, ref_length, iterations, threads):
    '''
    Python wrapper for simulateDepths C++ function.
    '''
    read_lengths_array = (c_int * len(read_lengths))(*read_lengths)
    ptr = C_LIB.simulateDepths(read_lengths_array, len(read_lengths), ref_length, iterations,
                               threads)
    return c_string_to_python_string(ptr)



# This function gets the mean and standard deviation of alignments between random sequences.
C_LIB.multipleSequenceAlignment.argtypes = [POINTER(c_char_p), # Sequences
                                            POINTER(c_char_p), # Qualities
                                            c_int,             # Sequence count
                                            c_int,             # Piece size
                                            c_int,             # Match score
                                            c_int,             # Mismatch score
                                            c_int,             # Gap open score
                                            c_int]             # Gap extension score             
C_LIB.multipleSequenceAlignment.restype = c_void_p

def multiple_sequence_alignment(full_length_sequences, full_length_qualities,
                                start_only_sequences, start_only_qualities,
                                end_only_sequences, end_only_qualities,
                                scoring_scheme, piece_size=1000):
    '''
    Python wrapper for multipleSequenceAlignment C++ function.
    piece_size is used for chopping up long sequences which would take too long to align in their
    entirety. For example, if the sequences were 5000 bp and piece_size is 1000 bp, the C++
    function will conduct the alignment in multiple overlapping 1000 bp pieces.
    '''
    full_count = len(full_length_sequences)
    start_count = len(start_only_sequences)
    end_count = len(end_only_sequences)

    full_length_qualities = fill_out_qualities(full_length_sequences, full_length_qualities)
    start_only_qualities = fill_out_qualities(start_only_sequences, start_only_qualities)
    end_only_qualities = fill_out_qualities(end_only_sequences, end_only_qualities)

    # Pad out start/end sequences/qualities with N/+.
    mean_full_length = int(round(sum([len(x) for x in full_length_sequences]) / \
                                 len(full_length_sequences)))

    for i, start_only_sequence in enumerate(start_only_sequences):
        missing_bases = mean_full_length - len(start_only_sequence)
        full_length_sequences.append(start_only_sequence + ('N' * missing_bases))
        full_length_qualities.append(start_only_qualities[i] + ('+' * missing_bases))
    for i, end_only_sequence in enumerate(end_only_sequences):
        missing_bases = mean_full_length - len(end_only_sequence)
        full_length_sequences.append(('N' * missing_bases) + end_only_sequence)
        full_length_qualities.append(('+' * missing_bases) + end_only_qualities[i])

    sequences = (c_char_p * len(full_length_sequences))(*full_length_sequences)
    qualities = (c_char_p * len(full_length_qualities))(*full_length_qualities)
    ptr = C_LIB.multipleSequenceAlignment(sequences, qualities, len(full_length_sequences),
                                          piece_size,
                                          scoring_scheme.match, scoring_scheme.mismatch,
                                          scoring_scheme.gap_open, scoring_scheme.gap_extend)
    result = c_string_to_python_string(ptr)
    result_parts = result.split(';')
    consensus = result_parts[0]
    all_scores = result_parts[1].split(',')
    full_length_scores = all_scores[:full_count]
    start_only_scores = all_scores[full_count:full_count + start_count]
    end_only_scores = all_scores[-end_count:]
    return consensus, full_length_scores, start_only_scores, end_only_scores


def fill_out_qualities(sequences, qualities):
    '''
    Given a list of sequences and qualities, this function fills out any missing qualities with '+'
    (PHRED+33 for 10% chance of error).
    '''
    while len(qualities) < len(sequences):
        qualities.append('')
    for i, seq in enumerate(sequences):
        qualities[i] += '+' * (len(seq) - len(qualities[i]))
        qualities[i] = qualities[i][:len(seq)]
    return qualities

