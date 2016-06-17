#!/usr/bin/env python
'''
This script prepares random sequences for the consensus_test.py script.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division

import random
import os
import sys
sys.dont_write_bytecode = True
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'lib'))
from misc import get_random_sequence, get_random_base

def main():
    length = int(sys.argv[1])
    mutation_rate = float(sys.argv[2])
    full_span_count = int(sys.argv[3])
    start_only_count = int(sys.argv[4])
    end_only_count = int(sys.argv[5])

    base_filename = 'msa_test'

    original_sequence = get_random_sequence(length)
    base_file = open(base_filename + '_original.fastq', 'w')
    base_file.write('>original\n')
    base_file.write(original_sequence + '\n')
    base_file.close()

    full_span_file = open(base_filename + '_full_span.fastq', 'w')
    for i in range(full_span_count):
        mutated_sequence = mutate_sequence(original_sequence, length * mutation_rate)
        full_span_file.write('@full_span_' + str(i+1) + '\n')
        full_span_file.write(mutated_sequence + '\n')
        full_span_file.write('+\n')
        full_span_file.write(get_random_qualities(len(mutated_sequence)) + '\n')
    full_span_file.close()

    start_only_file = open(base_filename + '_start_only.fastq', 'w')
    for i in range(start_only_count):
        mutated_sequence = mutate_sequence(original_sequence, length * mutation_rate)
        seq_length = random.randint(200, length)
        mutated_sequence = mutated_sequence[:seq_length]
        start_only_file.write('>start_only_' + str(i+1) + '\n')
        start_only_file.write(mutated_sequence + '\n')
        start_only_file.write('+\n')
        start_only_file.write(get_random_qualities(len(mutated_sequence)) + '\n')
    start_only_file.close()

    end_only_file = open(base_filename + '_end_only.fastq', 'w')
    for i in range(end_only_count):
        mutated_sequence = mutate_sequence(original_sequence, length * mutation_rate)
        seq_length = random.randint(200, length)
        mutated_sequence = mutated_sequence[-seq_length:]
        end_only_file.write('>end_only_' + str(i+1) + '\n')
        end_only_file.write(mutated_sequence + '\n')
        end_only_file.write('+\n')
        end_only_file.write(get_random_qualities(len(mutated_sequence)) + '\n')
    end_only_file.close()

def mutate_sequence(sequence, mutation_count):
    for _ in range(int(round(mutation_count))):
        pos = random.randint(0, len(sequence) - 1)
        mut_type = random.randint(0, 2)
        if mut_type == 0:
            sequence = sequence[:pos] + get_random_base() + sequence[pos+1:]
        elif mut_type == 1:
            sequence = sequence[:pos] + sequence[pos+1:]
        elif mut_type == 2:
            sequence = sequence[:pos] + get_random_base() + sequence[pos:]
    return sequence

def get_random_qualities(length):
    qualities = ''
    for _ in range(length):
        qualities += get_random_quality()
    return qualities

def get_random_quality():
    qual = random.randint(1, 20)
    return str(unichr(qual+33))

if __name__ == '__main__':
    main()
