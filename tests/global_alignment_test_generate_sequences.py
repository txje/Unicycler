
from __future__ import print_function
from __future__ import division

import random
import sys

def main():
    length_1 = int(sys.argv[1])
    length_2 = int(sys.argv[2])
    mutation_rate = float(sys.argv[3])

    sequence_1 = get_random_sequence(length_1)
    sequence_1_file = open('global_alignment_test_1.fasta', 'w')
    sequence_1_file.write('>sequence_1\n')
    sequence_1_file.write(sequence_1 + '\n')
    sequence_1_file.close()

    length_diff = length_1 - length_2

    mutation_count = int(round(mutation_rate * length_1))
    substitution_count = int(round(mutation_count / 3))
    indel_count = int(round(mutation_count * 2 / 3))

    if length_diff > indel_count:
        indel_count = length_diff

    insertion_count = int(round((indel_count - length_diff) / 2))
    deletion_count = length_diff + insertion_count

    if deletion_count < 0:
        insertion_count -= deletion_count
        deletion_count = 0

    print()
    print('target mutation count:', mutation_count)
    print('length diff:          ', length_diff)
    print()
    print('substitution count:   ', substitution_count)
    print('target indel count:   ', indel_count)
    print()
    print('insertion count:      ', insertion_count)
    print('deletion count:       ', deletion_count)

    sequence_2 = sequence_1
    sequence_2 = mutate_sequence(sequence_2, deletion_count, 'del')
    sequence_2 = mutate_sequence(sequence_2, substitution_count, 'sub')
    sequence_2 = mutate_sequence(sequence_2, insertion_count, 'ins')

    sequence_2_file = open('global_alignment_test_2.fasta', 'w')
    sequence_2_file.write('>sequence_2\n')
    sequence_2_file.write(sequence_2 + '\n')
    sequence_2_file.close()


def get_random_base():
    rand_int = random.randint(0, 3)
    if rand_int == 0:
        return 'A'
    elif rand_int == 1:
        return 'C'
    elif rand_int == 2:
        return 'G'
    elif rand_int == 3:
        return 'T'

def get_random_sequence(length):
    sequence = ''
    for _ in range(length):
        sequence += get_random_base()
    return sequence

def mutate_sequence(sequence, mutation_count, mut_type):
    for _ in range(mutation_count):
        pos = random.randint(0, len(sequence) - 1)
        if mut_type == 'sub':
            sequence = sequence[:pos] + get_random_base() + sequence[pos+1:]
        elif mut_type == 'del':
            sequence = sequence[:pos] + sequence[pos+1:]
        elif mut_type == 'ins':
            sequence = sequence[:pos] + get_random_base() + sequence[pos:]
    return sequence

if __name__ == '__main__':
    main()
