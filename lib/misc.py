'''
Miscellaneous function used by various scripts.

Author: Ryan Wick
email: rrwick@gmail.com
'''

import sys
import os
import subprocess
import random
import math


def float_to_str(num, decimals, max_num=0):
    '''
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    '''
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '%.' + str(decimals) + 'f'
        num_str = num_str % num
        parts = num_str.split('.')
        before_decimal = parts[0]
        after_decimal = parts[1]
        num_str = int_to_str(int(before_decimal)) + '.' + after_decimal
    if max_num > 0:
        max_str = float_to_str(max_num, decimals)
        num_str = num_str.rjust(len(max_str))
    return num_str

def int_to_str(num, max_num=0):
    '''
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    '''
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))

def check_files_and_programs(files, spades_path=None, graphmap_path=None):
    '''
    Checks to make sure all files in the list are present and either program, as needed.
    '''
    for file in files:
        check_file_exists(file)
    if graphmap_path:
        check_graphmap(check_graphmap)
    if spades_path:
        check_spades(spades_path)

def check_file_exists(filename): # type: (str) -> bool
    '''
    Checks to make sure the single given file exists.
    '''
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)

def quit_with_error(message): # type: (str) -> None
    '''
    Displays the given message and ends the program's execution.
    '''
    print('Error:', message, file=sys.stderr)
    sys.exit(1)

def check_graphmap(graphmap_path):
    '''
    Makes sure the GraphMap executable is available.
    '''
    process = subprocess.Popen(['which', graphmap_path], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = process.communicate()
    found_graphmap = bool(out) and not bool(err)
    if not found_graphmap:
        quit_with_error('could not find GraphMap at ' + graphmap_path + 
                        ', either fix path or run with --no_graphmap')

def check_spades(spades_path):
    '''
    Makes sure the SPAdes executable is available.
    '''
    process = subprocess.Popen(['which', spades_path], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = process.communicate()
    found_spades = bool(out) and not bool(err)
    if not found_spades:
        quit_with_error('could not find SPAdes at ' + spades_path)

    command = [spades_path, '-h']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    if not err.decode('utf-8'):
        quit_with_error('SPAdes was found but does not produce output (make sure to use '
                        '"spades.py" location, not "spades")')

def get_mean_and_st_dev(num_list):
    '''
    This function returns the mean and standard deviation of the given list of numbers.
    '''
    num = len(num_list)
    if num == 0:
        return None, None
    mean = sum(num_list) / num
    if num == 1:
        return mean, None
    sum_squares = sum((x - mean) ** 2 for x in num_list)
    st_dev = (sum_squares / (num - 1)) ** 0.5
    return mean, st_dev

def print_progress_line(completed, total, base_pairs=None, prefix=None):
    '''
    Prints a progress line to the screen using a carriage return to overwrite the previous progress
    line.
    '''
    progress_str = ''
    if prefix:
        progress_str += prefix
    progress_str += int_to_str(completed) + ' / ' + int_to_str(total)
    progress_str += ' (' + '%.1f' % (100.0 * completed / total) + '%)'
    if base_pairs is not None:
        progress_str += ' - ' + int_to_str(base_pairs) + ' bp'
    print('\r' + progress_str, end='', flush=True)

def get_nice_header(header):
    '''
    For a header with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    '''
    if is_header_spades_format(header):
        return 'NODE_' + header.split('_')[1]
    else:
        return header.split()[0]

def is_header_spades_format(contig_name):
    '''
    Returns whether or not the header appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    '''
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and \
           (contig_name_parts[0] == 'NODE' or contig_name_parts[0] == 'EDGE') and \
           contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'

def reverse_complement(seq):
    '''
    Given a DNA sequences, this function returns the reverse complement sequence.
    '''
    return ''.join([complement_base(seq[i]) for i in range(len(seq) - 1, -1, -1)])

def complement_base(base):
    '''
    Given a DNA base, this returns the complement.
    '''
    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'G':
        return 'C'
    if base == 'C':
        return 'G'
    if base == 'a':
        return 't'
    if base == 't':
        return 'a'
    if base == 'g':
        return 'c'
    if base == 'c':
        return 'g'
    forward = 'RYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'YRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]

def get_random_base():
    '''
    Returns a random base with 25% probability of each.
    '''
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
    '''
    Returns a random sequence of the given length.
    '''
    sequence = ''
    for _ in range(length):
        sequence += get_random_base()
    return sequence

def get_median(sorted_list):
    '''
    Returns the median of a list of numbers. Assumes the list has already been sorted.
    '''
    count = len(sorted_list)
    index = (count - 1) // 2
    if count % 2:
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0

def get_percentile(unsorted_list, percentile):
    '''
    Returns a percentile of a list of numbers. Doesn't assume the list has already been sorted.
    Implements the nearest rank method:
    https://en.wikipedia.org/wiki/Percentile#The_Nearest_Rank_method
    '''
    sorted_list = sorted(unsorted_list)
    fraction = percentile / 100.0
    rank = int(math.ceil(fraction * len(unsorted_list)))
    if rank == 0:
        return sorted_list[0]
    return sorted_list[rank-1]

def weighted_average(num_1, num_2, weight_1, weight_2):
    '''
    A simple weighted mean of two numbers.
    '''
    weight_sum = weight_1 + weight_2
    return num_1 * (weight_1 / weight_sum) + num_2 * (weight_2 / weight_sum)

def weighted_average_list(nums, weights):
    '''
    A simple weighted mean of a list of numbers.
    '''
    w_sum = sum(weights)
    return sum(num * (weights[i] / w_sum) for i, num in enumerate(nums))

def print_section_header(message, verbosity):
    '''
    Prints a header for std out, unless verbosity is zero, in which case it does nothing.
    '''
    if verbosity > 0:
        print()
        print(message)
        print('-' * len(message), flush=True)

def round_to_nearest_odd(num):
    '''
    Rounds a float to an odd integer.
    '''
    round_up = int(math.ceil(num))
    if round_up % 2 == 0:
        round_up += 1
    round_down = int(math.floor(num))
    if round_down % 2 == 0:
        round_down -= 1
    up_diff = round_up - num
    down_diff = num - round_down
    if up_diff > down_diff:
        return round_down
    else:
        return round_up

def get_compression_type(filename):
    '''
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    '''
    magic_dict = {'gz':  (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type
