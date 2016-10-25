"""
Miscellaneous function used by various scripts.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import sys
import os
import subprocess
import random
import math
import gzip
import argparse
import shutil
import re


def float_to_str(num, decimals, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if decimals == 0:
        return int_to_str(int(round(num)), max_num=max_num)
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
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))


def check_files_and_programs(files, args, spades_path=None, graphmap_path=None,
                             makeblastdb_path=None, tblastn_path=None, gene_db_path=None,
                             pilon_path=None, java_path=None, samtools_path=None,
                             bowtie2_path=None, bowtie2_build_path=None):
    """
    Checks to make sure all files in the list are present and either program, as needed.
    """
    for file in files:
        check_file_exists(file)
    if graphmap_path:
        check_graphmap(graphmap_path)
    if spades_path:
        check_spades(spades_path)
    if makeblastdb_path and tblastn_path:
        check_blast(makeblastdb_path, tblastn_path, gene_db_path)
    if pilon_path:
        check_pilon(pilon_path, java_path, samtools_path, bowtie2_path, bowtie2_build_path, args)


def check_file_exists(filename):  # type: (str) -> bool
    """
    Checks to make sure the single given file exists.
    """
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)


def quit_with_error(message):  # type: (str) -> None
    """
    Displays the given message and ends the program's execution.
    """
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def check_graphmap(graphmap_path):
    """
    Makes sure the GraphMap executable is available.
    """
    if shutil.which(graphmap_path) is None:
        quit_with_error('could not find GraphMap at ' + graphmap_path +
                        ', either fix path or run with --no_graphmap')


def check_spades(spades_path):
    """
    Makes sure the SPAdes executable is available.
    """
    if shutil.which(spades_path) is None:
        quit_with_error('could not find SPAdes at ' + spades_path)

    command = [spades_path, '-h']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()

    if not err.decode():
        quit_with_error('SPAdes was found but does not produce output (make sure to use '
                        '"spades.py" location, not "spades")')


def check_blast(makeblastdb_path, tblastn_path, gene_db_path):
    """
    Makes sure the BLAST executables are available.
    """
    if shutil.which(makeblastdb_path) is None:
        quit_with_error('could not find makeblastdb - either specify its location using '
                        '--makeblastdb_path or use --no_rotate to remove BLAST dependency')
    if shutil.which(tblastn_path) is None:
        quit_with_error('could not find tblastn - either specify its location using '
                        '--tblastn_path or use --no_rotate to remove BLAST dependency')
    if not os.path.isfile(gene_db_path):
        quit_with_error('could not find file: ' + gene_db_path +
                        '\neither specify a different start gene database using --start_genes '
                        'or use --no_rotate')


def check_pilon(pilon_path, java_path, samtools_path, bowtie2_path, bowtie2_build_path, args):
    """
    Makes sure the Pilon executable is available. Unlike the other tools, Pilon's target name may
    change based on the version (e.g. pilon-1.18.jar, pilon-1.19.jar, etc.). This function will
    therefore set args.pilon_path to the first matching jar file it finds.
    """
    if shutil.which(java_path) is None:
        quit_with_error('could not find java - either specify its location using '
                        '--java_path or use --no_pilon to remove Java dependency')
    if shutil.which(samtools_path) is None:
        quit_with_error('could not find samtools - either specify its location using '
                        '--samtools_path or use --no_pilon to remove Samtools dependency')
    if shutil.which(bowtie2_path) is None:
        quit_with_error('could not find bowtie2 - either specify its location using '
                        '--bowtie2_path or use --no_pilon to remove Bowtie2 dependency')
    if shutil.which(bowtie2_build_path) is None:
        quit_with_error('could not find bowtie2-build - either specify its location using '
                        '--bowtie2_build_path or use --no_pilon to remove Bowtie2 dependency')
    found_pilon_path = get_pilon_jar_path(pilon_path)
    if found_pilon_path:
        args.pilon_path = found_pilon_path
    else:
        quit_with_error('could not find pilon*.jar - either specify its location using '
                        '--pilon_path or use --no_pilon to remove Pilon dependency')


def get_pilon_jar_path(pilon_path):
    """
    Returns the path to pilon.jar. If the given path is correct, it just returns that, as an
    absolute path. Otherwise it tries to find it.
    """
    if pilon_path and os.path.isfile(pilon_path):
        return os.path.abspath(pilon_path)
    for directory in os.environ['PATH'].split(':'):
        try:
            path_files = [f for f in os.listdir(directory)
                          if os.path.isfile(os.path.join(directory, f))]
        except FileNotFoundError:
            path_files = []
        pilon_jars = [f for f in path_files if f.startswith('pilon') and f.endswith('.jar')]
        if pilon_jars:
            return os.path.join(directory, sorted(pilon_jars)[-1])  # return the latest version
    return None


def get_mean_and_st_dev(num_list):
    """
    This function returns the mean and standard deviation of the given list of numbers.
    """
    num = len(num_list)
    if num == 0:
        return None, None
    mean = sum(num_list) / num
    if num == 1:
        return mean, None
    sum_squares = sum((x - mean) ** 2 for x in num_list)
    st_dev = (sum_squares / (num - 1)) ** 0.5
    return mean, st_dev


def print_progress_line(completed, total, base_pairs=None, prefix=None, end_newline=False):
    """
    Prints a progress line to the screen using a carriage return to overwrite the previous progress
    line.
    """
    progress_str = ''
    if prefix:
        progress_str += prefix
    progress_str += int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'
    if base_pairs is not None:
        progress_str += ' - ' + int_to_str(base_pairs) + ' bp'
    end_char = '\n' if end_newline else ''
    print('\r' + progress_str, end=end_char, flush=True)


def get_nice_header(header):
    """
    For a header with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    """
    if is_header_spades_format(header):
        return 'NODE_' + header.split('_')[1]
    else:
        return header.split()[0]


def is_header_spades_format(contig_name):
    """
    Returns whether or not the header appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    """
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and \
        (contig_name_parts[0] == 'NODE' or contig_name_parts[0] == 'EDGE') and \
        contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'


def reverse_complement(seq):
    """
    Given a DNA sequences, this function returns the reverse complement sequence.
    """
    return ''.join([complement_base(seq[i]) for i in range(len(seq) - 1, -1, -1)])


def complement_base(base):
    """
    Given a DNA base, this returns the complement.
    """
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
    """
    Returns a random base with 25% probability of each.
    """
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
    """
    Returns a random sequence of the given length.
    """
    sequence = ''
    for _ in range(length):
        sequence += get_random_base()
    return sequence


def get_median(sorted_list):
    """
    Returns the median of a list of numbers. Assumes the list has already been sorted.
    """
    count = len(sorted_list)
    index = (count - 1) // 2
    if count % 2:
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1]) / 2.0


def get_percentile(unsorted_list, percentile):
    """
    Returns a percentile of a list of numbers. Doesn't assume the list has already been sorted.
    Implements the nearest rank method:
    https://en.wikipedia.org/wiki/Percentile#The_Nearest_Rank_method
    """
    return get_percentile_sorted(sorted(unsorted_list), percentile)


def get_percentile_sorted(sorted_list, percentile):
    """
    Same as the above function, but assumes the list is already sorted.
    """
    if not sorted_list:
        return 0.0
    fraction = percentile / 100.0
    rank = int(math.ceil(fraction * len(sorted_list)))
    if rank == 0:
        return sorted_list[0]
    return sorted_list[rank - 1]


def weighted_average(num_1, num_2, weight_1, weight_2):
    """
    A simple weighted mean of two numbers.
    """
    weight_sum = weight_1 + weight_2
    return num_1 * (weight_1 / weight_sum) + num_2 * (weight_2 / weight_sum)


def weighted_average_list(nums, weights):
    """
    A simple weighted mean of a list of numbers.
    """
    w_sum = sum(weights)
    if w_sum == 0.0:
        return 0.0
    else:
        return sum(num * (weights[i] / w_sum) for i, num in enumerate(nums))


def print_section_header(message, verbosity, last_newline=True):
    """
    Prints a header for stdout, unless verbosity is zero, in which case it does nothing.
    """
    if verbosity > 0:
        print('\n')
        print(bold_yellow_underline(message), end=('\n' if last_newline else ''), flush=True)


def round_to_nearest_odd(num):
    """
    Rounds a float to an odd integer.
    """
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
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
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


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    seq_file = open_func(filename, 'rt')
    first_char = seq_file.read(1)

    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is neither FASTA or FASTQ')


def get_num_agreement(num_1, num_2):
    """
    Returns a value between 0.0 and 1.0 describing how well the numbers agree.
    1.0 is perfect agreement and 0.0 is the worst.
    """
    if num_1 == 0.0 and num_2 == 0.0:
        return 1.0
    if num_1 < 0.0 and num_2 < 0.0:
        num_1 *= -1
        num_2 *= -1
    if num_1 * num_2 < 0.0:
        return 0.0
    return min(num_1, num_2) / max(num_1, num_2)


def flip_number_order(num_1, num_2):
    """
    Given two segment numbers, this function possibly flips them around. It returns the new numbers
    (either unchanged or flipped) and whether or not a flip took place. The decision is somewhat
    arbitrary, but it needs to be consistent so when we collect bridging read sequences they are
    always in the same direction.
    """
    if num_1 > 0 and num_2 > 0:
        flip = False
    elif num_1 < 0 and num_2 < 0:
        flip = True
    elif num_1 < 0:  # only num_1 is negative
        flip = abs(num_1) > abs(num_2)
    else:  # only num_2 is negative
        flip = abs(num_2) > abs(num_1)
    if flip:
        return (-num_2, -num_1), True
    else:
        return (num_1, num_2), False


def load_fasta(filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta file.
    """
    fasta_seqs = []
    fasta_file = open(filename, 'rt')
    name = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if not line:
            continue
        if line[0] == '>':  # Header line = start of new contig
            if name:
                fasta_seqs.append((name.split()[0], sequence))
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence))
    fasta_file.close()
    return fasta_seqs


def score_function(val, half_score_val):
    """
    For inputs of 0.0 and greater, this function returns a value between 0.0 and 1.0, approaching
    1.0 with large values. The half_score_val argument is the point at which the function returns
    0.5. If it's large the function approaches 1.0 more slowly, if it's small the function
    approaches 1.0 more quickly.
    """
    return 1.0 - (half_score_val / (half_score_val + val))


def strip_read_extensions(read_file_name):
    """
    This function removes extensions from a file name.
    """
    base_name = os.path.basename(read_file_name)
    name_parts = base_name.split('.')
    for i in range(2):
        if len(name_parts) > 1 and len(name_parts[-1]) <= 5:
            name_parts = name_parts[:-1]
    return '.'.join(name_parts)


def add_line_breaks_to_sequence(sequence, line_length):
    """
    Wraps sequences to the defined length.  All resulting sequences end in a line break.
    """
    seq_with_breaks = ''
    pos = 0
    while pos < len(sequence):
        seq_with_breaks += sequence[pos:pos+line_length] + '\n'
        pos += line_length
    return seq_with_breaks


class MyHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _split_lines(self, text, width):
        if text.startswith('B|') or text.startswith('R|'):
            text_lines = text[2:].splitlines()
            wrapped_text_lines = []
            for line in text_lines:
                if len(line) <= width:
                    wrapped_text_lines.append(line)
                else:
                    wrap_column = 2

                    # The bridging mode help text should wrap each line around to the column of
                    # the equals sign.
                    if text.startswith('B|'):
                        line_parts = line.split()
                        wrap_column += line.find('=')
                        join = ''
                        current_line = '  ' + line_parts[0]

                    # The other multi-option help texts should wrap an entire option at a time.
                    else:  # text.startswith('R|')
                        line_parts = line.split(', ')
                        join = ','
                        current_line = line_parts[0]
                    for part in line_parts[1:]:
                        if len(current_line) + len(join) + 1 + len(part) <= width:
                            current_line += join + ' ' + part
                        else:
                            wrapped_text_lines.append(current_line + join)
                            current_line = ' ' * wrap_column + part
                    wrapped_text_lines.append(current_line)
            return wrapped_text_lines
        else:
            return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        if text.startswith('R|'):
            return argparse.RawDescriptionHelpFormatter._fill_text(self, text[2:], width, indent)
        else:
            return argparse.HelpFormatter._fill_text(self, text, width, indent)


def print_table(table, alignments='', max_col_width=30, col_separation=3, indent=2,
                green_row=-1):
    """
    Args:
        table: a list of lists of strings (one row is one list, all rows must be the same length)
        alignments: a string of L and R, indicating the alignment for each row
        max_col_width: values longer than this will be wrapped
        col_separation: the number of spaces between columns
        indent: the number of spaces between the table and the left side of the terminal
        green_row: the index of a row in the table to make green

    Returns:
        nothing, just prints the table
    """
    alignments += 'L' * (len(table[0]) - len(alignments))  # Fill out with L, if incomplete
    col_widths = [0] * len(table[0])
    for row in table:
        col_widths = [min(max(col_widths[i], len_without_format(x)), max_col_width)
                      for i, x in enumerate(row)]
    separator = ' ' * col_separation
    indenter = ' ' * indent
    for i, row in enumerate(table):
        while any(x for x in row):
            aligned_row = []
            for value, col_width, alignment in zip(row, col_widths, alignments):
                if i == 0 or alignment == 'L':
                    aligned_row.append(value.ljust(col_width)[:col_width])
                else:
                    aligned_row.append(value.rjust(col_width)[:col_width])
            row_str = separator.join(aligned_row)
            if i == 0:
                row_str = bold_underline(row_str)
            if i == green_row:
                row_str = green(row_str)
            else:
                row_str = row_str.replace('PASS', green('PASS'))
                row_str = row_str.replace('FAIL', red('FAIL'))
            print(indenter + row_str, flush=True)
            row = [x[col_widths[j]:] for j, x in enumerate(row)]


def colour(text, text_colour):
    if text_colour.lower() == 'green':
        return green(text)
    elif text_colour.lower() == 'red':
        return red(text)
    else:
        return text


def green(text):
    return '\033[32m' + text + '\033[0m'


def bold_green(text):
    return '\033[1m' + '\033[32m' + text + '\033[0m'


def red(text):
    return '\033[31m' + text + '\033[0m'


def bold(text):
    return '\033[1m' + text + '\033[0m'


def dim(text):
    return '\033[2m' + text + '\033[0m'


def bold_underline(text):
    return '\033[1m' + '\033[4m' + text + '\033[0m'


def bold_yellow_underline(text):
    return '\033[1m' + '\033[93m' + '\033[4m' + text + '\033[0m'


def len_without_format(text):
    return len(re.sub('\033.*?m', '', text))


def get_all_files_in_current_dir():
    return [f for f in os.listdir('.') if os.path.isfile(f)]


def convert_fastq_to_fasta(fastq, fasta):
    if get_compression_type(fastq) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(fastq, 'rt') as fastq:
        with open(fasta, 'wt') as fasta:
            for line in fastq:
                name = line.strip()[1:].split()[0]
                sequence = next(fastq).strip()
                _ = next(fastq)
                _ = next(fastq)
                fasta.write('>' + name + '\n')
                fasta.write(sequence + '\n')


def print_verbosity(text, verbosity, min_verbosity):
    if verbosity >= min_verbosity:
        print(text, flush=True)
