from __future__ import print_function
from __future__ import division
import sys
import os
import subprocess


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
