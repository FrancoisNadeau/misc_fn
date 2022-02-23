#!/usr/bin/python3

from argparse import ArgumentParser
from docstring_parser import parse as dsparse

def get_desc(function_name: str) -> tuple:
    """
    Parse a function's docstring automatically
    
    Args:
        function_name: str
            Name of the function for which to get documentation from.

    Returns: tuple(desc, help_msgs)
        desc: str
            Concatenation of short and long function descriptions
        help_msgs: tuple(str)
            Tuple of strings representing each parameter's help message
        
    """

    parsed = dsparse(function_name.__doc__)
    help_msgs = tuple(prm.description for prm
                      in parsed.params)
    desc = '\n'.join([parsed.short_description,
                      parsed.long_description])
    return desc, help_msgs

def main():
    get_desc(sys.argv[1])

if __name__ == '__main__':
    main()

