#! /usr/bin/env python

import argparse as ap
import sys

from src.extract.extract import extract_parser
from src.model.model     import model_parser
from src.version import __version__

def main():
    """
    Main function to run the inSTRbility tool.
    """

    parser = ap.ArgumentParser(prog='inSTRbility',
                               add_help=False,
                               formatter_class=ap.RawTextHelpFormatter,
                               description='inSTRbility: A tool to estimate instability of STR loci from long read sequencing data')

    parser._action_groups.pop()
    print("inSTRbility - Analysing somatic instability at tandem repeat loci.\nDashnow Lab\n", file=sys.stderr)

    parser.add_argument('-h', '--help',    action='help',    help="Print help and exit")
    parser.add_argument('-v', '--version', action='version', help="Print version", version=f'inSTRbility {__version__}')
    
    subparsers = parser.add_subparsers(dest='command')

    model_parser(subparsers)
    extract_parser(subparsers)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()

    if args.command is None:
        print("Usage:")
        print("    inSTRbility [OPTIONS] <COMMAND>\n")
        print("Commands:")
        for name, sp in subparsers.choices.items():
            print(f"  {name:<9} {sp.description}")
        print("\nOptions:")
        print("  -h, --help     Print help")
        print("  -v, --version  Print version")
        sys.exit()


if __name__ == "__main__":
    main()
