import argparse
import logging
import coloredlogs
import os
import sys

from parallel_augustus import checks
from parallel_augustus import pipeline


def main():
    parser = argparse.ArgumentParser(
        prog="parallel_augustus",
        description="\n\nBreaks the input genome into chunks to give to Augustus in parallel",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    mandatory_args = parser.add_argument_group("Mandatory arguments")
    mandatory_args.add_argument(
        "--genome",
        "-g",
        action="store",
        dest="input_genome",
        help="Input genome file break into chunks",
        default=None,
        required=True,
    )
    mandatory_args.add_argument(
        "--output",
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory name",
        default="hapog_results",
        required=False,
    )
    mandatory_args.add_argument(
        "--chunks",
        "-c",
        action="store",
        dest="chunks",
        help="Number of chunks to divide the genome into.",
        default=8,
        type=int,
        required=False,
    )
    mandatory_args.add_argument(
        "--processes",
        "-p",
        action="store",
        dest="processes",
        help="Number of Augustus processes to launch in parallel",
        default=8,
        type=int,
        required=True,
    )

    mandatory_args.add_argument(
        "--extra",
        "-e",
        action="store",
        dest="augustus_params",
        type=str,
        nargs=argparse.REMAINDER,
        help="Parameters to pass to augustus between quotes and separated by spaces. Example: --extra '--species=human --protein=on'",
    )

    args = parser.parse_args()

    coloredlogs.install(
        level="DEBUG",
        fmt="%(asctime)s - %(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args.input_genome = os.path.abspath(args.input_genome)
    checks.run_checks(args.input_genome)
    pipeline.run(
        args.input_genome,
        args.output_dir,
        args.chunks,
        args.processes,
        args.augustus_params,
    )
