from Bio import SeqIO
import logging
import os


def run(genome: str, output_dir: str):
    create_directories(output_dir)
    genome_size = get_genome_size(genome)


def create_directories(output_dir):
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        logging.error(f"Directory {output_dir} already exists.")
        logging.error("Please remove it before launching parallel_augustus.")
        exit(1)
    except FileNotFoundError:
        logging.error(f"Path to {output_dir} does not exist.")
        exit(1)
    except PermissionError:
        logging.error(f"Unsufficient permissions to write output directory {output_dir}")
        exit(1)

    os.chdir(output_dir)

    try:
        os.mkdir("chunks")
        os.mkdir("augustus")
        os.mkdir("results")
    except:
        logging.error("Could not create subdirectories")


def get_genome_size(genome):
    logging.info(f"Parsing {genome}")
    cumul_size = 0
    with open(genome) as genome_file:
        for record in SeqIO.parse(genome_file, "fasta"):
            cumul_size += len(record.seq)
    logging.info(f"Parsed {cumul_size} bases")
    return cumul_size