from Bio import SeqIO
import logging
import numpy as np
import os
import warnings


def run(genome: str, output_dir: str, chunks: int):
    create_directories(output_dir)
    genome_size = get_genome_size(genome)
    create_chunks(genome, genome_size, chunks)


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
    logging.debug(f"Parsed {cumul_size} bases")
    return cumul_size


def create_chunks(genome, genome_size, chunks):
    logging.info(f"Trying to fragment input genome into {chunks} chunks")

    chunk_size = genome_size / chunks

    chunk_sizes = []
    current_chunk = 1
    current_chunk_size = 0
    current_chunk_file = open("chunks/chunks_1.fasta", "w")

    with open(genome) as g, warnings.catch_warnings():
        warnings.simplefilter("ignore")

        for record in SeqIO.parse(g, "fasta"):
            if current_chunk_size >= chunk_size and current_chunk != chunks:
                current_chunk_file.close()
                current_chunk_file = open(f"chunks/chunks_{current_chunk + 1}.fasta", "w")
                current_chunk += 1
                chunk_sizes.append(current_chunk_size)
                current_chunk_size = 0

            current_chunk_file.write(record.format("fasta"))
            current_chunk_size += len(record.seq)

        chunk_sizes.append(current_chunk_size)
        current_chunk_file.close()
    
    median_chunk_size = np.median(chunk_sizes)
    nb_chunks = len(chunk_sizes)
    logging.info(f"Created {nb_chunks} chunks with a median size of {int(median_chunk_size)} bases")
