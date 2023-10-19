from Bio import SeqIO
import glob
import logging
import numpy as np
import os
import subprocess
import time
import warnings


def run(genome: str, output_dir: str, chunks: int, processes: int, params: str):
    create_directories(output_dir)
    genome_size = get_genome_size(genome)
    create_chunks(genome, genome_size, chunks)
    launch_augustus(processes, params)
    concatenate_results()


def create_directories(output_dir: str):
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
        logging.error(
            f"Unsufficient permissions to write output directory {output_dir}"
        )
        exit(1)

    os.chdir(output_dir)

    try:
        os.mkdir("chunks")
        os.mkdir("augustus")
        os.mkdir("logs")
    except:
        logging.error("Could not create subdirectories")


def get_genome_size(genome: str):
    logging.info(f"Parsing {genome}")
    cumul_size = 0
    with open(genome) as genome_file:
        for record in SeqIO.parse(genome_file, "fasta"):
            cumul_size += len(record.seq)
    logging.debug(f"Parsed {cumul_size} bases")
    return cumul_size


def create_chunks(genome: str, genome_size: int, chunks: int):
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
                current_chunk_file = open(
                    f"chunks/chunks_{current_chunk + 1}.fasta", "w"
                )
                current_chunk += 1
                chunk_sizes.append(current_chunk_size)
                current_chunk_size = 0

            current_chunk_file.write(record.format("fasta"))
            current_chunk_size += len(record.seq)

        chunk_sizes.append(current_chunk_size)
        current_chunk_file.close()

    median_chunk_size = np.median(chunk_sizes)
    nb_chunks = len(chunk_sizes)
    logging.info(
        f"Created {nb_chunks} chunks with a median size of {int(median_chunk_size)} bases"
    )


def launch_augustus(processes: int, params: str):
    logging.info("Launching Augustus on each chunk")

    augustus_params = []
    if params:
        for p in params:
            param = p.split(" ")
            for pa in param:
                augustus_params.append(pa)

    procs = []
    for chunk in glob.glob("chunks/*.fasta"):
        chunk_prefix = chunk.split("/")[-1].replace(".fasta", "")

        cmd = ["augustus"]
        cmd.append(chunk)
        for p in augustus_params:
            cmd.append(p)

        with open("logs/augustus.cmds", "a") as cmd_file:
            print(" ".join(cmd), flush=True, file=cmd_file)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            procs.append(
                subprocess.Popen(
                    cmd,
                    stdout=open(f"augustus/{chunk_prefix}.gff", "w"),
                    stderr=open(f"logs/augustus_{chunk_prefix}.e", "w"),
                )
            )

        # Only launch a job if there is less than 'processes' running
        # Otherwise, wait for any to finish before launching a new one
        while len([p for p in procs if p.poll() is None]) >= int(processes):
            time.sleep(10)

    has_failed = False
    for p in procs:
        p.wait()

        return_code = p.returncode
        if return_code != 0:
            logging.error(
                f"ERROR: Augustus didn't finish successfully, exit code: {return_code}"
            )
            logging.error("Faulty command: %s" % (" ".join(p.args)))
            has_failed = True

    if has_failed:
        exit(1)

    logging.info("Augustus finished successfully")


def concatenate_results():
    out = open("augustus.gff", "w")

    for chunk in glob.glob("augustus/*.gff"):
        with open(chunk) as inf:
            for line in inf:
                if not line.startswith("#"):
                    line = line.rstrip("\n")
                    out.write(f"{line}\n")

    out.close()
    logging.info("Done. Results can be found in augustus.gff.")
