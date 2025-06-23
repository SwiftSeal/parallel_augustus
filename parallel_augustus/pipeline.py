from Bio import SeqIO
import glob
import logging
import os
import subprocess
import time
import warnings


def run(genome: str, output_dir: str, chunks: int, processes: int, params: str):
    create_directories(output_dir)
    split_fasta_entries(genome)
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
        os.mkdir("entries")
        os.mkdir("augustus")
        os.mkdir("logs")
    except Exception:
        logging.error("Could not create subdirectories")


def split_fasta_entries(genome: str):
    logging.info(f"Splitting {genome} into single-entry FASTA files")
    entry_count = 0
    with open(genome) as genome_file:
        for i, record in enumerate(SeqIO.parse(genome_file, "fasta"), 1):
            entry_path = f"entries/entry_{i}.fasta"
            with open(entry_path, "w") as entry_file:
                SeqIO.write(record, entry_file, "fasta")
            entry_count += 1
    logging.info(f"Created {entry_count} single-entry FASTA files in 'entries/'")


def launch_augustus(processes: int, params: str):
    logging.info("Launching Augustus on each entry file")

    augustus_params = []
    if params:
        for p in params:
            param = p.split(" ")
            for pa in param:
                augustus_params.append(pa)

    procs = []
    for entry in glob.glob("entries/*.fasta"):
        entry_prefix = entry.split("/")[-1].replace(".fasta", "")

        cmd = ["augustus"]
        cmd.append(entry)
        for p in augustus_params:
            cmd.append(p)

        with open("logs/augustus.cmds", "a") as cmd_file:
            print(" ".join(cmd), flush=True, file=cmd_file)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            procs.append(
                subprocess.Popen(
                    cmd,
                    stdout=open(f"augustus/{entry_prefix}.gff", "w"),
                    stderr=open(f"logs/augustus_{entry_prefix}.e", "w"),
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

    for entry in glob.glob("augustus/*.gff"):
        with open(entry) as inf:
            # Skip the header
            line = ""
            while line != "#\n":
                line = inf.readline()

            for line in inf:
                line = line.rstrip("\n")
                out.write(f"{line}\n")

    out.close()
    logging.info("Done. Results can be found in augustus.gff.")
