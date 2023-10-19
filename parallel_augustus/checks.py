import logging
import os
import shutil


def run_checks(genome: str):
    if not os.path.exists(genome):
        logging.error(f"{genome} not found")
        exit(1)

    check_dependencies()


def check_dependencies():
    missing_dependency = False
    tools = ["augustus"]
    for tool in tools:
        if not is_in_path(tool):
            logging.error(f"{tool} not found")
            missing_dependency = True
        else:
            logging.debug(f"Found {tool}")

    if missing_dependency:
        exit(-1)


def is_in_path(tool: str):
    return shutil.which(tool) is not None
