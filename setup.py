from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md")) as f:
    long_description = f.read()

setup(
    name="parallel_augustus",
    packages=["parallel_augustus"],
    version="1.0.1",
    license="CeCILL",
    description="Simple wrapper around Augustus to bring faster restitution times",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=["Benjamin Istace"],
    url="https://github.com/bistace/parallel_augustus",
    download_url="https://github.com/bistace/parallel_augustus",
    keywords=[
        "bioinformatics",
        "genomics",
        "genome",
        "annotation",
        "augustus"
    ],
    install_requires=["biopython", "coloredlogs"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        "Programming Language :: Python :: 3.8",
    ],
    entry_points={"console_scripts" : ["parallel_augustus = parallel_augustus.cli:main"]},
)