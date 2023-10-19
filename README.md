# Efficiently run Augustus in multi-core environments
`parallel_augustus` cuts the input genome in chunks to feed to Augustus. The number of chunks and Augustus processes launched in parallel are configurable.

# Installation
```
pip install --user parallel_augustus
```

# Requirements
  - Python >= 3.8
  - Augustus

# Execution
Before launching `parallel_augustus`, please make sure that Augustus is available in your path. If you run `augustus -h` and do not encounter an error, you are good to go.

`parallel_augustus` first divides the genome (`-g` argument) in the desired number of chunks (`-c` argument, it can be any number > 1). After that it launches `-p` processes of Augustus in parallel until there is no chunks left. You can pass parameters to Augustus via the `--extra` flag. At least a `--species=thing` is required by Augustus. 

_IMPORTANT_: the `--extra` argument should be the last one in the command line due to shortcomings in the `argparse` python module.

With all that in mind, a typical `parallel_augustus` command line will look like this:
```bash
parallel_augustus -g genome.fasta \
    -o augustus_results \
    -c 500 \
    -p 20 \
    --extra '--species=human'
```
This command creates the output directory `augustus_results`, divides the genome into 500 chunks and launches 20 processes of Augustus in parallel. At the end, results are concatenated into `augustus_results/augustus.gff`.