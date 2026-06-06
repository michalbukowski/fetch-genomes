## fetch-genomes

A set of Bash scripts that download all known genomic sequences from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/) for a given list of taxonomy IDs, using [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/) tools. Assembly metadata is saved as TSV files and genome sequences are stored as gzip-compressed FASTA files.

### 1. Prerequisites

The scripts rely on [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/) tools, which are installed and managed through a Conda environment. To set up the environment you need to have Conda installed, for instance [Miniconda](https://docs.conda.io/en/latest/miniconda.html). The scripts were tested on Ubuntu 24.04.4 LTS x86-64.

### 2. Setting up the environment

Create and activate the `ncbi-cli` Conda environment:
```bash
conda env create --file envs/ncbi-cli.yml
conda activate ncbi-cli
```

### 3. Input files

Before running the scripts, prepare the two input files located by default in the `input/` directory:
- `taxids.txt` &ndash; one [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) ID per line, identifying the taxa for which genome assemblies will be downloaded
- `fields.txt` &ndash; metadata column names to include in the output TSV files, one per line; the first name must be `accession` &ndash; a comprehensive list of all available field names is provided in the reference for [dataformat](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/tsv/dataformat_tsv_genome) tool

In the `input/` directory you will find example files for downloading genomes of _Staphylococcus agnetis_ (taxid 985762) and _Staphylococcus delphini_ (taxid 53344). You can test the scripts using them before changing anything &ndash; simply run them in the proper order with default argument values (see the following sections).

### 4. Directory structure and script files

The scripts use the following directory structure:
```
fetch-genomes/
├── envs/
│   └── ncbi-cli.yml
├── input/
│   ├── taxids.txt
│   └── fields.txt
├── output/
├── genomes/
├── fna/
├── genomes.zip
├── fetch_genomes.sh
├── check_md5sums.sh
└── restruct.sh
```
`fetch_genomes.sh`, `check_md5sums.sh`, and `restruct.sh` are the main scripts described in section 5. The `envs/` directory contains the Conda environment file for NCBI Datasets CLI. Input files with taxonomy IDs and metadata field names are located in `input/`. The `output/` directory is created by `fetch_genomes.sh` and contains assembly metadata TSV files and the list of accession numbers. `genomes.zip` is the intermediate dehydrated genome archive downloaded by `fetch_genomes.sh`. The `genomes/` directory is also created by `fetch_genomes.sh` and holds the rehydrated NCBI dataset with the downloaded genome sequences. You can use the `check_md5sums.sh` script for an optional but recommended verification of MD5 checksums of the downloaded genome files. The `fna/` directory is created by `restruct.sh` and provides a flat view of all genome FASTA files, each hard-linked from `genomes/`.

### 5. Scripts

| Script | Description |
|---|---|
| `fetch_genomes.sh` | Downloads assembly metadata from NCBI for each taxid listed in `input/taxids.txt` and saves it as per-taxid TSV files (`assemblies_taxid-*.tsv`) and one merged TSV file (`assemblies.tsv`) in the `output/` directory. Then downloads and rehydrates the corresponding genome sequences into the `genomes/` directory. |
| `check_md5sums.sh` | Verifies MD5 checksums of the rehydrated `.fna.gz` genome files against the expected values listed in `md5sum.txt` in the `genomes/` directory. Reports each mismatch and prints the total count of failed checks at the end. |
| `restruct.sh` | Restructures the NCBI dataset directory into a flat directory of genome FASTA files by creating hard links for all `.fna.gz` files from the `genomes/ncbi_dataset/data/` directory into the `fna/` directory, one file per assembly accession. |

### 6. Running the scripts

With the `ncbi-cli` environment active and the input files prepared, run the scripts in order from the repository directory.

Download assembly metadata and genome sequences:
```bash
./fetch_genomes.sh
```
Optional arguments:
- `--taxid FILE` &ndash; path to the file with taxids (default: `input/taxids.txt`)
- `--fields FILE` &ndash; path to the file with metadata field names (default: `input/fields.txt`)
- `--outdir DIR` &ndash; output directory for TSV files and the accession list (default: `output/`)
- `--genarch FILE` &ndash; path for the ZIP file with dehydrated genomes (default: `genomes.zip`)
- `--gendir DIR` &ndash; directory for the rehydrated genomes (default: `genomes/`)
- `--redo` &ndash; forces re-download of assembly metadata even when output files from a previous run already exist; by default, if `output/accessions.txt` is present, the script skips metadata download and proceeds directly to downloading genomes

Optionally, verify the integrity of the downloaded genome files:
```bash
./check_md5sums.sh
```
Optional argument:
- `--gendir DIR` &ndash; directory containing the `.fna.gz` genome files and the `md5sum.txt`  file (default: `genomes/`)

Restructure the downloaded genomes into a flat FASTA directory:
```bash
./restruct.sh
```
Optional arguments:
- `--gendir DIR` &ndash; rehydrated NCBI dataset directory to restructure (default: `genomes/`)
- `--flatdir DIR` &ndash; output directory for the hard-linked `.fna.gz` files; must not exist yet (default: `fna/`)

You do not need to pass any arguments to the scripts. If none are provided, the scripts will work correctly using the default values and the directory structure described in section 4.

Once all three scripts have been run and the genome sequences are collected in the flat `fna/` directory
(or the directory set with `--flatdir`), the intermediate archive and dataset directory can be safely removed:
```bash
rm -rf genomes.zip genomes
```
Note: `rm -rf` permanently deletes the files. Make sure the `restruct.sh` script completed successfully before removing them.