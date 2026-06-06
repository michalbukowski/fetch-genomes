#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Restructures the rehydrated NCBI dataset directory into a flat directory of
# genome FASTA files. For each assembly under <gendir>/ncbi_dataset/data/,
# it creates a hard link to the assembly's .fna.gz file in the output
# directory, naming each link <ACCESSION>_genomic.fna.gz. Exits with an error
# if the output directory already exists or if an assembly directory contains
# an unexpected number of .fna.gz files.
#
# USAGE:
#   ./restruct.sh [--gendir DIR] [--flatdir DIR]
#
# ARGUMENTS:
#   --gendir  DIR   Rehydrated NCBI dataset directory to restructure.
#                   Default: genomes
#   --flatdir DIR   Output directory for hard-linked .fna.gz genome files.
#                   Must not exist yet; created by the script.
#                   Default: fna
#
# OUTPUT:
#   <flatdir>/   Flat directory of hard-linked .fna.gz genome files,
#                one per assembly accession.

set -euo pipefail

# The source directory (rehydrated NCBI dataset directory) and
# the target one (flat output directory for genome FASTA files).
GEN_DIR="genomes"
FLAT_DIR="fna"

#-------------------------------------------------------------------------------
parse_args() {
    # Parses command-line arguments and overrides the corresponding global variables.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --gendir)
                GEN_DIR="$2"
                shift 2
                ;;
            --flatdir)
                FLAT_DIR="$2"
                shift 2
                ;;
            *)
                echo "Unknown argument: $1"
                exit 1
                ;;
        esac
    done
}

#-------------------------------------------------------------------------------
parse_args "$@"

if [[ -e "${FLAT_DIR}" ]]; then
    echo "The directory \"${FLAT_DIR}\" exists!"
    exit 1
fi

mkdir -p "${FLAT_DIR}"

asm_dirs=("${GEN_DIR}/ncbi_dataset/data/"*/)   # Collect all per-assembly subdirectories.
count="${#asm_dirs[@]}"
for i in "${!asm_dirs[@]}"; do
    asm_dir="${asm_dirs[$i]}"
    if [[ ! -d "$asm_dir" ]]; then
        continue   # Skip non-directory entries (e.g., metadata files).
    fi

    asm=$(basename "$asm_dir")   # Assembly accession (equals the directory name).
    echo -ne "\rProcessing $(( i + 1 ))/${count}: ${asm}\033[K"

    files=("${asm_dir}/${asm}_"*"_genomic.fna.gz")
    if [[ "${#files[@]}" -ne 1 ]]; then
        echo -e "\nWrong file number: ${#files[@]}"
        exit 1
    fi
    ln "${files[0]}" "${FLAT_DIR}/${asm}_genomic.fna.gz"   # Hard-link into the flat output directory.
done

echo -e "\nDone"
