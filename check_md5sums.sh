#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Verifies MD5 checksums of rehydrated FASTA genome files against the expected
# values listed in md5sum.txt. Only .fna files (read from their .gz archives)
# are checked; other file types are silently skipped. Reports each mismatch
# and prints a total count at the end.
#
# USAGE:
#   ./check_md5sums.sh [--gendir DIR]
#
# ARGUMENTS:
#   --gendir DIR   Directory containing the .fna.gz genome files and the
#                  md5sum.txt checksum file.
#                  Default: genomes
#
# OUTPUT:
#   Prints a "Wrong checksum: ..." line for each mismatch found, followed
#   by a summary with the total number of failed checks.

set -euo pipefail

# The target directory with rehydrated genomes.
GEN_DIR="genomes"

#-------------------------------------------------------------------------------
parse_args() {
    # Parses command-line arguments and overrides the corresponding global variables.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --gendir)
                GEN_DIR="$2"
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
echo "Checking files..."

parse_args "$@"

# Check for the existence of the file with MD5 checksums.
if [[ ! -f "$GEN_DIR/md5sum.txt" ]]; then
    echo "The directory \"$GEN_DIR\" does not contain the \"md5sum.txt\" file."
    exit 1
fi

# Count .fna entries in the checksums file and initialise the file and error counters.
sumcount=$(grep -c '_genomic\.fna' "$GEN_DIR/md5sum.txt" || true)
counter=1
incorr=0
while read line; do
    # Get the file name (part of the md5sum line after the last space),
    # and the file extension.
    file="$GEN_DIR/${line##*\ }"
    ext=${file##*.}
    
    # Skip non-FASTA files.
    if [[ $ext != "fna" ]]; then
        continue
    fi

    echo -ne "\rChecking ${counter}/${sumcount} ${file}\033[K"

    # Get the checksum (the part before the first space),
    # generate it for the downloaded file and strip the
    # file name off it.
    s1=${line%%\ *}                     
    s2=$(gunzip -c "${file}.gz" | md5sum)
    s2=${s2%%\ *}

    # Check whether the checksums are equal. If not, report it.
    if [[ $s1 != $s2 ]]; then
        echo -e "\nWrong checksum: $file $s1 $s2"
        incorr=$(( incorr + 1 ))
    fi

    # Increase the checksum counter.
    counter=$(( counter + 1 ))
done < "$GEN_DIR/md5sum.txt"

echo -e "\nThe check has been finished. The number of wrong MD5 checksums: ${incorr}."
