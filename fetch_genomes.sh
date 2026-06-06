#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license

# Downloads genomic sequences from the NCBI GenBank database for all known
# genome assemblies related to given taxonomy IDs (NCBI Taxonomy) using
# NCBI Datasets CLI tools. Assembly metadata is saved as TSV files, one per
# taxid and one merged file, written to the output directory. The script then
# downloads and rehydrates gzipped FASTA files with contig sequences. To fetch
# other file types (e.g., GFF3, GenBank), modify the datasets download command
# in download_dehydrated().
#
# USAGE:
#   ./fetch_genomes.sh [--taxid FILE] [--fields FILE] [--outdir DIR]
#                      [--genarch FILE] [--gendir DIR] [--redo]
#
# ARGUMENTS:
#   --taxid  FILE    Text file listing one NCBI Taxonomy ID per line.
#                    Default: input/taxids.txt
#   --fields FILE    Text file listing metadata column names, one per line.
#                    The first name must be "accession".
#                    Default: input/fields.txt
#   --outdir DIR     Directory for output TSV files and the accession list.
#                    Created automatically if it does not exist.
#                    Default: output
#   --genarch FILE   Path for the ZIP file with dehydrated genomes.
#                    Default: genomes.zip
#   --gendir DIR     Directory for the rehydrated genomes.
#                    Default: genomes
#   --redo           Force re-download of assembly metadata even when output
#                    files from a previous run already exist.
#
# OUTPUT:
#   <outdir>/assemblies_taxid-<TAXID>.tsv   Assembly metadata per taxid
#   <outdir>/assemblies.tsv                 Merged assembly metadata (all taxids)
#   <outdir>/accessions.txt                 Assembly accession numbers
#   <genarch>                               ZIP file with dehydrated genomes
#   <gendir>/                               Rehydrated gzipped FASTA files

set -euo pipefail

# Existing input files with taxid values and metadata column names.
TAXIDS_FILE="input/taxids.txt"
FIELDS_FILE="input/fields.txt"

# The directory for output files and a name for the output file that will
# contain assembly accessions to be downloaded (saved to the OUTDIR).
OUTDIR="output"
ACCNS_FILE="accessions.txt"

# Names for output files that will contain assembly metadata,
# per taxid and merged (saved to OUTDIR).
ASMS_TAXID_FILE="assemblies_taxid-%s.tsv"
ASMS_FINAL_FILE="assemblies.tsv"

# Should the assembly metadata be re-downloaded and
# the ACCNS_FILE re-created (TRUE)? If FALSE,
# the existing ACCNS_FILE will be used, if it exists.
REDO=FALSE

# Paths for the dehydrated genome archive and the rehydrated genome directory.
GEN_ARCH="genomes.zip"
GEN_DIR="genomes"

#-------------------------------------------------------------------------------
parse_args() {
    # Parses command-line arguments and overrides the corresponding global variables.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --outdir)
                OUTDIR="$2"
                shift 2
                ;;
            --taxid)
                TAXIDS_FILE="$2"
                shift 2
                ;;
            --fields)
                FIELDS_FILE="$2"
                shift 2
                ;;
            --genarch)
                GEN_ARCH="$2"
                shift 2
                ;;
            --gendir)
                GEN_DIR="$2"
                shift 2
                ;;
            --redo)
                REDO=TRUE
                shift
                ;;
            *)
                echo "Unknown argument: $1"
                exit 1
                ;;
        esac
    done
}

create_outdir() {
    # Creates OUTDIR if absent; exits if the path already exists as a non-directory.

    if [[ ! -e "$OUTDIR" ]]; then
        mkdir -p "$OUTDIR"
    else
        if [[ ! -d "$OUTDIR" ]]; then
            echo "The path \"$OUTDIR\" exists and is not a directory"
            exit 1
        fi
    fi
}

clean_jsonl() {
    # Sanitises string values in JSONL: collapses whitespace runs to ", " and strips quotes/edge spaces.
    jq -c 'walk(if type == "string" then gsub("[\r\n\t]+"; ", ")
           | gsub("\""; "") | sub("^ +"; "") | sub(" +$"; "")
           else . end)'
}

get_asmsumm () {
    # Downloads assembly metadata from NCBI for each taxid, writes per-taxid TSV files,
    # collects accession numbers, and merges everything into a single TSV.
    echo "Checking whether the first column name is \"accession\""
    read first_field < "$FIELDS_FILE" || true
    if [[ "${first_field}" != "accession" ]]; then
        echo "The first field in \"$FIELDS_FILE\" must be \"accession\""
        exit 1
    fi
    
    echo "Reading taxids"
    taxids=($(paste -sd " " "$TAXIDS_FILE"))
    if [[ "${#taxids[@]}" -eq 0 ]]; then
        echo "No taxids found in \"$TAXIDS_FILE\""
        exit 1
    fi

    rm -f "$OUTDIR/$ACCNS_FILE"
    
    for taxid in "${taxids[@]}"; do
        echo "Downloading assemblies summary for taxid: \"${taxid}\""
        printf -v asms_taxid_file "$OUTDIR/$ASMS_TAXID_FILE" "$taxid"
        datasets summary genome taxon "$taxid"  \
            --assembly-source GenBank           \
            --as-json-lines                     \
        |                                       \
        clean_jsonl                             \
        |                                       \
        dataformat tsv genome                   \
            --fields "$(paste -sd, "$FIELDS_FILE")" \
          > "$asms_taxid_file"
        
        echo "Extracting assembly accession numbers for taxid: \"${taxid}\""
        tail -n +2 "$asms_taxid_file" \
        |                             \
        while read -a fields; do
            echo "${fields[0]}"
        done >> "$OUTDIR/$ACCNS_FILE"
    done
    
    echo "Joining assemblies summary files"
    printf -v asms_taxid_file "$OUTDIR/$ASMS_TAXID_FILE" "${taxids[0]}"
    head -1 "$asms_taxid_file" > "$OUTDIR/$ASMS_FINAL_FILE"
    for taxid in "${taxids[@]}"; do
        printf -v asms_taxid_file "$OUTDIR/$ASMS_TAXID_FILE" "${taxid}"
        tail -n +2 "$asms_taxid_file" >> "$OUTDIR/$ASMS_FINAL_FILE"
    done
}

download_dehydrated () {
    # Downloads a dehydrated genome archive for all accessions and unzips it locally.
    echo "Downloading dehydrated genomes"
    datasets download genome accession    \
        --dehydrated                      \
        --inputfile "$OUTDIR/$ACCNS_FILE" \
        --filename  "$GEN_ARCH"
    unzip "$GEN_ARCH" -d "$GEN_DIR"
}

rehydrate_genomes () {
    # Fetches the actual sequence data for the dehydrated archive and stores it gzip-compressed.
    echo "Rehydrating genomes"
    datasets rehydrate          \
        --gzip                  \
        --directory "$GEN_DIR"
}

#-------------------------------------------------------------------------------
parse_args "$@"
create_outdir

if [[ ! -f "$OUTDIR/$ASMS_FINAL_FILE" || ! -f "$OUTDIR/$ACCNS_FILE" || "$REDO" == TRUE ]]; then
    get_asmsumm
fi

download_dehydrated
rehydrate_genomes

echo "All operations have been successfully completed"
