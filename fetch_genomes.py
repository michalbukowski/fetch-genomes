#!/usr/bin/env python3

import os, sys, argparse
import requests, urllib
import hashlib
from time import sleep
import pandas as pd


assembly_summary = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
summary_cols = 'assembly_accession taxid assembly_level asm_name ftp_path'.split()
summary_copy = 'assembly_summary_copy.tsv'
assembly_formats = {
    'fna'  : 'genomic.fna.gz',
    'gbff' : 'genomic.gbff.gz',
    'gff'  : 'genomic.gff.gz',
    'rna'  : 'rna_from_genomic.fna.gz',
    'cds'  : 'cds_from_genomic.fna.gz',
    'prot' : 'translated_cds.faa.gz'
}
assembly_levels = {
    'chr'  : 'Chromosome',
    'scff' : 'Scaffold',
    'cmpl' : 'Complete Genome',
    'ctg'  : 'Contig'
}
md5sums_fname = 'md5checksums.txt'
esearch_retmax = 100000
esearch_path = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + \
               'db=taxonomy&term=txid{taxid}[orgn]&retmode=json&retmax={retmax}&retstart={retstart}'


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--assembly-summary', type=str, default=None,
        metavar='file_path', help='Path to a custom local file in TSV format ' +
            'with assembly summary data, default: assembly summary fetched ' +
           f'from NCBI, "{assembly_summary}"')

    parser.add_argument('-c', '--summary-copy', type=str, default=summary_copy,
        metavar='file_path', help='Path where to save assembly summary for ' +
            f'chosen taxids in TSV fromat, default: {summary_copy}')

    parser.add_argument('-t', '--taxids', type=int, nargs='*', default=None,
        metavar='taxid', help='IDs of taxonomical branches to retrive genomic ' +
                              'sequences for, default: all existing(!)')

    parser.add_argument('-l', '--assembly-levels', type=str, nargs='+', default=None,
        choices=assembly_levels.keys(), metavar='level',
        help='Indicates genomes of which assembly level will be downloaded: ' +
             'chromosome (chr), scaffold (scff), complete (cmpl), ' +
             'contig (ctg), default: all levels')

    parser.add_argument('-o', '--output-dir', type=str, default='genomes',
        metavar='dir_path',
        help='Path to the directory for downloaded files, dafault: "genomes"')

    parser.add_argument('-f', '--formats', type=str, nargs='+', default=['fna'],
        choices=assembly_formats.keys(), metavar='format',
        help='Indicates files of which formats will be downloaded: ' +
             'genomic sequences in nucleotide fasta format (fna), ' +
             'genomic sequences in GenBank format (gbff), ' +
             'annotation table (gff), ' +
             'RNA sequences in nucleotide fasta format (rna), ' +
             'coding sequences (CDS) in nucleotide fasta format (cds), ' +
             'translations of CDS in protein fasta format (prot), ' +
             'default: only fna')

    parser.add_argument('-n', '--non-interactive', action='store_true',
        help='Do not ask questions and overwrite existing data ' +
             '(be absolutely sure what you do)')
    
    parser.add_argument('-s', '--summary-only', action='store_true',
        help='For given taxids or all, only download summary')

    args = parser.parse_args()

    return args


def interrogate(msg):
    ans = ''
    while ans != 'yes' and ans != 'no':
        ans = input(msg + ' (yes/no)\n')
    if ans == 'yes':
        return True
    else:
        return False


def setup_env(args):

    if args.assembly_summary is None and args.taxids is None:
        ans = True if args.non_interactive else \
              interrogate('The script will run at default values for parameters ' +
                          '--assembly-summary and --taxids, which means that ' +
                          'all genomic sequences from GenBank will be downloaded. ' +
                          'You may narrow down the number of geneomes by providing ' +
                          'either your own assembly summary file or specific taxids. ' +
                          'Do you really want to download all genomic sequences?')
        if not ans:
            return False
        else:
            print_fn('[WARNING] proceeding to download ALL genomic sequences from GenBank!')

    if os.path.exists(args.summary_copy):
        if os.path.isfile(args.summary_copy):
            ans = True if args.non_interactive else \
                  interrogate(f'Assembly summary copy "{args.summary_copy}" exists. ' +
                      'Do you want to overwrite?')
            if not ans:
                return False
        else:
            print_fn(f'[ERROR] Assembly summary copy path "{args.summary_copy}" ' +
                'points to an existing directory')
            return False

    if not args.summary_only and os.path.exists(args.output_dir):
        if os.path.isdir(args.output_dir):
            ans = True if args.non_interactive else \
                  interrogate(f'Output directory "{args.output_dir}" exists. ' +
                      'Genomes will be saved alongside existing data. ' +
                      'Do you want to continue?')
            if not ans:
                return False
        else:
            print_fn(f'[ERROR] Output directory path "{args.output_dir}" ' +
                'points to an existing file')
            return False
    os.makedirs(args.output_dir, exist_ok=True)

    return True


def fetch_taxids(taxids):

    if taxids is None:
        return [], '[INFO] No taxid provided, assembly summary will not be filtered'

    all_taxids = taxids.copy()

    for taxid in taxids:

        taxid_chunk = []

        res_chunk = [None]
        retstart = 0

        while len(res_chunk) > 0:

            res = requests.get(
                esearch_path.format(taxid=taxid, retmax=esearch_retmax, retstart=retstart),
                timeout=60
            )
            if res.status_code != 200:
                return None, f'[ERROR] A problem occured while running NCBI esearch for taxid {taxid}'

            res_json = res.json()
            if 'esearchresult' not in res_json or 'idlist' not in res_json['esearchresult']:
                return None, f'[ERROR] Unexpected result format returned by NCBI efetch for taxid {taxid}'

            res_chunk = [ int(taxid) for taxid in res_json['esearchresult']['idlist'] ]

            taxid_chunk.extend(res_chunk)
            retstart += esearch_retmax

            sleep(0.5)

        if len(taxid_chunk) == 0:
            return None, f'[ERROR] NCBI efetch returned no results for taxid {taxid}'

        all_taxids.extend(taxid_chunk)

    return all_taxids, f'[INFO] Fetched total number of {len(all_taxids)} taxids'


def fetch_summary(summary_path):
    if summary_path is None:
        try:
            res = urllib.request.urlopen(assembly_summary, timeout=60)
            summary_df = pd.read_csv(res, skiprows=1, index_col=None, sep='\t')
            summary_df.rename(columns={'# assembly_accession' : 'assembly_accession'}, inplace=True)
        except:
            return None, f'[ERROR] Assembly summary cannot be fetched from "{assembly_summary}"'

        return summary_df, '[INFO] Fetched assembly summary of {} rows and {} columns'.format(*summary_df.shape)

    else:
        try:
            summary_df = pd.read_csv(summary_path, index_col=None, sep='\t')
        except:
            return None, f'[ERROR] Assembly summary cannot be loaded from "{summary_path}"'

        if not set(summary_cols).issubset( set(summary_df.columns) ):
            return None, '[ERROR] Assembly summary does not contain required columns: ' + ', '.join(summary_cols)

        return summary_df, '[INFO] Loaded assembly summary of {} rows and {} columns'.format(*summary_df.shape)


def filter_taxids(summary_df, taxids):

    if len(taxids) == 0:
        return summary_df, f'[INFO] No taxid provided, all {summary_df.shape[0]} assemblies will be processed'

    summary_df = summary_df[ summary_df['taxid'].isin(taxids) ]

    if summary_df.shape[0] > 0:
        msg = f'[INFO] There is {summary_df.shape[0]} assemblies for the provided taxids'
    else:
        msg = '[WARNING] No assemblies in the summary for the provided taxids'

    return summary_df, msg


def filter_levels(summary_df, levels):

    if levels is None:
        return summary_df, f'[INFO] No assembly levels provided, all {summary_df.shape[0]} assemblies will be processed'
    else:
        levels = [ assembly_levels[level] for level in levels ]

    summary_df = summary_df[ summary_df['assembly_level'].isin(levels) ]

    if summary_df.shape[0] > 0:
        msg = f'[INFO] There is {summary_df.shape[0]} assemblies for the provided assembly levels'
    else:
        msg = '[WARNING] No assemblies in the summary for the provided assembly levels'

    return summary_df, msg


def save_summary(summary_df, fpath):

    try:
        summary_df.to_csv(fpath, index=False, sep='\t')
    except:
        return None, f'[ERROR] Cannot save a copy of filtered summary to "{fpath}"'

    return True, f'[INFO] Filtered summary successfully saved to "{fpath}"'


def fetch_genomes(summary_df, formats, output_dir):

    not_found = 0
    fetched   = 0
    existing  = 0

    summary_df.reset_index(drop=True, inplace=True)

    for index, (assembly_accession, ftp_path) in summary_df[
        'assembly_accession ftp_path'.split()
    ].iterrows():
    
        if ftp_path.startswith('https://'):
            ftp_path = 'ftp://' + ftp_path[8:]
        
        pos = ftp_path.rfind('/')
        asm_full_name = ftp_path[pos+1:]
        
        done = [False] * len(formats)
        for i, fmt in enumerate(formats):
            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathout = f'{output_dir}/{fnamein}'
            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    done[i] = True
        if all(done):
            existing += len(formats)
            yield f'[INFO] Skipping {assembly_accession}, already fetched'
            continue

        yield f'\n[INFO] Fetching files for assembly {assembly_accession} ' + \
              f'({index+1}/{summary_df.shape[0]})...'

        try:
            res = urllib.request.urlopen(ftp_path, timeout=60)
            lines = res.read().decode().rstrip().split('\n')
            flist = [ line.split()[-1] for line in lines ]
        except:
            yield f'[ERROR] Cannot fetch file list from "{ftp_path}"'
            yield f'[WARNING] Skipping assembly {assembly_accession}...'
            continue

        yield f'[INFO] There is {len(flist)} files at "{ftp_path}"'

        full_path = f'{ftp_path}/{md5sums_fname}'

        try:
            res = urllib.request.urlopen(full_path, timeout=60)
            md5sums = res.read().decode().rstrip().split('\n')
        except:
            yield f'[ERROR] Info on MD5 checksums cannot be fetched from "{full_path}"'
            yield f'[WARNING] Skipping assembly {assembly_accession}...'
            continue

        md5sums = [ line.split() for line in md5sums ]
        md5sums = { line[1].lstrip('./') : line[0] for line in md5sums }

        yield f'[INFO] MD5 checksums for {assembly_accession} successfully fetched'

        old_fetched = fetched

        for fmt in formats:

            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathin = f'{ftp_path}/{fnamein}'

            if not fnamein in flist:
                yield f'[ERROR] No such file for {assembly_accession} assembly: "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                not_found += 1
                continue

            fpathout = f'{output_dir}/{fnamein}'

            if not fnamein in md5sums:
                yield f'[ERROR] Cannot find MD5 checksum for {assembly_accession} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                continue

            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    yield f'[INFO] The output path "{fpathout}" exists and is a file, considered done'
                    yield f'[INFO] Skipping {assembly_accession} assembly file: "{fpathin}"'
                    existing += 1
                else:
                    yield f'[ERROR] The output path "{fpathout}" exists and is not a file'
                    yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                continue

            try:
                res = urllib.request.urlopen(fpathin, timeout=60)
                content = res.read()
            except:
                yield f'[ERROR] {assembly_accession} assembly file cannot be fetched from: "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                continue

            yield f'[INFO] {assembly_accession} assembly file "{fpathin}" successfully fetched'

            md5sum = hashlib.md5(content).hexdigest()

            if md5sum == md5sums[fnamein]:
                yield f'[INFO] Correct MD5 checksum ({md5sum}) for {assembly_accession} assembly file: "{fpathin}"'
            else:
                yield f'[ERROR] Incorrect MD5 checksum ({md5sum}) ' + \
                      f'for {assembly_accession} assembly file ({md5sums[fnamein]}): "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                continue

            tmpfpathout = f'{output_dir}/.{fnamein}'
            try:
                with open(tmpfpathout, 'wb') as f:
                    f.write(content)
            except:
                yield f'[ERROR] Cannot save to "{fpathout}" the {assembly_accession} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                continue
                
            try:
                os.rename(tmpfpathout, fpathout)
            except:
                yield f'[ERROR] Cannot save to "{fpathout}" the {assembly_accession} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {assembly_accession} assembly file: "{fpathin}"'
                os.remove(tmpfpathout)
            else:
                yield f'[INFO] {assembly_accession} assembly file "{fpathin}" successfully saved to "{fpathout}"'
                fetched += 1

    total = summary_df.shape[0] * len(formats)
    left  = total-existing-not_found-fetched
    yield f'\n[INFO] Fetched {fetched} files out of {total} inferred (not found on site: {not_found})'
    if left > 0:
        yield f'\n[WARNING] {left} files are still to be fetched'
    else:
        yield '[INFO] All files have been successfully fetched'
    yield '[INFO] Fetching genomes has been completed'

def main():
    args = parse_args()
    
    res = setup_env(args)
    if not res:
        sys.exit('[INFO] Exiting...')
    
    extra_msg = ' from NCBI (it may take a while...)' \
                if args.assembly_summary is None else ''
    print(f'[INFO] Fetching assembly summary{extra_msg}', flush=True)
    summary_df, msg = fetch_summary(args.assembly_summary)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    print('[INFO] Fetching taxids...', flush=True)
    taxids, msg = fetch_taxids(args.taxids)
    print(msg, flush=True)
    if taxids is None:
        sys.exit(1)
    
    summary_df, msg = filter_taxids(summary_df, taxids)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    summary_df, msg = filter_levels(summary_df, args.assembly_levels)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    status, msg = save_summary(summary_df, args.summary_copy)
    print(msg, flush=True)
    if status is None:
        sys.exit(1)
    elif args.summary_only:
        print('[INFO] Only assembly summary requested, exiting...', flush=True)
        sys.exit(0)
    
    for msg in fetch_genomes(summary_df, args.formats, args.output_dir):
        print(msg, flush=True)


if __name__ == '__main__':
    main()

