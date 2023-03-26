#!/usr/bin/env python3
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# USAGE, example (for more see ./fetch_genomes.py -h):
# Download genomes for taxid 1279 (Staphylococcus) and 1350 (Enterococcus) and
# all subtaxa to a default directory (genomes) in the current location:
# ./fetch_genomes.py -t 1279 1350
# Resume previous downolading based on saved filtered assembly summary:
# ./fetch_genomes.py -a assembly_summary_copy.tsv
# Retrive filtered assembly summary only:
# ./fetch_genomes.py -t 1279 1350 -s
# You may find desirable taxids here: https://www.ncbi.nlm.nih.gov/taxonomy

import os, sys, argparse
import requests, urllib
import hashlib
from time import sleep
import pandas as pd

# A path to a TSV file on NCBI server that contains info on genomic assemblies.
assembly_summary = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt'
# Columns that must be present in the summary file.
summary_cols = 'assembly_accession taxid assembly_level asm_name ftp_path'.split()
# The default path where to save a filtered summary copy.
summary_copy = 'assembly_summary_copy.tsv'
# Data that can be obtained from NCBI GenBank for a given genomic assembly,
# see parse_args function for more information.
assembly_formats = {
    'fna'  : 'genomic.fna.gz',
    'gbff' : 'genomic.gbff.gz',
    'gff'  : 'genomic.gff.gz',
    'rna'  : 'rna_from_genomic.fna.gz',
    'cds'  : 'cds_from_genomic.fna.gz',
    'prot' : 'translated_cds.faa.gz'
}
# Possible genomic assembly levels to choose from, see parse_args function
# for more information.
assembly_levels = {
    'chr'  : 'Chromosome',
    'scff' : 'Scaffold',
    'cmpl' : 'Complete Genome',
    'ctg'  : 'Contig'
}
# The default directory where to save downloaded genomes.
gen_dir = 'genomes'
# The filename with MD5 checksums in remote assembly directories.
md5sums_fname = 'md5checksums.txt'
# The maximum number of records to fetch from NCBI Taxonomy.
esearch_retmax = 100000
# A request template to NCBI Taxonomy database.
esearch_path = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + \
               'db=taxonomy&term=txid{taxid}[orgn]&retmode=json&'            + \
               'retmax={retmax}&retstart={retstart}'

def parse_args():
    '''Parse arguments:
       -a, ‑‑assembly‑summary -- a path to a custom local file in TSV format that
           contains information on assemblies that are to be downloaded,
           default: assembly summary will be fetched from NCBI GenBank FTP site
       -c, ‑‑summary-copy -- a path to a TSV file where to save the filtered
           assembly summary for chosen taxids in TSV format,
           default: _assembly_summary_copy.tsv_
       -t, ‑‑taxids -- space-separated IDs of taxa to retrive genomic sequences for,
           default: all existing(!)
       -l, ‑‑assembly-levels -- space-separated assembly levels that will be taken
           into consideration: chromosome (chr), scaffold (scff), complete (cmpl),
           contig (ctg), default: all levels
       -o, ‑‑output-dir -- a path to the directory for downloaded genomes,
           dafault: _genomes_
       -f, ‑‑formats -- formats of data to be downloaded: genomic sequences in
           nucleotide fasta format (fna), genomic sequences in GenBank format (gbff),
           annotation table (gff), RNA sequences in nucleotide fasta format (rna),
           coding sequences (CDS) in nucleotide fasta format (cds), translations
           of CDS in protein fasta format (prot), default: fna
       -n, ‑‑non-interactive -- do not ask questions and overwrite existing data
           (be absolutely sure what you do)
       -s, ‑‑summary-only -- for given taxids or all, only download assembly summary
    '''
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '-a', '--assembly-summary', type=str, default=None, metavar='file_path',
        help='A path to a custom local file in TSV format that contains' +
             ' information on assemblies that are to be downloaded, default: ' +
            f'assembly summary will be fetched from NCBI, "{assembly_summary}"'
    )
    parser.add_argument(
        '-c', '--summary-copy', type=str, default=summary_copy, metavar='file_path',
        help='A path to a TSV file where to save the filtered assembly ' +
            f'summary for chosen taxids in TSV format, default: {summary_copy}'
    )
    parser.add_argument(
        '-t', '--taxids', type=int, nargs='*', default=None, metavar='taxid',
        help='Space-separated IDs of taxa to retrive genomic sequences for, ' +
             'default: all existing(!)'
    )
    parser.add_argument(
        '-l', '--assembly-levels', type=str, nargs='+', default=None,
        choices=assembly_levels.keys(), metavar='level',
        help='Space-separated assembly levels that will be taken into' +
            'consideration: chromosome (chr), scaffold (scff), ' +
            'complete (cmpl), contig (ctg), default: all levels'
    )
    parser.add_argument(
        '-o', '--output-dir', type=str, default=gen_dir, metavar='dir_path',
        help='Path to the directory for downloaded genomes, dafault: "genomes"'
    )
    parser.add_argument(
        '-f', '--formats', type=str, nargs='+', default=['fna'],
        choices=assembly_formats.keys(), metavar='format',
        help='Formats of data to be downloaded: ' +
             'genomic sequences in nucleotide fasta format (fna), ' +
             'genomic sequences in GenBank format (gbff), ' +
             'annotation table (gff), ' +
             'RNA sequences in nucleotide fasta format (rna), ' +
             'coding sequences (CDS) in nucleotide fasta format (cds), ' +
             'translations of CDS in protein fasta format (prot), ' +
             'default: only fna'
    )
    parser.add_argument(
        '-n', '--non-interactive', action='store_true',
        help='Do not ask questions and overwrite existing data ' +
             '(be absolutely sure what you do)'
    )
    parser.add_argument(
        '-s', '--summary-only', action='store_true',
        help='For given taxids or all, only download assembly summary'
    )
    
    args = parser.parse_args()
    return args

def interrogate(msg):
    '''In interactive mode, shows a yes/no question (msg) and retrive an answer,
       not used when non-interatvie mode (-n, --non-interactive) is on. Arguments:
       msg -- a message/question to be communicated
    '''
    ans = ''
    while ans != 'yes' and ans != 'no':
        ans = input(msg + ' (yes/no)\n')
    if ans == 'yes':
        return True
    else:
        return False

def setup_env(args):
    '''Shows information on planned actions and ask for confirmation, unless
       non-interatvie mode (-n, --non-interactive) is on. Arguments:
       args -- the object returned by parse_args function containg values
               of command line arguments
    '''
    # No taxid or pre-filtered assembly summary is provided, ask whether to
    # download all genomes from NCBI GenBank.
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
            print('[WARNING] proceeding to download ALL genomic sequences from GenBank!')
    
    # Check wheteher indicated assembly summary path exists, if not or it is
    # not a file, show an error message and return False.
    if args.assembly_summary is not None:
        if not os.path.exists(args.assembly_summary):
            print(f'[ERROR] Assembly summary path "{args.assembly_summary}" ' +
                'does not exist')
            return False
        elif not os.path.isfile(args.assembly_summary):
            print(f'[ERROR] Assembly summary path "{args.assembly_summary}" ' +
                'points to an existing directory')
            return False
    
    # Confirm wheteher to overvire an existing assembly summary or show
    # an error message and return False if the path point to a directory.
    if os.path.exists(args.summary_copy):
        if os.path.isfile(args.summary_copy):
            ans = True if args.non_interactive else \
                  interrogate(f'Assembly summary copy "{args.summary_copy}" exists. ' +
                      'Do you want to overwrite?')
            if not ans:
                return False
        else:
            print(f'[ERROR] Assembly summary copy path "{args.summary_copy}" ' +
                'points to an existing directory')
            return False
    
    # Confirm whether to download genomes to an existing directory, show
    # an error message and return False if the path points to a file.
    if not args.summary_only and os.path.exists(args.output_dir):
        if os.path.isdir(args.output_dir):
            ans = True if args.non_interactive else \
                  interrogate(f'Output directory "{args.output_dir}" exists. ' +
                      'Genomes will be saved alongside existing data. ' +
                      'Do you want to continue?')
            if not ans:
                return False
        else:
            print(f'[ERROR] Output directory path "{args.output_dir}" ' +
                'points to an existing file')
            return False
    
    # Create the output diretory for download genomes, if does not exist, unless
    # only summary is to be fetched.
    if not args.summary_only:
        os.makedirs(args.output_dir, exist_ok=True)
    return True


def fetch_taxids(taxids):
    '''Fetches IDs for all subtaxa for provided taxids from NCBI Taxonomy in chunks
       of esearch_retmax (hardcoded for maximum size of 100000) by sending
       HTTPS requests with GET method and obtaining responses in JSON format,
       sleeps 0.5 sec after each request to avoid being blocked by NCBI server.
       Arguments:
       taxids -- a list of taxids of interest, may be empty if all taxids ought
                 to be taken into account
    '''
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
    '''Fetches assembly summary form NCBI GenBank FTP server if a path to an
       exisiting one is not given, returns None and an error message when
       any problem arises. Arguments:
       summary_path -- a path to a custom local file in TSV format that
                       contains information on assemblies that are to be downloaded,
                       if None, the summary will be fetch from NCBI server
    '''
    if summary_path is None:
        try:
            res = urllib.request.urlopen(assembly_summary, timeout=60)
            summary_df = pd.read_csv(res, skiprows=1, index_col=None, sep='\t')
            summary_df.rename(columns={'# assembly_accession' : 'assembly_accession'}, inplace=True)
        except:
            return None, f'[ERROR] Assembly summary cannot be fetched from "{assembly_summary}"'
        else:
            return summary_df, '[INFO] Fetched assembly summary of {} rows and {} columns'.format(*summary_df.shape)
    else:
        try:
            summary_df = pd.read_csv(summary_path, index_col=None, sep='\t')
        except:
            return None, f'[ERROR] Assembly summary cannot be loaded from "{summary_path}"'
            
        if not set(summary_cols).issubset( set(summary_df.columns) ):
            return None, '[ERROR] Assembly summary does not contain required columns: ' + ', '.join(summary_cols)
        else:
            return summary_df, '[INFO] Loaded assembly summary of {} rows and {} columns'.format(*summary_df.shape)

def filter_taxids(summary_df, taxids):
    '''Having a complete list of requensted taxids and taxids for all subtaxa,
       selects desirable rows from the assembly summary DataFrame. Arguments:
       summary_df -- a Pandas DataFrame with the assembly summary
       taxids     -- a list of taxids of interest
    '''
    if len(taxids) == 0:
        return summary_df, f'[INFO] No taxid provided, all {summary_df.shape[0]} assemblies will be processed'
        
    summary_df = summary_df[ summary_df['taxid'].isin(taxids) ]
    
    if summary_df.shape[0] > 0:
        msg = f'[INFO] There is {summary_df.shape[0]} assemblies for the provided taxids'
    else:
        msg = '[WARNING] No assemblies in the summary for the provided taxids'
        
    return summary_df, msg


def filter_levels(summary_df, levels):
    '''Having assembly levels indicated, selects desirable rows from the assembly
       summary DataFrame. Arguments:
       summary_df -- a Pandas DataFrame with the assembly summary
       levels     -- a list of assembly levels of interest
    '''
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
    '''Saves the filtered assembly summary. Arguments:
       summary_df -- a Pandas DataFrame with the final assembly summary
       fpath      -- a path where the summary ought to be saved as a TSV file
    '''
    try:
        summary_df.to_csv(fpath, index=False, sep='\t')
    except:
        return None, f'[ERROR] Cannot save a copy of filtered summary to "{fpath}"'
    else:
        return True, f'[INFO] Filtered summary successfully saved to "{fpath}"'

def fetch_genomes(summary_df, formats, output_dir):
    '''Having filtered assembly summary, fetches finally selected genomes by sending
       FTP requests via urllib. Arguments:
       summary_df -- a Pandas DataFrame with the final assembly summary
       formats    -- formats of data to be retrieved
       output_dir -- a directory for the data to be saved to
    '''
    not_found = 0
    fetched   = 0
    existing  = 0
    
    # Iterate over assembly accession numbers and corresponding FTP paths
    # in the assembly summary DataFrame, if an error appears, show a message
    # and proceed to next genome.
    summary_df.reset_index(drop=True, inplace=True)
    for index, (asm_acc, ftp_path) in summary_df[
        'assembly_accession ftp_path'.split()
    ].iterrows():
        if ftp_path.startswith('https://'):
            ftp_path = 'ftp://' + ftp_path[8:]
        pos = ftp_path.rfind('/')
        asm_full_name = ftp_path[pos+1:]
        
        # Check whether all requested files are already downloaded, if so,
        # yield a proper message and continue to next iteration/genome.
        done = [False] * len(formats)
        for i, fmt in enumerate(formats):
            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathout = f'{output_dir}/{asm_acc}_{suffix}'
            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    done[i] = True
        if all(done):
            existing += len(formats)
            yield f'[INFO] All files requested for {asm_acc} exist and are files, considered done'
            yield f'[INFO] Skipping {asm_acc}, already fetched'
            continue
        yield f'\n[INFO] Fetching files for assembly {asm_acc} ' + \
            f'({index+1}/{summary_df.shape[0]})...'
        
        # Fetch file list form the genome directory, if unsuccessful, yield
        # a proper message and continue to next iteration/genome.
        try:
            res = urllib.request.urlopen(ftp_path, timeout=60)
            lines = res.read().decode().rstrip().split('\n')
            flist = [ line.split()[-1] for line in lines ]
        except KeyboardInterrupt as e:
            raise e
        except:
            yield f'[ERROR] Cannot fetch file list from "{ftp_path}"'
            yield f'[WARNING] Skipping assembly {asm_acc}...'
            continue
        yield f'[INFO] There is {len(flist)} files at "{ftp_path}"'
        
        # Fetch the file with MD5 sums for genome files, if unsuccessful, yield
        # a proper message and continue to next iteration/genome.
        full_path = f'{ftp_path}/{md5sums_fname}'
        try:
            res = urllib.request.urlopen(full_path, timeout=60)
            md5sums = res.read().decode().rstrip().split('\n')
        except KeyboardInterrupt as e:
            raise e
        except:
            yield f'[ERROR] Info on MD5 checksums cannot be fetched from "{full_path}"'
            yield f'[WARNING] Skipping assembly {asm_acc}...'
            continue
        md5sums = [ line.split() for line in md5sums ]
        md5sums = { line[1].lstrip('./') : line[0] for line in md5sums }
        yield f'[INFO] MD5 checksums for {asm_acc} successfully fetched'
        
        # Iterate over per-genome requested files (data formats) and
        # fetch those files.
        old_fetched = fetched
        for fmt in formats:
            suffix  = assembly_formats[fmt]
            fnamein = f'{asm_full_name}_{suffix}'
            fpathin = f'{ftp_path}/{fnamein}'
            
            # Look up wheter request file exists on the server, if not, yield
            # a proper message and continue to next iteration/genome.
            if not fnamein in flist:
                yield f'[ERROR] No such file for {asm_acc} assembly: "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                not_found += 1
                continue
            
            # Continue to next iteration/genome if requested file already exists
            # in the ouput directory, yeild a proper message if the local path
            # points to a directory.
            fpathout = f'{output_dir}/{asm_acc}_{suffix}'
            if os.path.exists(fpathout):
                if os.path.isfile(fpathout):
                    yield f'[INFO] The output path "{fpathout}" exists and is a file, considered done'
                    yield f'[INFO] Skipping {asm_acc} assembly file: "{fpathin}", already fetched'
                    existing += 1
                else:
                    yield f'[ERROR] The output path "{fpathout}" exists and is not a file'
                    yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                continue
            
            # Check if there is a MD5 sum for the file to be downloaded,
            # if not, yeild an error message and continue to next iteration.
            if not fnamein in md5sums:
                yield f'[ERROR] Cannot find MD5 checksum for {asm_acc} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                continue
            
            # Fetch the content of the file to be downloaded, yield a proper
            # message if it goes wrong and continue to next iteration/genome.
            try:
                res = urllib.request.urlopen(fpathin, timeout=60)
                content = res.read()
            except:
                yield f'[ERROR] {asm_acc} assembly file cannot be fetched from: "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                continue
            yield f'[INFO] {asm_acc} assembly file "{fpathin}" successfully fetched'
            
            # Generate an MD5 sum for the downloaded content and compare to the one
            # retrieved from the server, if these do not agree, yield a proper
            # message and continue to next iteration/genome.
            md5sum = hashlib.md5(content).hexdigest()
            if md5sum == md5sums[fnamein]:
                yield f'[INFO] Correct MD5 checksum ({md5sum}) for {asm_acc} assembly file: "{fpathin}"'
            else:
                yield f'[ERROR] Incorrect MD5 checksum ({md5sum}) ' + \
                      f'for {asm_acc} assembly file ({md5sums[fnamein]}): "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                continue
            
            # Try to save content to a temporary file, if unsuccessful, yield
            # a proper message and continue to next iteration/genome.
            tmpfpathout = f'{output_dir}/.{asm_acc}_{suffix}'
            try:
                with open(tmpfpathout, 'wb') as f:
                    f.write(content)
            except:
                yield f'[ERROR] Cannot save to "{fpathout}" the {asm_acc} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                continue
            
            # Try to rename the temporary file to give it the final name,
            # if unsuccessful, yield a proper message.
            try:
                os.rename(tmpfpathout, fpathout)
            except:
                yield f'[ERROR] Cannot save to "{fpathout}" the {asm_acc} assembly file: "{fpathin}"'
                yield f'[WARNING] Skipping {asm_acc} assembly file: "{fpathin}"'
                os.remove(tmpfpathout)
            else:
                yield f'[INFO] {asm_acc} assembly file "{fpathin}" successfully saved to "{fpathout}"'
                fetched += 1
    
    # Show the summary, especially how many files already existed or were
    # successfully fetch as well as fetching of how many failed and require
    # a rerun.
    total = summary_df.shape[0] * len(formats)
    left  = total-existing-not_found-fetched
    yield f'\n[INFO] Fetched {fetched} files out of {total} inferred ' + \
          f'(already existing: {existing}, not found on site: {not_found})'
    if left > 0:
        yield f'\n[WARNING] {left} files are still to be fetched'
    else:
        yield '[INFO] All files have been successfully fetched'
    yield '[INFO] Fetching genomes has been completed'

def main():
    '''The entry poin function that executes all stages one by one, receive messages
       from stage-executing functions and prints them, if a critical error araises,
       shows a message and exits.
    '''
    # Parse command line arguments.
    args = parse_args()
    
    # Setup all variables, check up the environment, a function that
    # interactacts with a user unless non-interatvie mode is on
    # (-n, --non-interactive),
    # the only function that pronts messages on its own.
    res = setup_env(args)
    if not res:
        sys.exit('[INFO] Exiting...')
    
    # Fetch assembly summary from NCBI GenBank FTP server, if a path to
    # an existing one is not provided.
    extra_msg = ' from NCBI (it may take a while...)' \
                if args.assembly_summary is None else ''
    print(f'[INFO] Fetching assembly summary{extra_msg}', flush=True)
    summary_df, msg = fetch_summary(args.assembly_summary)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    # Fetch taxids for subtaxa of provided taxids, if the list is empty,
    # fetch_taxids function handles it properly.
    print('[INFO] Fetching taxids...', flush=True)
    taxids, msg = fetch_taxids(args.taxids)
    print(msg, flush=True)
    if taxids is None:
        sys.exit(1)
    
    # Filter assembly summary DataFrame to keep genomes belonging to taxa
    # contained in the fetch taxid list, if the list is empty
    # filter_taxids function handles it properly.
    summary_df, msg = filter_taxids(summary_df, taxids)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    # Filter assembly summary DataFrame to keep genomes of chosen assembly levels.
    summary_df, msg = filter_levels(summary_df, args.assembly_levels)
    print(msg, flush=True)
    if summary_df is None:
        sys.exit(1)
    
    # Save filtered summary, i.e. summary that contain only those genomes that
    # are to be fetched, exit if only summary has been requested, wihtout
    # fetching any data.
    status, msg = save_summary(summary_df, args.summary_copy)
    print(msg, flush=True)
    if status is None:
        sys.exit(1)
    elif args.summary_only:
        print('[INFO] Only assembly summary requested, exiting...', flush=True)
        sys.exit(0)
    
    # Iteratively fetch selected genomes using fetch_genomes generator that
    # yields messages on the progress.
    for msg in fetch_genomes(summary_df, args.formats, args.output_dir):
        print(msg, flush=True)

# Entry point.
if __name__ == '__main__':
    main()

