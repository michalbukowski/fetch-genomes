## fetch-genomes

### Requirements
The script should work in all Python 3 environents with Pandas library installed. Below I provide versions for which I tested the script:
- Python (3.11)
- Pandas (1.5.3)

### Short description

The script allows for downloading genomes from [NCBI GenBank FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/genbank), based on the content of [`assembly_summary_genbank.txt`](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt) or any TSV file that provides desirable data on genome assemblies in the following colums: `assembly_accession`, `taxid`, `assembly_level`, asm_name and `ftp_path`. For more information see [`README.txt`](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/README.txt) at the NCBI GenBank FTP site.

You can run it without any command line options or provide one or more `taxid` values to narrow down the number of genomes you want to retrieve. Taxa IDs of your interest you may find in [NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy). Genomes of requested taxa and all subordinate subtaxa will be retrieved.

### Example usage
- Download genomes for taxid 1279 (genus _Staphylococcus_) and 1350 (genus _Enterococcus_) and all subtaxa, i.e. all genomes assigned to the genera as well as all subordinate species, subspecies etc. Save the genomes to a default directory (`genomes`) in the current location:
```bash
./fetch_genomes.py -t 1279 1350
```
- If you have problems with network connection, you may rerun the script until all genomes are successfully retrived, i.e. when you see in the end a message saying: `[INFO] All files have been successfully fetched`. Simply resume previous downolading or retry to download skipped genomes based on saved filtered assembly summary from a previous search (existing files will _not_ be redownloaded):
```bash
./fetch_genomes.py -a assembly_summary_copy.tsv
```
- Retrive filtered assembly summary only, i.e. assembly summary on genomes belonging to requested taxa, without downloading anything else in order to examine it or modify before use:
```bash
./fetch_genomes.py -t 1279 1350 -s
```

### Command line options
| Option | Use |
|---|---|
| `‑a`,&nbsp;`‑‑assembly‑summary` | A path to a custom local file in TSV format that contains information on assemblies that are to be downloaded, default: assembly summary will be fetched from NCBI GenBank FTP site|
| `‑c`,&nbsp;`‑‑summary-copy` | A path to a TSV file where to save the filtered assembly summary for chosen taxids in TSV format, default: _assembly_summary_copy.tsv_|
| `‑t`,&nbsp;`‑‑taxids` | Space-separated IDs of taxa to retrive genomic sequences for, default: all existing(!) |
| `‑l`,&nbsp;`‑‑assembly-levels` | Space-separated assembly levels that will be taken into consideration: chromosome (`chr`), scaffold (`scff`), complete (`cmpl`), contig (`ctg`), default: all levels |
| `‑o`,&nbsp;`‑‑output-dir` | A path to the directory for downloaded genomes, dafault: _genomes_ |
| `‑f`,&nbsp;`‑‑formats` | Formats of data to be downloaded: genomic sequences in nucleotide fasta format (`fna`), genomic sequences in GenBank format (`gbff`), annotation table (`gff`), RNA sequences in nucleotide fasta format (`rna`), coding sequences (`CDS`) in nucleotide fasta format (`cds`), translations of CDS in protein fasta format (`prot`), default: fna |
| `‑n`,&nbsp;`‑‑non-interactive` | Do not ask questions and overwrite existing data (be absolutely sure what you do) |
| `‑s`,&nbsp;`‑‑summary-only` | For given taxids or all, only download assembly summary |

