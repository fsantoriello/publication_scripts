#!/usr/bin/env python3

"""
Created on April 18, 2025, updated with Claude Code Opus 4.7 on April 21, 2026

@author: Francis J Santoriello

This script accepts a nucleotide fasta file and (optionally) a corresponding protein
fasta file. It then:
    1. Runs VIBRANT to identify phage sequences
    2. For --ncbi inputs: filters the companion protein fasta down to proteins whose
       parent contig survived VIBRANT
    3. For --phagedb inputs (chvd/gov2/imgvr/mgv/gpd): reformats VIBRANT's
       protein output into an NCBI-style header

This script returns:
    1. All VIBRANT output files
    2. all_phage_prots.faa (or whatever --output names it)
"""

import os
import sys
import argparse
import textwrap
import subprocess
import time
from Bio import SeqIO

### User Input ###
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
                prog='extract_vibrant_prots.py',
                formatter_class=MyFormatter,
                description=textwrap.dedent('''\
                    extract_vibrant_prots.py

                    Created on April 18, 2025

                    @author: Francis J Santoriello

                    This script accepts a nucleotide fasta file. It then:
                        1. Runs VIBRANT to identify phage sequences
                        2. Filters or reformats the associated protein fasta

                    This script returns:
                        1. All VIBRANT output files
                        2. A filtered / reformatted protein fasta
                    '''),
                usage='%(prog)s [options] <nt_fasta>\n\n '
                      'Dependencies: python 3, VIBRANT, biopython\n')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='# of parallel processes for VIBRANT.')
parser.add_argument('-p', '--path', required=True,
                    help='Path to directory containing VIBRANT_run.py')
parser.add_argument('-d', '--database', required=True,
                    help='Path to VIBRANT databases.')
parser.add_argument('--ncbi', action='store_true',
                    help='Indicates phage sequences will be compared to pre-made db '
                         'as output by phage_search_ncbi.py')
parser.add_argument('--prot_fasta',
                    help='Protein fasta file for post-processing. Required with --ncbi.')
# Per-database delimiter used by reformat_vibrant_faa to parse the contig
# accession out of VIBRANT's tab-delimited header.
DB_SPLITTERS = {'chvd': '@', 'gov2': 'length', 'imgvr': '|',
                'mgv': '_', 'gpd': '_VIRSorter'}
parser.add_argument('--phagedb', choices=sorted(DB_SPLITTERS),
                    help='Source phage database - selects the header delimiter for '
                         'the reformat path. Required unless --ncbi is set.')
parser.add_argument('-o', '--output', default='all_phage_prots.faa',
                    help='Name for output fasta file.')
parser.add_argument('--skip-vibrant', action='store_true',
                    help='Skip the VIBRANT run (useful when re-running downstream '
                         'steps against existing VIBRANTout/ output).')
parser.add_argument('nt_fasta', help='Nucleotide fasta file to search for phages.')

args = parser.parse_args()

# --- Argument validation ---
# --ncbi needs the matching protein fasta; the reformat path needs exactly one db flag.
if args.ncbi and not args.prot_fasta:
    parser.error("--prot_fasta is required when --ncbi is set")
if not args.ncbi and args.phagedb is None:
    parser.error("--phagedb is required unless --ncbi is set")


### Variables ###
seq_file = args.nt_fasta

VIBRANT_PATH = args.path
VIBRANT_DB = args.database
OUTPUT = args.output

# VIBRANT names its output dirs after the basename of the input (no extension),
# so strip both the directory and the .fna suffix to build the expected output path.
seqs = os.path.splitext(os.path.basename(seq_file))[0]
VIBRANT_outpath = f"VIBRANTout/VIBRANT_{seqs}/VIBRANT_phages_{seqs}/{seqs}"


### Modules ###
def run(cmd):
    """Run a shell command and raise CalledProcessError on non-zero exit."""
    subprocess.run(cmd, shell=True, check=True)


def extract_corr_prots(nuc_fasta, prot_fasta, output_path):
    """
    Filter prot_fasta down to records whose parent nt accession is present in
    nuc_fasta. The nt accession is parsed from the FASTA description, expected
    to contain `[<accession> | ...]` as written by phage_search_ncbi.py.
    Args:
       nuc_fasta:   Nucleotide fasta (VIBRANT's phages_combined.fna).
       prot_fasta:  Companion protein fasta produced by phage_search_ncbi.py.
       output_path: Destination fasta file for the filtered proteins.
    """
    # Collect surviving contig ids from the nt fasta.
    with open(nuc_fasta, "r") as seqfile:
        nt_accs = {record.id for record in SeqIO.parse(seqfile, 'fasta')}

    # Walk the protein fasta and keep any record whose bracketed accession is in nt_accs.
    with open(prot_fasta, 'r') as prot_db, open(output_path, 'w') as outfile:
        for record in SeqIO.parse(prot_db, 'fasta'):
            acc_start = record.description.find('[')
            acc_end = record.description.find(' |')
            if acc_start == -1 or acc_end == -1:
                continue  # header doesn't match the expected "[acc | ...]" layout
            acc_id = record.description[acc_start + 1:acc_end]

            if acc_id in nt_accs:
                outfile.write(f">{record.description}\n{record.seq}\n")


def reformat_vibrant_faa(vib_prot_fasta, splitter, output_path):
    """
    Rewrite VIBRANT protein FASTA headers into an NCBI-style format. The incoming
    description is expected to be tab-delimited; column 0 carries the contig id
    (with a per-db delimiter we split on) and columns 3-4 carry the protein id
    and product annotation.
    Args:
       vib_prot_fasta: VIBRANT phages_combined.faa.
       splitter:       Per-database delimiter used to recover the contig accession.
       output_path:    Destination fasta file for reformatted records.
    """
    with open(vib_prot_fasta, "r") as protfile, open(output_path, 'w') as outfile:
        for record in SeqIO.parse(protfile, 'fasta'):
            annos = record.description.split('\t')
            if len(annos) < 5:
                # Non-prodigal / unexpected header layout - skip rather than crash.
                print(f"Skipping {record.id}: expected >=5 tab fields, got {len(annos)}",
                      file=sys.stderr)
                continue

            nt_acc = annos[0].split(splitter)[0]
            gene_num = annos[0].split('_')[-1]
            prot_id = annos[3]
            new_entry_id = f"{gene_num}_{prot_id}_{nt_acc}"

            new_entry_product = annos[4]

            # Strip the trailing _N gene index to recover the organism/contig label.
            org_index = annos[0].rfind('_')
            org = annos[0][:org_index] if org_index != -1 else annos[0]
            new_entry_desc = f"[{org}]"

            header = f">{new_entry_id} {new_entry_product} {new_entry_desc}"
            outfile.write(f"{header}\n{record.seq}\n")


### Main ###

if __name__ == "__main__":
    if not args.skip_vibrant:
        print("Running VIBRANT...")
        start_time = time.time()
        run(f"python {os.path.join(VIBRANT_PATH, 'VIBRANT_run.py')} "
            f"-i {seq_file} -folder VIBRANTout "
            f"-t {args.threads} -no_plot -d {VIBRANT_DB}")
        elapsed_time = time.time() - start_time
        print(f"VIBRANT complete. Elapsed time: {elapsed_time:.4f} seconds")
    else:
        print("Skipping VIBRANT run (--skip-vibrant).")

    if args.ncbi:
        # Filter the companion protein fasta to proteins whose contig survived VIBRANT.
        phage_seqs = f"{VIBRANT_outpath}.phages_combined.fna"
        print("Filtering proteins...")
        extract_corr_prots(phage_seqs, args.prot_fasta, OUTPUT)
    else:
        # Reformat VIBRANT's own protein fasta into NCBI-style headers.
        phage_prots = f"{VIBRANT_outpath}.phages_combined.faa"
        splitter = DB_SPLITTERS[args.phagedb]
        print("Reformatting VIBRANT protein fasta...")
        reformat_vibrant_faa(phage_prots, splitter, OUTPUT)

    print("Process completed.")
