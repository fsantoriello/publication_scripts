#!/usr/bin/env python3

"""
Created on May 2, 2025, , updated with Claude Code Opus 4.7 on April 21, 2026

@author: Francis J Santoriello

This script accepts a nucleotide fasta file and a corresponding protein fasta file.
It then:
    1. Generates a subset of the protein fasta containing only the sequences
       corresponding to sequences in the nt fasta

This script returns:
    1. filtered_prots.faa
"""

import argparse
import textwrap
from Bio import SeqIO

### User Input ###
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
                prog='fetch_final_prots.py',
                formatter_class=MyFormatter,
                description=textwrap.dedent('''\
                    fetch_final_prots.py

                    Created on May 2, 2025

                    @author: Francis J Santoriello

                    This script accepts a nucleotide fasta file and a corresponding
                    protein fasta file. It then:
                    1. Generates a subset of the protein fasta containing only the
                       sequences corresponding to sequences in the nt fasta

                    This script returns:
                    1. filtered_prots.faa
                    '''),
                usage='%(prog)s [options] <nt_fasta> <prot_fasta>\n\n '
                      'Dependencies: python 3, biopython\n')
parser.add_argument('--ncbi', action='store_true',
                    help='Indicates ncbi-formatted headers (bracketed accession '
                         'followed by " |"). Otherwise closing "]" is used.')
parser.add_argument('-o', '--output', default='filtered_prots.faa',
                    help='Name for output fasta file.')
parser.add_argument('nt_fasta', help='Nucleotide fasta file to filter against.')
parser.add_argument('prot_fasta', help='Protein fasta file for filtering.')

args = parser.parse_args()


### Variables ###
seq_file = args.nt_fasta
prot_file = args.prot_fasta
OUTPUT = args.output


### Modules ###
def parent_accession(description, ncbi):
    """
    Extract the parent nucleotide accession from a protein FASTA description.
    For NCBI-style headers the accession sits between '[' and ' |'; for the
    in-house/reformatted headers it sits between '[' and ']'. Returns None if
    the expected delimiters are missing.
    """
    start = description.find('[')
    end = description.find(' |') if ncbi else description.find(']')
    if start == -1 or end == -1 or end <= start:
        return None
    return description[start + 1:end]


def extract_prots(nuc_fasta, prot_fasta, output_path, ncbi):
    """
    Filter prot_fasta down to records whose parent nucleotide accession (parsed
    from the FASTA description) is present in nuc_fasta.
    Args:
       nuc_fasta:   A fasta file of nucleotide sequences.
       prot_fasta:  A fasta file of corresponding protein sequences.
       output_path: Destination fasta file for the filtered proteins.
       ncbi:        Use NCBI-style bracket parsing ('[acc | ...]') when True.
    """
    # Collect surviving contig ids from the nt fasta as a set for O(1) lookups.
    with open(nuc_fasta, "r") as seqfile:
        nuc_ids = {record.id for record in SeqIO.parse(seqfile, 'fasta')}

    # Walk the protein fasta and keep any record whose parent accession is in nuc_ids.
    with open(prot_fasta, 'r') as prot_db, open(output_path, 'w') as outfile:
        for record in SeqIO.parse(prot_db, 'fasta'):
            acc = parent_accession(record.description, ncbi)
            if acc is None or acc not in nuc_ids:
                continue
            outfile.write(f">{record.description}\n{record.seq}\n")


### Main ###

if __name__ == "__main__":
    print("Filtering sequences...")
    extract_prots(seq_file, prot_file, OUTPUT, args.ncbi)
    print("Process completed.")
