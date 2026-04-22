#!/usr/bin/env python3

'''
Created on April 21, 2025, updated with Claude Code Opus 4.7 on April 21, 2026

@author: Francis J Santoriello

This script accepts a protein fasta file, a corresponding nucleotide fasta, and HMM
profiles. It then:
    1. Runs HMM searches for all HMM profiles in the working directory against the custom
       protein database
    2. Extracts all other proteins from the same contig as the hit into a new protein
       fasta (all proteins from all contigs with hits)
    3. Extracts the nucleotide sequence corresponding to each hit into a new nucleotide
       fasta file.

This script returns:
    1. {inputdb}_{hmm_profile}.faa
    2. {inputdb}_{hmm_profile}.fna
'''

import os
import sys
import argparse
import subprocess
import textwrap
from Bio import SeqIO
from Bio import SearchIO

### User Input ###
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
                prog='HMM_search_phage_dbs.py',
                formatter_class=MyFormatter,
                description=textwrap.dedent('''\
                    HMM_search_phage_dbs.py

                    Created on April 21, 2025

                    @author: Francis J Santoriello

                    This script accepts a protein fasta file, a corresponding nucleotide fasta,
                    and HMM profiles. It then:
                    1. Runs HMM searches for all HMM profiles in the working directory against
                       the custom protein database
                    2. Extracts all other proteins from the same contig as the hit into a new
                       protein fasta (all proteins from all contigs with hits)
                    3. Extracts the nucleotide sequence corresponding to each hit into a new
                       nucleotide fasta file.

                    This script returns:
                    1. {inputdb}_{hmm_profile}.faa
                    2. {inputdb}_{hmm_profile}.fna
                    '''),
                usage='%(prog)s [options] <prot_fasta> <nt_fasta> <HMM_profiles>\n\n '
                      'Dependencies: python 3, HMMER, biopython\n')
parser.add_argument('-l', '--len', type=int, default=15000,
                    help='Minimum nucleotide sequence length for filtering')
parser.add_argument('-e', '--evalue', type=float, default=1e-8,
                    help='hmmsearch E-value threshold')
parser.add_argument('prot_fasta', help='Protein fasta file, database for hmmsearch')
parser.add_argument('nt_fasta',
                    help='Nucleotide fasta file, database for filtering based on hmmsearch')
parser.add_argument('hmm', nargs='+', help='HMM profiles for hmmsearch')

args = parser.parse_args()


### Variables ###
protdb = args.prot_fasta
nucdb = args.nt_fasta


### Modules ###
def run(cmd):
    """Run a shell command and raise CalledProcessError on non-zero exit."""
    subprocess.run(cmd, shell=True, check=True)


def parent_contig(seq_id):
    """
    Strip the trailing `_N` CDS index (prodigal-style IDs) to recover the parent
    contig accession. Returns the original id unchanged if no underscore is found.
    """
    idx = seq_id.rfind('_')
    return seq_id[:idx] if idx != -1 else seq_id


def read_hit_contigs(hmmsearch_table):
    """
    Parse an hmmsearch --tblout file and return the set of parent contig ids
    corresponding to every hit target.
    """
    hit_contigs = set()
    for record in SearchIO.parse(hmmsearch_table, 'hmmer3-tab'):
        for hit in record.hits:
            hit_contigs.add(parent_contig(hit.id))
    return hit_contigs


def extract_hmm_hits_to_faa(hmmsearch_table, prot_fasta, output):
    '''
    Takes an hmmsearch hit table and extracts all proteins on the hit contig from the
    queried protein database into a new fasta file.
    '''
    hit_contigs = read_hit_contigs(hmmsearch_table)

    with open(prot_fasta, 'r') as prot_handle, open(output, 'w') as outfile:
        for record in SeqIO.parse(prot_handle, 'fasta'):
            rec_contig = parent_contig(record.id)
            if rec_contig not in hit_contigs:
                continue

            # Prodigal-style FASTA headers put the product annotation in the 5th
            # '#'-delimited field of the description. Fall back gracefully for
            # other sources.
            annos = record.description.split(' # ')
            product = annos[4] if len(annos) >= 5 else ''
            header = f">{record.id} {product} [{rec_contig}]"
            outfile.write(f"{header}\n{record.seq}\n")


def extract_hmm_hits_to_fna(hmmsearch_table, nt_fasta, output, min_len):
    '''
    Takes an hmmsearch hit table and extracts the nucleotide sequence of each hit contig
    (above min_len) into a new fasta file.
    '''
    hit_contigs = read_hit_contigs(hmmsearch_table)

    with open(nt_fasta, 'r') as nuc_handle, open(output, 'w') as outfile:
        for record in SeqIO.parse(nuc_handle, 'fasta'):
            if len(record.seq) < min_len:
                continue
            if record.id not in hit_contigs:
                continue
            outfile.write(f">{record.id}\n{record.seq}\n")


### Main ###

if __name__ == '__main__':
    # Drop the .faa suffix from the protein db path to use as an output prefix.
    protdb_short = os.path.splitext(protdb)[0]

    for hmm in args.hmm:
        # Name outputs after the profile basename rather than the full path.
        hmm_reform = os.path.basename(hmm)
        hmm_short = os.path.splitext(hmm_reform)[0]
        hmm_table = f"{protdb_short}_{hmm_short}_hits_table"

        # Run hmmsearch against the protein database with the user-specified E-value.
        print(f"Running hmmsearch for {hmm_reform}...")
        try:
            run(f"hmmsearch --tblout {hmm_table} -E {args.evalue} {hmm} {protdb}")
        except subprocess.CalledProcessError as e:
            print(f"hmmsearch failed for {hmm_reform}: {e}", file=sys.stderr)
            continue

        # Write one protein fasta and one nucleotide fasta per HMM profile.
        out_faa = f"{protdb_short}_{hmm_short}.faa"
        extract_hmm_hits_to_faa(hmm_table, protdb, out_faa)

        out_fna = f"{protdb_short}_{hmm_short}.fna"
        extract_hmm_hits_to_fna(hmm_table, nucdb, out_fna, args.len)

    print("Process completed.")
