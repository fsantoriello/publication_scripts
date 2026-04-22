#!/usr/bin/env python3

"""
Created on February 06, 2025, updated with Claude Code Opus 4.7 on April 21, 2026

@author: Francis J Santoriello

This script accepts a csv file of protein accessions. It then:
    1. Retrieves all unique identical protein accessions
    2. Retrieves the associated GenBank contig ids for each unique protein accession
    3. Retrieves the genbank files.
    4. Extracts all CDS features from these files and combines them into a single protein
       fasta file.

This script returns:
    1. all_seqs.fna
    2. all_prots.faa
"""

import os
import time
import argparse
import textwrap
import csv
from Bio import SeqIO
from Bio import Entrez

cwd = os.getcwd()

# NCBI rate limits: 3 req/s without api key, 10 req/s with. Throttle between
# per-protein IPG fetches so large input lists don't get rate-limited.
IPG_THROTTLE_SECONDS = 0.4

### User Input ###
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(
                prog='linearphage_search_ncbi.py',
                formatter_class=MyFormatter,
                description=textwrap.dedent('''\
                    linearphage_search_ncbi.py

                    Created on February 06, 2025

                    @author: Francis J Santoriello

                    This script accepts a csv file of protein accessions. It then:
                       1. Retrieves all unique identical protein accessions
                       2. Retrieves the associated GenBank contig ids for each unique
                          protein accession
                       3. Retrieves the genbank files.
                       4. Extracts all CDS features from these files and combines them
                          into a single protein fasta file.

                    This script returns:
                       1. all_seqs.fna
                       2. all_prots.faa
                    '''),
                usage='%(prog)s [options] <prot_id_csv_file> <email>\n\n '
                      'Dependencies: python 3, biopython\n')
parser.add_argument('-l', '--len', type=int, default=15000,
                    help='Minimum sequence length for filtering')
parser.add_argument('-o', '--outdir', default=cwd,
                    help='Output path for generated files')
parser.add_argument('--api',
                    help='NCBI API Key for increased query rate limit')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Print raw IPG records as they are fetched')
parser.add_argument('filename',
                    help='CSV file containing list of protein ids')
parser.add_argument('email',
                    help='Email address for use with Entrez')


args = parser.parse_args()

### Variables ###
Entrez.email = args.email
if args.api is not None:
    Entrez.api_key = args.api

os.makedirs(args.outdir, exist_ok=True)

# Read all rows of the CSV and flatten into a single list of accessions, so
# the input can be one-per-line or a single comma-separated row.
with open(args.filename, newline='') as f:
    reader = csv.reader(f)
    protein_accessions = [acc.strip() for row in reader for acc in row if acc.strip()]
if not protein_accessions:
    raise SystemExit(f"No protein accessions found in {args.filename}")

all_prots = os.path.join(args.outdir, "all_prots.faa")
all_seqs = os.path.join(args.outdir, "all_seqs.fna")

### Modules ###
def all_identical_protein_ids_to_faa(input_protein_ids):
    """
    Takes a list of protein IDs, retrieves all unique identical protein IDs in IPG,
    fetches the associated GenBank nucleotide record for each, extracts every CDS
    feature, and writes the combined nucleotide and protein fasta files.
    Args:
       input_protein_ids: A list of protein IDs.
    """
    # Maps nucleotide accession -> organism label, populated while parsing IPG records
    # and reused when writing FASTA headers.
    org_dict = {}

    # Retrieve the Identical Protein Group (IPG) record for each user-entered
    # protein accession and parse out one nucleotide accession per unique source.
    for prot_id in input_protein_ids:
        ipg_record = Entrez.efetch(id=prot_id, db='ipg').read().decode()
        if args.verbose:
            print(ipg_record)

        # IPG records are TSV; skip the header row then split each row on tabs.
        nt_acc_clean = [row.split("\t") for row in ipg_record.split("\n")[1:] if row]

        # Precompute the set of assembly accessions already covered by an INSDC
        # row (column 9). Used to skip RefSeq rows that duplicate an INSDC row.
        insdc_assemblies = {row[9] for row in nt_acc_clean
                            if len(row) > 9 and row[1] == 'INSDC'}

        # Prefer INSDC rows; keep RefSeq rows only when not duplicated by INSDC.
        for acc in nt_acc_clean:
            if len(acc) <= 9:
                continue
            if acc[1] == 'INSDC':
                if acc[2]:
                    org_dict[acc[2]] = f"{acc[8]} {acc[9]}"
                    print(f"{acc[2]} : {acc[8]} {acc[9]} appended.")
                else:
                    print("No associated accession id.")
            elif acc[1] == 'RefSeq':
                if acc[9] in insdc_assemblies:
                    print(f"{acc[9]} passed as duplicate.")
                elif acc[2]:
                    org_dict[acc[2]] = f"{acc[8]} {acc[9]}"
                    print(f"{acc[2]} : {acc[8]} {acc[9]} appended.")
                else:
                    print("No associated accession id.")

        # Stay under NCBI's request/sec limit on large input lists.
        time.sleep(IPG_THROTTLE_SECONDS)

    print(f"Retrieved {len(org_dict)} nucleotide accession ids.")

    # Retry settings for transient NCBI failures during chunk fetch/parse.
    max_retries = 3
    retry_delay = 2  # seconds

    with open(all_seqs, 'w') as seq_file, open(all_prots, 'w') as prot_file:
        # Break nt accession ids into chunks of 9999 - max allowable uids by efetch.
        nt_accs = list(org_dict)
        chunked_nt_accs = [nt_accs[i:i + 9999] for i in range(0, len(nt_accs), 9999)]

        print("Fetching genbank files in chunks...")

        for i, chunk in enumerate(chunked_nt_accs):
            for attempt in range(1, max_retries + 1):
                handle = None
                try:
                    handle = Entrez.efetch(id=chunk, db='nuccore',
                                           rettype='gbwithparts', retmode='text')
                    print(f"Chunk {i+1}: fetching {len(chunk)} sequences")

                    processed = 0
                    for record in SeqIO.parse(handle, "genbank"):
                        if len(record.seq) < args.len:
                            print(f"{record.id} sequence too short.")
                            continue

                        # Build everything in memory first - nothing touches disk
                        # until the record is fully parsed, so an exception can't
                        # desync the nucleotide and protein files.
                        org_label = org_dict.get(record.id, "unknown")
                        nt_entry = f">{record.id} [{org_label}]\n{record.seq}\n"

                        prot_entries = []
                        feat_total = 0
                        for feat in record.features:
                            feat_total += 1
                            if feat.type != "CDS":
                                continue
                            # Require all three qualifiers to build a useful header.
                            if not all(q in feat.qualifiers for q in
                                       ('protein_id', 'product', 'translation')):
                                print("Passed for missing qualifiers.")
                                continue
                            # '?' in the translation indicates an ambiguous residue
                            # from an undefined codon - skip these entries.
                            if '?' in feat.qualifiers['translation'][0]:
                                print(f"{feat.qualifiers['protein_id'][0]} passed for "
                                      f"undefined translation.")
                                continue
                            header = (f">{feat.qualifiers['protein_id'][0]} "
                                      f"{feat.qualifiers['product'][0]} "
                                      f"[{record.id} | {org_label}]")
                            prot_entries.append(
                                f"{header}\n{feat.qualifiers['translation'][0]}\n")

                        # Commit both files only after the record is fully parsed.
                        seq_file.write(nt_entry)
                        prot_file.writelines(prot_entries)

                        processed += 1
                        remaining = len(chunk) - processed
                        print(f"Chunk {i+1}: appended {len(prot_entries)} CDS of "
                              f"{feat_total} features from {record.id}. "
                              f"({remaining} remaining).")
                    break  # chunk succeeded, exit retry loop

                except Exception as e:
                    print(f"Attempt {attempt} to retrieve or parse genbank data "
                          f"failed: {e}")
                    if attempt < max_retries:
                        print(f"Retrying in {retry_delay} seconds...\n")
                        time.sleep(retry_delay)
                    else:
                        print(f"Chunk {i+1}: all attempts failed, skipping.")
                finally:
                    if handle is not None:
                        handle.close()

        print("All nucleotide and all protein fasta files complete.")


### Main ###

if __name__ == "__main__":
    all_identical_protein_ids_to_faa(protein_accessions)
    print("Process completed.")
