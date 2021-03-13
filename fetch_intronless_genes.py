#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 12:39:33 2020

@author: andresgarcíagarcía
@edited: gabrielemilioherreraoropeza
@supervisors: katiaaviñapadilla AND maribelhernandezrosales

"""

import sys, os, errno, subprocess

def install(package):
    print(f"Installing {package}...", end = "\t")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    print("DONE")

try:
    import ensembl_rest
except ImportError:
    install("ensembl-rest")
    import ensembl_rest

try:
    import pandas as pd
except ImportError:
    install("pandas")
    import pandas as pd
               

def isFile(string):
    if os.path.isfile(string) == False:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), string)
    else:
        pass
    
def single_exon_isoforms(gene_info):
    return ','.join([transcript['id'] for transcript in gene_info['Transcript']
                             if len(transcript['Exon']) == 1])

# --- File Verification and Import

isFile("fetch_genes.py")

path2intronDB = sys.argv[2]
isFile(path2intronDB)

import fetch_genes

# --- START ---

species = (sys.argv[1]
            if (len(sys.argv) > 1)
            else 'gallus gallus')

print("\nSpecies:", species)

client = ensembl_rest.EnsemblClient()

# --- Fetch information of the s in Gallus gallus
print('\nFetching chromosome data...', end = "   ")

chrom_data, genome_info = fetch_genes.chromosomes_info(species = species)

print("DONE.")

# --- Fetch the ids of all the genes in the Gallus Gallus genome
print('\nFetching all gene IDs...')

all_genes = sum((fetch_genes.genes_in_chrom(species, chrom, chrom_data[chrom]['length']) 
                    for chrom in chrom_data), 
                []
            )


# --- Get the list of exons for all genes
print('\nFetching gene info...')

fields = ['id', 'display_name', 'description', 'biotype', 'start', 'end', 'seq_region_name']

genes = fetch_genes.get_info(all_genes, 
                             (lambda gene_info : {
                                 "exons_n": len({exon['id'] 
                                                 for transcript in gene_info['Transcript']
                                                 for exon in transcript['Exon']}),
                                 "transcripts_n": len(gene_info['Transcript']),
                                 "single_exon_isoforms": single_exon_isoforms(gene_info),
                                 **{field:gene_info[field] for field in fields if field in gene_info}
                                 }
                             )
        )

# --- Filter those which only have one exon and protein-coding ---
one_exon_genes = {gene_id for gene_id, gene_info in genes.items() if gene_info['exons_n'] == 1 and gene_info['biotype'] == 'protein_coding' and gene_info['transcripts_n'] == 1}
multiple_exon_genes = {gene_id for gene_id, gene_info in genes.items() if gene_info['exons_n'] > 1 and gene_info['biotype'] == 'protein_coding'}

# --- Report the results
print(f"\nOf the {len(all_genes)} genes in {species}, "
      f"{len(multiple_exon_genes)} are protein-coding multi-exon genes, "
      f"{len(one_exon_genes)} are protein-coding single-exon genes,")

# --- Identify intronless genes by checking whether single exon genes have introns in the UTR regions

df_introns = pd.read_csv(path2intronDB, sep = '\t')
lst_introns_utr = df_introns['transcript_id'].to_list()
intronless_genes = {gene_id for gene_id in one_exon_genes if not genes[gene_id]['single_exon_isoforms'] in lst_introns_utr}
ui_single_exon_genes = {gene_id for gene_id in one_exon_genes if not gene_id in intronless_genes}

print(f"of which {len(intronless_genes)} are intronless genes,"
      f"and {len(ui_single_exon_genes)} are ui-single-exon genes")

# --- Write the list of the single exon gene IDs
import csv

with open(f"../data/out/{species}-ui_single_exon_genes.tsv", 'w', newline='') as segs_out, \
     open(f"../data/out/{species}-multi_exon_genes.tsv", 'w', newline='') as megs_out, \
     open(f"../data/out/{species}-intronless_genes.tsv", 'w', newline='') as igs_out:

        fieldnames = fields + ['exons_n', 'transcripts_n', 'single_exon_isoforms']
        seg_writer = csv.DictWriter(segs_out, fieldnames=fieldnames, dialect='excel-tab')
        meg_writer = csv.DictWriter(megs_out, fieldnames=fieldnames, dialect='excel-tab')
        ig_writer = csv.DictWriter(igs_out, fieldnames=fieldnames, dialect='excel-tab')

        seg_writer.writeheader()
        meg_writer.writeheader()
        ig_writer.writeheader()
        for gene_id, gene_info in genes.items():
            if gene_id in intronless_genes:
                ig_writer.writerow(gene_info)            
            if gene_id in ui_single_exon_genes:
                seg_writer.writerow(gene_info)
            elif gene_id in multiple_exon_genes:
                meg_writer.writerow(gene_info)