#!/usr/bin/env python

import sys
import re

# Ensure proper usage
if len(sys.argv) != 2:
    print("Usage: {} <file>".format(sys.argv[0]))
    sys.exit(1)

file_name = sys.argv[1]

try:
    with open(file_name, 'r') as fh:
        genes = {}
        transcripts = {}
        
        # Process each line in the GFF-like file
        for line in fh:
            if not line.startswith('#'):
                columns = line.strip().split('\t')
                chrom, start, end = columns[0], columns[3], columns[4]

                # Capture gene information
                if columns[2] == 'gene':
                    gene_name, gene_type = None, None
                    match = re.search(r'ID=([^;]+);.*Note=([^;]+)', columns[8])
                    if match:
                        gene_name, gene_type = match.group(1), match.group(2)

                    # Initialize a new gene entry
                    if gene_name not in genes:
                        genes[gene_name] = {
                            "gene_start": start,
                            "gene_end": end,
                            "chr": chrom,
                            "biotype": gene_type,
                            "exon_chromstart": set(),  # Use sets to ensure uniqueness
                            "exon_chromend": set()
                        }

                # Capture mRNA (transcript) information
                elif columns[2] == 'mRNA':
                    transcript_name, gene_name = None, None
                    match = re.search(r'ID=([^;]+);Parent=([^;]+)', columns[8])
                    if match:
                        transcript_name, gene_name = match.group(1), match.group(2)
                        transcripts[transcript_name] = gene_name  # Map transcript to its parent gene

                # Capture exon information associated with transcripts
                elif columns[2] == 'exon':
                    transcript_name = None
                    match = re.search(r'Parent=([^;]+)', columns[8])
                    if match:
                        transcript_name = match.group(1)

                    if transcript_name in transcripts:
                        gene_name = transcripts[transcript_name]
                        if gene_name in genes:
                            # Add exon start and end, ensuring no duplicates
                            genes[gene_name]["exon_chromstart"].add(start)
                            genes[gene_name]["exon_chromend"].add(end)

except FileNotFoundError:
    print("File not found: {}".format(file_name))
    sys.exit(1)

# Output the data in the correct format for topr
print("\t".join(["chrom", "gene_start", "gene_end", "gene_symbol", "biotype", "exon_chromstart", "exon_chromend"]))
for gene, data in genes.items():
    if not gene:
        continue
    # Sort and join exon positions for readability
    exon_starts = ",".join(sorted(data.get("exon_chromstart", []), key=int))
    exon_ends = ",".join(sorted(data.get("exon_chromend", []), key=int))
    print("\t".join([data["chr"], data["gene_start"], data["gene_end"], gene, data["biotype"], exon_starts, exon_ends]))
