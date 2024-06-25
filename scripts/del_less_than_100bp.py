#!/bin/python

import sys

def filter_genes(input_file, output_file, min_length=100):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        write_gene = False
        gene_sequence = ""
        for line in infile:
            if line.startswith(">"):
                if write_gene and len(gene_sequence) >= min_length:
                    outfile.write(gene_header)
                    outfile.write(gene_sequence + "\n")
                gene_header = line
                gene_sequence = ""
                write_gene = True
            else:
                gene_sequence += line.strip()
        # Check the last gene
        if write_gene and len(gene_sequence) >= min_length:
            outfile.write(gene_header)
            outfile.write(gene_sequence + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python XXX.py <input gene catalog> <output gene catalog>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    filter_genes(input_file, output_file)

