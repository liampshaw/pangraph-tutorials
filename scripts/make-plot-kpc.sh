#!/bin/bash
baseFile=$1

# Find KPC containing block
echo "Running blast to find KPC-2..."
makeblastdb -in "$baseFile".fa -dbtype 'nucl'
geneBlock=$(blastn -query data/kpc/kpc2.fa -db "$baseFile".fa -outfmt 6 | cut -f 2)

# Prepare files
echo "Preparing files..."
python scripts/prepare-pangraph-gfa.py "$baseFile".gfa

# Bandage plot
echo "Plotting Bandage graph plot..."
/Applications/Bandage.app/Contents/MacOS/Bandage image "$baseFile".gfa.coloured.gfa "$baseFile".gfa.png --height 3000 --width 5000 --colour custom

# R plot
echo "Plotting linearized blocks..."
Rscript scripts/plot-blocks.R "$baseFile".gfa.blocks.csv $geneBlock "$baseFile".gfa.png "$baseFile".blocks.pdf
