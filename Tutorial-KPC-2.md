# KPC-2

In this tutorial we'll use `pangraph` to explore a small dataset of contigs which all contain a beta-lactamase gene KPC-2 (or KPC-3 - there is 1 SNP difference).

*Dependencies:* python (biopython), pangraph, R (ggplot2, gggenes, cowplot), blast

## Dataset

The contigs are in `data/kpc2-contigs.fa`. These have been extracted from the dataset detailed in David et al. - see `data/kpc/Table-S2-David-EUSCAPE.csv`.

The names of the contigs give the `contig|study_ID` e.g.

```
CABGBR010000009.1|EuSCAPE_ES089
```

We could run pangraph on the full contigs (n=44). However, first we extract a consistent region around the KPC gene. We extract 10kb upstream and 5kb downstream and only keep those contigs that are long enough for this purpose.

```
python scripts/extract-region.py --gene data/kpc/kpc2.fa \
				--input data/kpc/kpc-contigs.fa \
				--output data/kpc/kpc-contigs-u10k-d5k.fa \
				--upstream 10000 \
				--downstream 5000 \
				--complete
```

We can now run pangraph on these extracted regions (n=34). This should take <1 minute to run.

```
pangraph build data/kpc/kpc2-contigs-u10k-d5k.fa > data/kpc/pangraph_kpc_u10k_d5k.json
pangraph export --edge-minimum-length 0 data/kpc/pangraph_kpc_u10k_d5k.json -p pangraph_kpc_u10k_d5k -o data/kpc/
```

Pangraph gives us three forms of output:

* `data/kpc/pangraph_kpc_u10k_d5k.json` - json file storing the whole graph
* `data/kpc/pangraph_kpc_u10k_d5k.gfa` - graph in [GFA](http://gfa-spec.github.io/GFA-spec/GFA1.html) format
* `data/kpc/pangraph_kpc_u10k_d5k.fa` - multifasta containing the consensus sequences of the pangenome blocks


We first want to find the block that contains the KPC gene.

```
makeblastdb -in data/kpc/pangraph_kpc_u10k_d5k.fa -dbtype 'nucl'
geneBlock=$(blastn -query data/kpc/kpc2.fa -db data/kpc/pangraph_kpc_u10k_d5k.fa -outfmt 6 | cut -f 2)
```

In the pre-computed data, the gene is on the block with ID `HALKKORSXJ` (this is just a random string).

We then convert the gfa into a csv that stores each contig as a linear path of blocks with their start and end positions.

```
python scripts/prepare-pangraph-gfa.py data/kpc/pangraph_kpc_u10k_d5k.gfa
```

This makes three output files in `data/kpc`

* `${input}.blocks.csv` - dataset of genome and block start/end positions
* `${input}.colours.csv` - blocks with colours (hex codes). Blocks that only appear once in the dataset are coloured grey, others get random colours.
* `${input}.coloured.gfa` - a gfa with the block colours added as an extra field

We can then plot a graph with Bandage. Using the `--colour custom` flag will use the block colours we have added.

```
Bandage image data/kpc/pangraph_kpc_u10k_d5k.gfa.coloured.gfa data/kpc/pangraph_kpc_u10k_d5k.gfa.png --height 400 --width 7000 --colour custom
```

We can now plot the blocks in a linear fashion with an R script, using the same colours.

```
Rscript scripts/plot-blocks.R data/kpc/pangraph_kpc_u10k_d5k.gfa.blocks.csv $geneBlock data/kpc/pangraph_kpc_u10k_d5k.gfa.png data/kpc/pangraph_kpc_plot.pdf
```

The script `make-plot-kpc.sh` puts all these together:

```
./scripts/make-plot-kpc.sh data/kpc/pangraph_kpc_u10k_d5k
```

N.B. Because the block colours are random, there is a specified seed in the `readGFA` function of `scripts/prepare-pangraph-gfa.py` to make this reproducible. You can alter that seed if the colours aren't working for you.
