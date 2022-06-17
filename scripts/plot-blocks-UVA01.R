require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(vegan)
require(ggrepel)

args = commandArgs(trailingOnly=TRUE)
# Arguments: 
# 1 - gfa csv
# 2 - png from bandage
# 3 - output pdf

# Genome blocks approach
genome.blocks <- read.csv(args[1], header=T, stringsAsFactors = F)
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)

block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block)

# Easier labels for blocks
block.names <- unique(genome.blocks$block)
new.block.names <- letters[1:length(block.names)]
names(new.block.names) <- block.names
genome.blocks$block.simple <- new.block.names[genome.blocks$block]

# Choose a random anchor block present in all plasmids at copy of one
genomes.per.block <- sapply(unique(genome.blocks$block), function(x) length(unique(genome.blocks$genome[which(genome.blocks$block==x)])))
genomes.per.block <- genomes.per.block[names(block.counts)]
anchor.block <- names(genomes.per.block)[which(genomes.per.block==length(unique(genome.blocks$genome)))][1]

anchor.block.locations <- genome.blocks[which(genome.blocks$block==anchor.block),c("start", "end")]
rownames(anchor.block.locations) <- genome.blocks[which(genome.blocks$block==anchor.block),"genome"]
anchor.block.locations$total.length <- sapply(rownames(anchor.block.locations), function(x) max(genome.blocks[which(genome.blocks$genome==x),"end"]))
transformAnnotationsBlocks <- function(genome, default_offset=0){
  offset <- (default_offset-anchor.block.locations[genome, c("start")]) 
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))

genome.blocks$genome <- ordered(genome.blocks$genome,
                                levels=sort(as.character(unique(genome.blocks$genome)), decreasing = TRUE))

genome.blocks$centre <- (genome.blocks$new.start+genome.blocks$new.end)/2

p.blocks <- ggplot(genome.blocks, aes(xmin = new.start, xmax = new.end, forward = forward, y = as.factor(genome), fill = block)) +
  geom_gene_arrow(arrow_body_height = unit(5, "mm"), 
                  arrowhead_height = unit(5, "mm"),
                  arrowhead_width = unit(0, "mm"), size=0) +
  geom_text_repel(aes(x = centre, y = genome, label = block.simple),
                  nudge_y = 0.5,size=2)+
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  ggtitle("Linear representation")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")
p.graph <- ggdraw() + draw_image(args[2])+
  ggtitle("Graph representation")+
  theme(plot.title=element_text(hjust=0.5))

pdf(args[3], height=5, width=10)
cowplot::plot_grid(p.blocks, p.graph, nrow=1, align='h')
dev.off()
