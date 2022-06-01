library(ggplot2)
library(cowplot)
library(gggenes)
args = commandArgs(trailingOnly=TRUE)
# 1 - gfa csv
# 2 - name of block with gene in (for centering)
# 3 - png from bandage
# 4 - output pdf

# Genome blocks approach
genome.blocks <- read.csv(args[1], header=T, stringsAsFactors = F)
#colnames(genome.blocks) <- c("genome", "block", "strand", "start", "end", "colour")
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x, 
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)


KPC.block.locations <- genome.blocks[which(genome.blocks$block==args[2]),c("start", "end")]
rownames(KPC.block.locations) <- genome.blocks[which(genome.blocks$block==args[2]),"genome"]
transformAnnotationsBlocks <- function(genome, default_offset=9926){
  offset <- default_offset-KPC.block.locations[genome, c("start")]
  return(offset)
}

genome.blocks$new.start <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["start"])+transformAnnotationsBlocks(x["genome"]))-10001
genome.blocks$new.end <- apply(genome.blocks, MARGIN=1, function(x) as.numeric(x["end"])+transformAnnotationsBlocks(x["genome"]))-10001

# Order the genomes by similarity?
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)
genome.blocks$genome <- ordered(genome.blocks$genome, 
                                levels=names(genome.paths))



# Subset to unique paths
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]
# And store the number of reps
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) table(genome.paths)[x])
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1, 
                                        "", paste0("n=", x)))


p.blocks <- ggplot(genome.blocks.unique, aes(xmin = new.start, xmax = new.end, forward = forward, y = genome, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"), 
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(0.1, "mm")) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome, labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearized blocks")+
  theme(plot.title=element_text(hjust=0.5))
#graph_file <- system.file("extdata", "../graph-10k-5k.jpg", package = "cowplot")
p.graph <- ggdraw() + draw_image(args[3])+
  ggtitle("Graph representation")+
  theme(plot.title=element_text(hjust=0.5))

pdf(args[4], height=5, width=10)
cowplot::plot_grid(p.blocks, p.graph, nrow=1, align='h')
dev.off()
