library(dada2)
library(phyloseq)
library(ape)
library(ggplot2)
library(ampvis2)
library(grid)
library(gridExtra)
library(vegan)
library(tidyverse)
library(patchwork)
library(cowplot)

# before starting the dada2 workflow, samples must be demultiplexed, primers/adapters are removed, and the F and R files contain reads in matching order
# this preprocessing script works through a time-series of samples from engineered bioreactors, amplified the 16S region using the V3-V4 primers/region, and was sequenced with the Illumina 2x300 chemistry (incorrectly, was supposed to be 2x250, but will be more to throw out)
# forward primer: CCTACGGGNGGCWGCAG
# reverse primer: GACTACHVGGGTATCTAATCC
# also check that all fastq files and databases aren't in the cloud, will need to pull them down

# primer lengths
primerF <- "CCTACGGGNGGCWGCAG"
primerR <- "GACTACHVGGGTATCTAATCC"
lengthF <- nchar(primerF)
lengthR <- nchar(primerR)

# path and files setup
path <- "data/16S_tags/rawData"
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names=TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# inspect quality profiles of a couple forward/reverse reads
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# filter and quality trim
filtFs <- file.path("data/16S_tags/cleanData", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("data/16S_tags/cleanData", paste0(sample.names,"_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,215), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, trimLeft = c(lengthF, lengthR))

# check quality again
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

# error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# sample inference for obtaining sequence variants
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# can use the dada function to pool sequences to inform from multiple samples, but can address that at a later time

# merge paired reads to get the full denoised dataasets
# most reads should merge, if not getting a lot of merge paired-reads check trimming parameters
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
# can get specific length of expected amplicon length
# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# looking at number of reads that went through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# taxa assigned with silva, general and exact species
taxa <- assignTaxonomy(seqtab.nochim, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/R1PopDynamics/databases/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
taxa <- addSpecies(taxa, "/Users/emcdaniel/Desktop/McMahon-Lab/EBPR-Projects/R1R2/R1PopDynamics/databases/silva_species_assignment_v138.1.fa.gz")
# taxa assigned with GTDB taxonomy, general and exact species
# gtdb_taxa <- assignTaxonomy(seqtab.nochim, "databases/GTDB_bac-arc_ssu_r86.fa.gz", multithread=FALSE)
# gtdb_taxa <- addSpecies(gtdb_taxa, "databases/GTDB_dada2_assignment_species.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
taxa.print

####################################
# create a phyloseq object
####################################
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out) <- c("timepoint")
rownames(samples.out) <- rownames(seqtab.nochim)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samples.out), tax_table(taxa))

# merge with sample information
sample_info <- read.csv("metadata/ReactorA_Sample_P_Info.csv")
colnames(sample_info)[1] <- "timepoint"
sample_metadata <- left_join(samples.out, sample_info)
row.names(sample_metadata) <- sample_metadata$timepoint

# phyloseq object with metadata
ps2 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(sample_metadata), tax_table(taxa))

# prune phyloseq object for just Accumulibacter and plot ASVs 
acc_ps <- subset_taxa(ps2, Genus=="Candidatus Accumulibacter")
acc_tree <- rtree(ntaxa(acc_ps), rooted=TRUE, tip.label=taxa_names(acc_ps))
plot(acc_tree)

# ampvis
#source the phyloseq_to_ampvis2() function from the gist
#devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")
ampvis2_obj <- phyloseq_to_ampvis2(ps2)
accumulibacter_OTUs <- amp_subset_taxa(ampvis2_obj, tax_vector=c("Candidatus Accumulibacter"), normalise=TRUE)

# plots of 16S dynamics over time-series and abundance
genus_heatmap <- amp_heatmap(ampvis2_obj, tax_aggregate = "Genus", group_by="operation_day", tax_show=12, plot_values = FALSE, plot_colorscale = "sqrt", tax_add="Phylum", plot_legendbreaks=c(0,5,10,20, 30, 40)) + theme(axis.text.y=element_text(face=c("italic"), size=5), legend.position=c("right"), axis.text.x=element_text(angle=0), axis.ticks.x = element_blank()) + scale_x_discrete(position="bottom")

genus_heatmap

genus_boxplot <- amp_boxplot(ampvis2_obj, tax_show = 12, tax_add="Phylum")

acc_ASVs_heatmap <- amp_heatmap(accumulibacter_OTUs, tax_aggregate="OTU", tax_show=5, plot_values=FALSE, group_by="timepoint", plot_colorscale = "sqrt", normalise=FALSE, plot_legendbreaks=c(0,5,10,15,20,25))+ scale_y_discrete(labels=c("C. Accumulibacter ASV5", "C. Accumulibacter ASV4", "C. Accumulibacter ASV3", "C. Accumulibacter ASV2", "C. Accumulibacter ASV1")) + theme(axis.text.y = element_text(face="italic", size=8))

genus_heatmap
acc_ASVs_heatmap

# shannon diversity 
shannon_plot <- plot_richness(ps2, x="operation_day", measure="Shannon") + scale_x_continuous(expand=c(0,0), limits=c(0,62), breaks=seq(0,62,2)) + ylab('Shannon\n Alpha Diversity\n') + xlab("Operation Day") + theme_bw() + theme(axis.title.x=element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=7), strip.background=element_blank(), strip.text.x=element_blank(), plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))

shannon_plot

plot_richness(ps2, x="sample", measure="Simpson") + ylab('Shannon\n Alpha Diversity\n') + theme_bw() + theme(axis.title.x=element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=7), strip.background=element_blank(), strip.text.x=element_blank(), plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))

plot_richness(ps2, x="sample", measure="Fisher") + ylab('Shannon\n Alpha Diversity\n') + theme_bw() + theme(axis.title.x=element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=7), strip.background=element_blank(), strip.text.x=element_blank(), plot.title=element_text(size=12, face="bold"), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6))

amplicon_grid <- plot_grid(genus_heatmap, shannon_plot, ncol = 1, align=c("v"), axis = "l", labels=c("A","B"))

amplicon_grid

ggsave("figures/Abigail-genus-heatmap.png", genus_heatmap, width=10, height=5, units=c("in"))
ggsave("figures/Abigail-genus-boxplot.png", genus_boxplot, width=9, height=4, units=c("in"))
ggsave("figures/Abigail-Acc-ASVs.png", acc_ASVs_heatmap, width=10, height=5, units=c("in"))
