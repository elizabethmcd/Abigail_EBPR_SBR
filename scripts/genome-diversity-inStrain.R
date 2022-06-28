library(tidyverse)
library(patchwork)
library(cowplot)
library(ggpubr)

# Abigail MAGs inStrain genome info - diversity within samples in the time-series 

genome_info_path <- "results/inStrain/genome_info/"
files <- dir(genome_info_path, pattern="*_genome_info.tsv")
ab_genome_info <- data_frame(filename = files) %>%
  mutate(file_contents = map(filename, ~ read_tsv(file.path(genome_info_path, .)))
  ) %>%
  unnest()

# filter
ab_genome_info_filtered <- ab_genome_info %>% 
  filter(coverage > 10 & breadth > 0.9) %>% 
  separate(filename, into=c("org", "sample"), sep="-vs-") %>% 
  select(genome, sample, coverage, breadth, nucl_diversity, r2_mean, d_prime_mean)
ab_genome_info_filtered$sample <- gsub(".IS_genome_info.tsv", "", ab_genome_info_filtered$sample)

# merge with metadata info

div_plot <- left_join(ab_genome_info_filtered, abigail_metadata) %>% 
  separate(taxonomy, into=c("taxonomy", "phylum"), sep="p__") %>%
  mutate(phylum = gsub(";.*", "", phylum)) %>%
  ggplot(aes(x=r2_mean, y=nucl_diversity)) + geom_point(aes(color=phylum), size=2.5) + facet_wrap(~ sample) + theme_bw() + ylab("Nucleotide Diversity π") + xlab(expression(r ^2)) + labs(color=c("Phylum")) + theme(legend.position = c("bottom"))


# merge with metadata info for group name 
colnames(metadata)[1] <- c("genome")
div_metadata <- left_join(ab_genome_info_filtered, metadata) %>% 
  select(genome, sample, nucl_diversity, r2_mean, group)

# plot to match colors of the bar plot of relative abundance 
facet_labels <- c("Operation Day 43", "Operation Day 60")
names(facet_labels) <- c("Abigail-2021-03-17", "Abigail-2021-04-03")

div_plot <- div_metadata %>% 
  ggplot(aes(x=r2_mean, y=nucl_diversity)) +
  geom_point(aes(color=group), size=3) + 
  scale_color_brewer(palette = "Set3") +
  facet_wrap(~ sample, labeller = labeller(sample = facet_labels)) +
  theme_bw() +
  ylab("Nucleotide Diversity π") +
  xlab(expression(r ^2)) +
  labs(color=c("Lineage")) 

label1.1 <- data.frame(r2_mean=0.5, nucl_diversity=.018, lab='Clade IC', sample = factor("Abigail-2021-03-17", levels=c("Abigail-2021-03-17", "Abigail-2021-04-03")))

label1.2 <- data.frame(r2_mean=0.5, nucl_diversity=.0175, lab='Clade IC', sample = factor("Abigail-2021-04-03", levels=c("Abigail-2021-03-17", "Abigail-2021-04-03")))

div_plot_labels <- div_plot + 
  geom_label(data = label1.1, label="Clade IC") +
  geom_label(data = label1.2, label="Clade IC")

# Grid of 16S data, relative abundance, and diversity 
bins_grid <- plot_grid(bar_plot, div_plot_labels, ncol = 2, align="h", axis="b", labels=c("A", "B"))
bins_grid

abigail_bins_grid <- ggarrange(bar_plot, div_plot_labels, labels = c("A", "B"), common.legend=TRUE, legend = "bottom")
abigail_bins_grid

abigail_grid <- plot_grid(amplicon_grid, bins_grid, ncol=1)
abigail_grid

ggsave("figures/abigail-bins-div-grid.png", abigail_bins_grid, width=30, height=20, units=c("cm"))

# combine with qPCR data and Accumulibacter ASVs