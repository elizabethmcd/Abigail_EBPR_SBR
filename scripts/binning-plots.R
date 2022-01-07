library(tidyverse)

checkm_stats <- read_tsv("results/binning/final-bins-checkm-stats.tsv", col_names = FALSE)
colnames(checkm_stats) <- c("bin", "lineage", "completeness", "contamination", "size", "contigs", "gc")
gtdb_tax <- read_tsv("results/binning/new-bins-GTDB-tax.tsv", col_names = FALSE)
colnames(gtdb_tax) <- c("bin", "taxonomy")

bin_table <- left_join(checkm_stats, gtdb_tax) %>% 
  select(bin, taxonomy, completeness, contamination, size, contigs, gc)

write.csv(bin_table, "results/binning/Abigail-bins-info.csv", quote = FALSE, row.names = FALSE)

# Metadata and relative abundance
metadata <- read.csv("results/binning/Abigail-bins-info.csv")

relative_abundance <- read_tsv("results/binning/abigail-relative-abundance.txt") %>% 
  filter(Genome != "unmapped")
colnames(relative_abundance) <- c("bin", "A-2021-03-17", "A-2021-04-03")

abund_metadata <- left_join(metadata, relative_abundance) %>% 
  select(bin, group, `A-2021-03-17`, `A-2021-04-03`) %>% 
  pivot_longer(!c(bin, group), names_to="date", values_to="abundance")

abund_avg <- left_join(metadata, relative_abundance) %>% 
  select(bin, group, `A-2021-03-17`, `A-2021-04-03`) %>%
  mutate(avg = (`A-2021-03-17` + `A-2021-04-03`) / 2 ) %>% 
  select(bin, group, avg)

bin_table <- left_join(metadata, relative_abundance)
bin_table$avg <- (bin_table$`A-2021-03-17` + bin_table$`A-2021-04-03`) / 2

write.csv(bin_table, "results/binning/Abigail-bins-table.csv", quote=FALSE, row.names = FALSE)

bar_plot <- abund_metadata %>% ggplot(aes(x=factor(date), y=abundance, fill=group)) + geom_bar(stat="identity", color="black") + scale_fill_brewer(palette = "Set3") + scale_y_continuous(expand=c(0,0), breaks=seq(0, 45, 5)) + ylab("Relative Abundance") + xlab("Sample Date") + theme_bw()

avg_plot <- abund_avg %>% ggplot(aes(x=reorder(bin, -avg), y=avg, fill=group)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Set3") + ylab("Average Relative Abundance") + xlab("Genome") + scale_y_continuous(expand=c(0,0), breaks=seq(0, 12, 2)) + theme_classic() + theme(axis.text.x=element_blank())

ggsave("figures/abigail-bins-sample-barplot.png", bar_plot, width=7, height=9, units=c("in"))
ggsave("figures/abigail-bins-avg-rel-abund.png", avg_plot, width=11, height=5, units=c("in"))
