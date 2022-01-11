library(tidyverse)

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

left_join(ab_genome_info_filtered, abigail_metadata) %>% 
  separate(taxonomy, into=c("taxonomy", "phylum"), sep="p__") %>%
  mutate(phylum = gsub(";.*", "", phylum)) %>%
  ggplot(aes(x=r2_mean, y=nucl_diversity)) + geom_point(aes(color=phylum), size=2.5) + facet_wrap(~ sample) + theme_bw()
