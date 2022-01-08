library(tidyverse)

profiles_path <- "results/inStrain/quick_profiles/"
files <- dir(profiles_path, pattern='*.csv')
abigail_profiles <- data_frame(filename = files) %>% 
  mutate(file_contents=map(filename, ~ read.csv(file.path(profiles_path, .)))
  ) %>% 
  unnest()

abigail_profiles_table <- abigail_profiles %>% 
  filter(breadth > 0.9 & coverage > 10) %>% 
  mutate(sample = gsub("-genomeCoverage.csv", "", filename)) %>% 
  select(genome, sample, coverage, breadth)

abigail_metadata <- read.csv("results/binning/Abigail-bins-table.csv") %>% 
  select(bin, taxonomy, completeness, contamination)
colnames(abigail_metadata) <- c("genome", "taxonomy", "completeness", "contamination")

abigail_profiles_info <- left_join(abigail_profiles_table, abigail_metadata) %>% 
  filter(completeness > 80 & contamination < 10)

abigail_pairings <- abigail_profiles_info %>% 
  select(genome, sample) %>% 
  mutate(genome = paste(genome, ".fa", sep="")) %>% 
  mutate(sample = paste(sample, "-spRep.sorted.bam", sep=""))

write_tsv(abigail_pairings, "metadata/Abigail-inStrain-queues.txt", col_names = FALSE)
