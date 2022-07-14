library(tidyverse)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)
library(cowplot)

# Abigail reactor chemical metadata and Accumulibacter clade quantification 


# reactor metadata 

reactor_metadata <- read.csv("metadata/ReactorA_Sample_P_Info.csv")

p_data <- reactor_metadata %>% 
  select(operation_day, P_release, P_final) %>% 
  gather(key="variable", value="value", -operation_day)

p_plot <- p_data %>% 
  ggplot(aes(x=operation_day, y=value)) +
  geom_line(aes(color=variable), size=2) +
  scale_color_manual(values=c("#0FA1D8", "#4D3F83"), labels=c("Total P Remaining \n at end of Aerobic Phase \n", "Total P Release \n at end of Anaerobic Phase")) +
  xlab("Operation Day") +
  ylab("Phosphorus (mg/L) \n") +
  scale_x_continuous(breaks=seq(0,65,5)) +
  theme_classic() +
  theme(legend.title=element_blank(), legend.position=c(.30, .88), axis.title.x=element_blank())
p_plot

# qPCR data 
qpcr_data <- read.csv("data/ppk1-qpcr-results.csv")
qpcr_data$operation_day <- as.factor(qpcr_data$operation_day)

qpcr_plot <- qpcr_data %>% 
  ggplot(aes(x=operation_day, y=copies, fill=clade)) + 
  geom_bar(stat="identity", color="black", size=1, position="dodge") +
  geom_errorbar(aes(ymin = (copies - std_dev), ymax = (copies + std_dev)), width=.3, position=position_dodge(1)) +
  scale_y_log10(limits=c(1,1e6), expand=c(0,0), labels = function(x) format(x, scientific = TRUE)) +
  scale_x_discrete() +
  xlab("Operation Day") +
  ylab("Copies/ng DNA") +
  scale_fill_brewer() +
  theme_classic() +
  theme(legend.position=c("top"))
qpcr_plot

# grids
reactor_grid <- plot_grid(p_plot, qpcr_plot, ncol=2, labels=c("A", "B"))

amplicon_new_grid <- plot_grid(genus_heatmap, ncol=1, labels=c("C"))

abigail_grid <- plot_grid(reactor_grid, amplicon_new_grid, ncol=1, rel_widths=c(1.2,2))
abigail_grid

ggsave("figures/abigail-reactor-grid.png", abigail_grid, width=25, height=15, units=c("cm"))

# sampling dotplot
sampling_metadata <- read.csv("metadata/methods_operation_dates.csv")
sampling_metadata$Method <- factor(sampling_metadata$Method, levels=c("16S_sequencing", "ppk1_qPCR", "metagenomes"))

sampling_chart <- sampling_metadata %>% 
  ggplot(aes(x = Operation_Day, y=fct_rev(Method))) +
  geom_beeswarm(groupOnX=FALSE) +
  geom_point(aes(size=1.2, color=Method)) +
  scale_x_continuous(breaks=seq(0,65,5)) +
  theme_bw() +
  scale_y_discrete(labels=c("Metagenomic \n Sequencing", "ppk1 qPCR", "16S rRNA \n Amplicon Sequencing"))  +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_text(face="bold", size=8), legend.position=c("none")) +
  scale_color_brewer(palette=c("Paired"))

sampling_chart
                     

p_grid <- plot_grid(sampling_chart, p_plot, ncol=1, align="v", rel_heights=c(.8,2), axis="l")
p_grid

p_clade_grid <- plot_grid(p_grid, qpcr_plot, ncol=2, labels=c("A","B"), rel_widths=c(2,1.5), rel_heights=c(2,1.5))
full_grid <- plot_grid(p_clade_grid, genus_heatmap, ncol=1, labels=c("", "C"), rel_heights=c(2,1.5))
full_grid
ggsave("figures/Abigail-metadata-grid.png", full_grid, width=35, height=25, units=c("cm"))


abigail_grid_v1 <- plot_grid(sampling_chart, p_plot, genus_heatmap, ncol=1, labels=c("A", "B", "C"), rel_heights=c(.8,1.8,2.2), align=c("v"), axis="l")
abigail_grid_v1
ggsave("figures/Abigail-metadata-grid-v1.png", abigail_grid_v1, width=25, height=20, units=c("cm"))

# grids with metagenome information as well 
plot_grid(abigail_grid_v1, bins_grid, ncol=2)
