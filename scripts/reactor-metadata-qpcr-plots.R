library(tidyverse)

# Abigail reactor chemical metadata and Accumulibacter clade quantification 


# reactor metadata 

reactor_metadata <- read.csv("metadata/ReactorA_Sample_P_Info.csv")

p_data <- reactor_metadata %>% 
  select(operation_day, P_release, P_final) %>% 
  gather(key="variable", value="value", -operation_day)

p_plot <- p_data %>% 
  ggplot(aes(x=operation_day, y=value)) +
  geom_line(aes(color=variable), size=2) +
  scale_color_manual(values=c("#0FA1D8", "#4D3F83"), labels=c("Total P Remaining at \n end of Aerobic Phase \n", "Total P Release at \n end of Anaerobic Phase")) +
  xlab("Operation Day") +
  ylab("Phosphorus (mg/L)") +
  theme_classic() +
  theme(legend.title=element_blank(), legend.position=c("top"))

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

# grid
reactor_grid <- plot_grid(p_plot, qpcr_plot, ncol=2, labels=c("A", "B"))

amplicon_new_grid <- plot_grid(genus_heatmap, shannon_plot, ncol=1, labels=c("C", "D"))

abigail_grid <- plot_grid(reactor_grid, amplicon_new_grid, ncol=1, rel_heights=c(1.3,2))
abigail_grid

ggsave("figures/abigail-reactor-grid.png", abigail_grid, width=25, height=20, units=c("cm"))
