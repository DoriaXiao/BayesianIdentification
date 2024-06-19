setwd("Simulation_stucy/Sim_result")

library(readxl)
library(ggplot2)
library(patchwork)

Sim_results <- read_excel("Sim_results.xlsx", 
                          sheet = "Sheet4")

colnames(Sim_results)
Sim_results$Class_probabilities <- factor(Sim_results$Class_probabilities, levels=c('D2','D6','D10'),
                                          labels=c('D2','D6','D10'))

Sim_results$Standard_deviations <- factor(Sim_results$Standard_deviations, levels=c('C5','N100','N50','N5','N1'),
                                          labels=c('C5','N100','N50','N5','N1'))
Sim_results$`Rhat>1.10` <- as.numeric(Sim_results$`Rhat>1.10`)
Sim_results$Stuck_sequence<- as.numeric(Sim_results$Stuck_sequence)
Sim_results$Stuck_chain <- as.numeric(Sim_results$Stuck_chain)
Sim_results$`Miniscule-class` <- as.numeric(Sim_results$`Miniscule-class`)
Sim_results$Total <- as.numeric(Sim_results$Total)

plot1 <- ggplot(Sim_results) + 
  aes(x = Standard_deviations, group = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = " ",
       x = "Standard deviations",
       color = "Problematic behavior") + 
  geom_line(aes(y = `Rhat>1.10`, color =  "Rhat>1.10",
                linetype = "Rhat>1.10"), size = 1, alpha = 1) +
  
  geom_line(aes(y = Total, color =  "Stuck sequence and/or miniscule class",
                linetype = "Stuck sequence and/or miniscule class"), size = 1, alpha = 1) +
  scale_color_manual(
    values = c("Rhat>1.10" = "black", 
               "Stuck sequence and/or miniscule class" = "black")
  ) +
  
  scale_linetype_manual(
    values = c("Rhat>1.10" = "dashed", 
               "Stuck sequence and/or miniscule class" = "solid"),
    guide =  FALSE
  ) + 
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",             # Adjust legend position
    panel.grid.major.y = element_blank(),  # Remove major grid lines on y-axis
    panel.grid.minor.y = element_blank(),  # Remove minor grid lines on y-axis
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(linewidth = 0.5),  # Adjust axis line size
    strip.background = element_blank(),    # Remove strip background
    panel.spacing = unit(0, "lines"),      # Remove spacing between facets
    
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  ) +
  facet_grid(~Class_probabilities, scales = "free_y") +
  guides(color =  guide_legend(title = "Problematic behavior", 
                               override.aes = list(linetype = c("dashed", "solid"))))+ 
  my_custom_theme()


plot2 <- ggplot(Sim_results) + 
  aes(x = Standard_deviations, group = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = "Percentage of problematic behavior",
       x = "Standard deviations",
       color = "Problematic behavior") + 
  geom_line(aes(y = Stuck_sequence, color = "Stuck sequence (20+ iterations)", 
                linetype = "Stuck sequence (20+ iterations)"), size = 1, alpha = 1) +
  geom_line(aes(y = Stuck_chain, color = "Stuck chain",  
                linetype = "Stuck chain"), size = 1, alpha = 1) +
  geom_line(aes(y = `Miniscule-class`, color = "Miniscule class", 
                linetype = "Miniscule class"), size = 1,  alpha = 1) +
  scale_color_manual(
    values = c("Stuck sequence (20+ iterations)" = "black", 
               "Stuck chain" = "black",
               "Miniscule class" = "black")
  ) +
  scale_linetype_manual(
    values = c("Stuck sequence (20+ iterations)" = "solid",
               "Stuck chain" = "dashed", 
               "Miniscule class" = "dotdash"),
    guide =  FALSE
  ) + 
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",             # Adjust legend position
    panel.grid.major.y = element_blank(),  # Remove major grid lines on y-axis
    panel.grid.minor.y = element_blank(),  # Remove minor grid lines on y-axis
    panel.border = element_blank(),        # Remove panel border
    axis.line = element_line(linewidth = 0.5),  # Adjust axis line size
    strip.text = element_blank(),
    strip.background = element_blank(),    # Remove strip background
    panel.spacing = unit(0, "lines"),      # Remove spacing between facets
    
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    axis.title.y = element_text(size = 14, hjust = - 20)
  ) +
  facet_grid(~Class_probabilities, scales = "free_y") +
  guides(color = guide_legend(title = "Problematic behavior", 
                              override.aes = list(linetype = c("solid", "dashed", "dotdash")))) + 
  my_custom_theme()





plot1 + plot2 + plot_layout(ncol = 1)
######################################################################################
Sim_results <- read_excel("Sim_results.xlsx", 
                          sheet = "Sheet3")


colnames(Sim_results)



Sim_results$Class_probabilities <- factor(Sim_results$Class_probabilities, levels=c('D2','D6','D10'),
                                          labels=c('D2','D6','D10'))

Sim_results$Standard_deviations <- factor(Sim_results$Standard_deviations, levels=c('C5','N100','N50','N5','N1'),
                                          labels=c('C5','N100','N50','N5','N1'))
Sim_results$Smaller_step_sizes <- factor(Sim_results$Smaller_step_sizes, levels=c('0','1'),
                                         labels=c('Default step size', 'Smaller step size'))


Sim_results$`Rhat>1.10` <- as.numeric(Sim_results$`Rhat>1.10`)
Sim_results$`Stuck-chain_20` <- as.numeric(Sim_results$`Stuck-chain_20`)
Sim_results$Peresist_stuck_chain <- as.numeric(Sim_results$Peresist_stuck_chain)
Sim_results$`Miniscule-class` <- as.numeric(Sim_results$`Miniscule-class`)
Sim_results$either_happen <- as.numeric(Sim_results$either_happen)

plot1 <- ggplot(Sim_results) + 
  aes(x = Standard_deviations, y = `Rhat>1.10`, color = Class_probabilities, group = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = expression("Proportion of 4-chain batches with " * hat(R) > 1.10), x = " ",
       color = "Class probability priors") + 
  geom_line(linewidth = 1, alpha = 0.5) +  # Adjust size and alpha for visibility
  theme_minimal(base_size = 14) +  # Adjust base font size
  theme(legend.position = "right",  # Left-align the legend to the top-left corner
        legend.justification = c(0, 1),  # Left-align the legend
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),      # Remove panel border
        axis.line = element_line(linewidth = 0.5)) +  # Adjust axis line size
  facet_grid(~Smaller_step_sizes)





plot2 <- ggplot(Sim_results) + 
  aes(x = Standard_deviations, group = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = "Percentage of chains with problematic behavior(s) ", x = "Standard Deviations",
       color = "Problematic behavior") + 
  geom_line(aes(y = `Stuck-chain_20`, color = "Chain stuck (20+ iterations)", linetype ="Chain stuck (20+ iterations)"), size = 1, alpha = 0.5) +
  geom_line(aes(y = Peresist_stuck_chain, color = "Persistent stuck chain",  linetype = "Persistent stuck chain"), size = 1, alpha = 0.5) +
  geom_line(aes(y = `Miniscule-class`, color = "Miniscule-class (excluding Persistent stuck) chain",
                linetype = "Miniscule-class (excluding Persistent stuck) chain"), size = 1,  alpha = 0.5) +
  geom_line(aes(y = either_happen, color = "Chain stuck (20+) or Miniscule-class, or both",
                linetype ="Chain stuck (20+) or Miniscule-class, or both"), size = 1.5, alpha = 0.5) +
  #geom_point(aes(y = `Stuck-chain`, color = "Stuck-chain"), size = 2) +
  scale_color_manual(values = c("Chain stuck (20+ iterations)" = "blue", 
                                "Persistent stuck chain" = "blue",
                                "Miniscule-class (excluding Persistent stuck) chain" = "red", 
                                "Chain stuck (20+) or Miniscule-class, or both" = "black")) +
  
  scale_linetype_manual(values = c( "Chain stuck (20+ iterations)" ="dashed",
                                    "Persistent stuck chain" = "solid", 
                                    "Miniscule-class (excluding Persistent stuck) chain" = "dotted", 
                                    "Chain stuck (20+) or Miniscule-class, or both" = "solid"),
                        guide = guide_legend(title = "Problematic behavior")) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",       # Adjust legend position
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),      # Remove panel border
        axis.line = element_line(linewidth = 0.5),  # Adjust axis line size
        strip.text = element_blank()) +
  #facet_grid(Smaller_step_sizes ~ Class_probabilities)
  facet_grid(Class_probabilities ~ Smaller_step_sizes)

plot1 + plot2 + plot_layout(ncol = 1)


ggplot(Sim_results) + 
  aes(x = Standard_deviations, group = Class_probabilities, color = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = "", x = " ") + 
  geom_line(aes(y = `Rhat>1.05`), size = 1.5, alpha = 0.5) +
  geom_line(aes(y = `Rhat>1.10`), size = 1, linetype = "dashed", alpha = 0.5) +
  #geom_point(aes(y = `Stuck-chain`), size = 2, shape = 16) +  # Use shape 16 (a filled circle) for Stuck-chain
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5)) +
  facet_grid(~Smaller_step_sizes)  # Separate panel for Stuck-chain


ggplot(Sim_results) + 
  aes(x = Standard_deviations, group = Class_probabilities, color = Class_probabilities) +
  scale_x_discrete(labels = abbreviate) + 
  labs(y = "Number of problematic chains", x = "Standard Deviations") + 
  geom_line(aes(y = `Stuck-chain`), size = 1.5, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5)) +
  facet_grid(~Smaller_step_sizes)  



