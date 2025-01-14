library(dplyr)
library(ggplot2)
library(tidyverse)


data <- data.frame(
  Gene = c("Pou5f1","Sox2","Klf4","Nanog","Myc"),
  rep1_0h = c(0, 112.022, 50.0957, 0, 0.156294),
  rep2_0h = c(0,134.472,106.338, 0,0),
  rep3_0h = c(0,200.848,63.0133,0,0),
  rep4_0h = c(0,205.872,82.2461,0,0.591259),
  rep5_0h = c(0,85.2485,3.25782,0,1.94681),
  rep6_0h = c(0,96.5367,3.35867,0,2.53887),
  rep1_4h = c(0,157.747, 9.49828, 0,17.33),
  rep2_4h = c(0,134.676,7.09814,0,19.661),
  rep1_16h = c(0,242.342, 64.6267,0.010199, 29.177), 
  rep2_16h =c(0,169.761,41.9704,0,30.3398), 
  rep1_24h = c(0,186.44,12.2367,0,17.9721), 
  rep2_24h =c(0,214.789,9.14078,0,18.3353), 
  rep1_36h =c(0,147.602,30.0493,0,9.25767), 
  rep2_36h =c(0,76.7887,15.6479,0,8.07092)
)

data$Gene <- factor(data$Gene, levels = c("Pou5f1", "Klf4", "Sox2", "Myc", "Nanog"))
data_long <- data %>%
  pivot_longer(cols = starts_with("rep"), 
               names_to = c("Replicate", "Time"), 
               names_pattern = "rep([0-9]+)_(\\d+h)", 
               values_to = "Expression")

data_long$Time <- factor(data_long$Time, levels = c("0h", "4h", "16h", "24h", "36h"))  # Adjust as needed
head(data_long)
data_summary <- data_long %>%
  group_by(Gene, Time) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE = sd(Expression) / sqrt(n()),  # Standard error
    Num_Replicates = n(),            # Number of replicates
    .groups = "drop"
  )

custom_colors <- c("Pou5f1" = "#3c2268", "Klf4" = "#ffae42",
                   "Sox2" = "#1d697c", "Myc" = "#d62728",
                   "Nanog" = "#9467bd")

head(data_summary) 

figure_name <- paste("mouse", "light.pdf", sep="_")
pdf(file =figure_name)
ggplot(data_summary, aes(x = Time, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # Stacked bar plot
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                position = position_dodge(width = 0.8), width = 0.2, color="grey") +  # Error bars
  geom_point(data = data_long, aes(x = Time, y = Expression, group = Gene),
             position = position_dodge(width = 0.8), size = 1, shape = 21, fill = "black") +  # Black points
  scale_fill_manual(values = custom_colors) +  # Custom bar colors
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300)) +
  labs( y = "Expression Level",
  fill = "Gene") +  # Custom fill legend for genes
  theme(axis.title.x = element_text(size = 20, face ="bold"),axis.title.y = element_text(size = 20, face ="bold"), axis.text.x = element_text(size = 20, face = "bold",angle = 45, hjust = 1),axis.text.y = element_text(size = 16,face = "bold"),
         axis.line.x = element_line(color = "black", size = 1), axis.line.y = element_line(color = "black", size = 1), legend.title = element_text(size = 16),  # Increase legend title font size
        legend.text = element_text(size = 20),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank())
dev.off()


####MOUSE NMDA 
data <- data.frame(
  Gene = c("Pou5f1","Sox2","Klf4","Nanog","Myc"),
  rep1_3h = c(0,225.527,56.1544,0,68.1088),
  rep2_3h = c(0,219.965,57.2214,0,78.5589),
  rep1_6h = c(0,129.759,130.246,0,52.0037),
  rep2_6h = c(0,141.794,196.757,0,77.1658),
  rep1_12h = c(0,136.028,27.8501,0,25.2171),
  rep2_12h = c(0,132.83,46.4528,0,30.1787),
  rep1_36h = c(0,121.545,5.854,0,8.67286), 
  rep2_36h = c(0,149.531,7.90845,0,10.6282),
  rep1_48h = c(0,97.1806,46.9563,0,79.7746),
  rep2_48h =c(0,97.9513,45.9081,0,10.923)
)

data$Gene <- factor(data$Gene, levels = c("Pou5f1", "Klf4", "Sox2", "Myc", "Nanog"))
data_long <- data %>%
  pivot_longer(cols = starts_with("rep"),
               names_to = c("Replicate", "Time"),
               names_pattern = "rep([0-9]+)_(\\d+h)",
               values_to = "Expression")
data_long$Time <- factor(data_long$Time, levels = c("3h", "6h", "12h", "36h", "48h"))  # Adjust as needed
head(data_long) 

data_summary <- data_long %>%
  group_by(Gene, Time) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE = sd(Expression) / sqrt(n()),  # Standard error
    Num_Replicates = n(),            # Number of replicates
    .groups = "drop"
  )


head(data_summary)
figure_name <- paste("mouse", "NMDA.pdf", sep="_")
pdf(file =figure_name)
ggplot(data_summary, aes(x = Time, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # Stacked bar plot
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                position = position_dodge(width = 0.8), width = 0.2, color="grey") +  # Error bars
  geom_point(data = data_long, aes(x = Time, y = Expression, group = Gene),
              position = position_dodge(width = 0.8), size = 1, shape = 21, fill = "black") +  # Black points
  scale_fill_manual(values = custom_colors) +  # Custom bar colors
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300)) +
  labs( y = "Expression Level",
  fill = "Gene") +  # Custom fill legend for genes
  theme(axis.title.x = element_text(size = 20, face ="bold"),axis.title.y = element_text(size = 20, face ="bold"), axis.text.x = element_text(size = 20, face = "bold",angle = 45, hjust = 1),axis.text.y = element_text(size = 16,face = "bold"),
         axis.line.x = element_line(color = "black", size = 1), axis.line.y = element_line(color = "black", size = 1), legend.title = element_text(size = 16),  # Increase legend title font size
        legend.text = element_text(size = 20),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank())

dev.off()




####Zebra Light 

custom_colors <- c("Pou5f1" = "#3c2268", "#1d697c", "Klf4" = "#ffae42",
                   "Sox2" = "#1d697c", "Myca" = "#d62728",
                   "Mycb" = "#899499", "Nanog" = "#9467bd")

data <- data.frame(
  Gene = c("Pou5f1","Sox2","Klf4","Nanog","Myca","Mycb"),
  rep1_0h = c(0, 77.5122, 2.13073,0, 18.348,81.4667),
  rep2_0h = c(1.06552,72.3649,3.04191,0,28.936,180.669),
  rep1_36h = c(0,68.0319,0.106719,0,12.3245,76.4843),
  rep2_36h = c(0.531904,102.122,0.0584342,0,22.4201,94.9254),
  rep3_0h = c(0,55.0572,0,0,8.35569,54.379),
  rep4_0h = c(0,53.4635,0,0,17.7717,97.9861),
  rep1_4h = c(0,120.271,0,0,17.8017,28.521),
  rep2_4h = c(0,42.8008,0,0,3.06707,77.9155),
  rep3_4h = c(0,63.2706,0.17922,0,13.2854,135.264),
  rep4_4h =c(0,55.7446, 0,0,3.02615,87.6132),
  rep1_10h = c(0,67.8266,0,0,0.110469,56.6311), 
  rep2_10h = c(0,57.2962,0,0,24.5277,40.8655), 
  rep3_10h = c(0,55.6711,0,0,27.3402,37.8433), 
  rep1_20h = c(0,80.8635,0,0,4.733592,187.907), 
  rep2_20h = c(0,86.986,0,0,2.42872,94.4524), 
  rep3_36h = c(0,67.5705,0,0,8.3521,114.934), 
  rep4_36h = c(0,72.9253,0,0,13.9685,153.914) )

data$Gene <- factor(data$Gene, levels = c("Pou5f1", "Klf4", "Sox2", "Myca", "Mycb", "Nanog"))

data_long <- data %>%
  pivot_longer(cols = starts_with("rep"),
               names_to = c("Replicate", "Time"),
               names_pattern = "rep([0-9]+)_(\\d+h)",
               values_to = "Expression")


data_long$Time <- factor(data_long$Time, levels = c("0h", "4h", "10h", "20h", "36h"))  # Adjust as needed
head(data_long)


table(data_long$Gene, data_long$Time)

data_summary <- data_long %>%
  group_by(Gene, Time) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE = sd(Expression) / sqrt(n()),  # Standard error
    Num_Replicates = n(),            # Number of replicates
    .groups = "drop"
  )

head(data_summary)
figure_name <- paste("zebra", "light.pdf", sep="_")
pdf(file =figure_name)
ggplot(data_summary, aes(x = Time, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # Stacked bar plot
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                position = position_dodge(width = 0.8), width = 0.2,color ="grey") +  # Error bars
  geom_point(data = data_long, aes(x = Time, y = Expression, group = Gene),
             position = position_dodge(width = 0.8), size = 1, shape = 21, fill = "black") +  # Black points
  scale_fill_manual(values = custom_colors) +  # Custom bar colors
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300)) +
  labs( y = "Expression Level",
  fill = "Gene") +  # Custom fill legend for genes
  theme(axis.title.x = element_text(size = 20, face ="bold"),axis.title.y = element_text(size = 20, face ="bold"), axis.text.x = element_text(size = 20, face = "bold",angle = 45, hjust = 1),axis.text.y = element_text(size = 16,face = "bold"),
         axis.line.x = element_line(color = "black", size = 1), axis.line.y = element_line(color = "black", size = 1), legend.title = element_text(size = 16),  # Increase legend title font size
        legend.text = element_text(size = 20),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank())
dev.off()


###Zebra NMDA 
data <- data.frame(
  Gene = c("Pou5f1","Sox2","Klf4","Nanog","Myca", "Mycb"),
  rep1_4h = c(0,90.0989,0,0,22.7341,115.529),
  rep2_4h = c(0,116.623,0,0,51.2476,260.316),
  rep1_10h = c(0,51.9135,0,0,16.013,143.885),
  rep2_10h = c(0,27.0225,0,0,13.7105,126.237),
  rep1_20h = c(0,86.0086,0.0284913,0,17.3453,118.4),
  rep2_20h = c(0,88.0794,0,0,17.0313,90.4823),
  rep3_20h = c(0,86.964,0,0,19.0615,89.3967),
  rep1_36h = c(0,45.7665,0,0,21.4406,140.483),
  rep2_36h = c(0,58.1092,0,0,57.6705,165.693),
  rep4_20h =c(0,30.6387,0.947452,0,26.489,141.845),
  rep5_20h = c(0.263547,66.4953,0.0343238,0,26.5617,190.982),
  rep3_36h = c(0.0120214, 46.7264,0.549743,0,36.8728,192.926),
  rep4_36h = c(0,62.2988,0.844783,0,34.8028,165.256) 
)

data$Gene <- factor(data$Gene, levels = c("Pou5f1", "Klf4", "Sox2", "Myca", "Mycb", "Nanog"))

data_long <- data %>%
  pivot_longer(cols = starts_with("rep"),
               names_to = c("Replicate", "Time"),
               names_pattern = "rep([0-9]+)_(\\d+h)",
               values_to = "Expression")

data_long$Time <- factor(data_long$Time, levels = c("4h", "10h", "20h", "36h"))  # Adjust as needed
head(data_long)


data_long_filtered <- data_long %>%
  filter(Expression != 0)


data_summary <- data_long %>%
  group_by(Gene, Time) %>%
  summarise(
    Mean_Expression = mean(Expression),
    SE = sd(Expression) / sqrt(n()),  # Standard error
    Num_Replicates = n(),            # Number of replicates
    .groups = "drop"
  )
head(data_summary)
figure_name <- paste("zebra", "NMDA.pdf", sep="_")
pdf(file =figure_name)
ggplot(data_summary, aes(x = Time, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +  # Stacked bar plot
  geom_errorbar(aes(ymin = Mean_Expression - SE, ymax = Mean_Expression + SE),
                position = position_dodge(width = 0.8), width = 0.2, color="grey") +  # Error bars
  geom_point(data = data_long, aes(x = Time, y = Expression, group = Gene),
             position = position_dodge(width = 0.2), size = 1, shape = 21, fill = "black") +  # Black points
  scale_fill_manual(values = custom_colors) +  # Custom bar colors
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300)) +
  labs( y = "Expression Level",
  fill = "Gene") +  # Custom fill legend for genes
  theme(axis.title.x = element_text(size = 20, face ="bold"),axis.title.y = element_text(size = 20, face ="bold"), axis.text.x = element_text(size = 20, face = "bold",angle = 45, hjust = 1),axis.text.y = element_text(size = 16,face = "bold"),
         axis.line.x = element_line(color = "black", size = 1), axis.line.y = element_line(color = "black", size = 1), legend.title = element_text(size = 16),  # Increase legend title font size
        legend.text = element_text(size = 20),  panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.border = element_blank())
dev.off()


