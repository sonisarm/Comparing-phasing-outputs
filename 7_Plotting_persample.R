##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Plotting ouput of Switch tool.          ###
##########################################################

# Load libraries
library(ggplot2)
library(tidyverse)

# Construct a default theme for plotting:
default.theme <- theme_bw() +
  theme(plot.title = element_text(size = 8, colour = "black", face = "bold", hjust = 0),
        axis.title = element_text(size = 8, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.line = element_line(size = 0.25, colour = "black"),
        legend.title = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 7, colour = "black"),
        legend.key.size = unit(0.3, "cm"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, colour = "black", face = "bold"))

# Set Directory
# CHANGE ACCORDING TO YOUR DIRECTORY #
setwd('/Users/soniasarmiento/Documents/Phasing')

# Load switch outputs, one per phasing type (e.g. read-base, pedigree (trio + multiple trio))
# CHANGE ACCORDING TO YOUR FILES #
phasing1 <- read.delim('outputexample.sample.switch.txt.gz', sep='', header=F) #e.g. load readbase
phasing2 <- read.delim('outputexample.sample.switch.txt.gz', sep='', header=F) #e.g. load single trio
phasing3 <- read.delim('outputexample.sample.switch.txt.gz', sep='', header=F) #e.g. load multiple trio

#Add descriptive columns to the database
  # Phasing 1 #
phasing1$type <- 'phasing1'
phasing1$scaffold <- 'scaffold 14'
colnames(phasing1)[1]<- 'indv'
colnames(phasing1)[4]<- 'SER'

  # Phasing 2 #
phasing2$type <- 'phasing2'
phasing2$scaffold <- 'scaffold 14'
colnames(phasing2)[1]<- 'indv'
colnames(phasing2)[4]<- 'SER'

  # Phasing 3 #
phasing3$type <- 'phasing3'
phasing3$scaffold <- 'scaffold 14'
colnames(phasing3)[1]<- 'indv'
colnames(phasing3)[4]<- 'SER'


# Bind the datasets
df <- rbind(phasing1, phasing2, phasing3)
df <- df %>% group_by(type)

# get the mean per phasing type
summary_df <- df %>% group_by(type) %>% summarize(m=mean(SER))


# PLOT SER per sample #
jpeg(paste0('Results_Sample_Switch',win,'.jpg'), width=12, height=6, unit='in',quality = 1000,res=800, antialias = "cleartype")
ggplot(df) +
  geom_line(aes(x=scaffold, y=SER, color=type), size = 0.5) + 
  xlab("Scaffold") + ylab("Switch error rate (%)") +
  default.theme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5)) +
  theme(legend.position="top") +
  scale_color_manual(values=c("#00AFBB", "#E7B800", '#FF4500')) +   # as many color as phasing types
  geom_text(data=summary_df,
            aes(x=scaffold, y=m, label=round(m,2), group=type, color=type),
            position = position_dodge(width = 1.3), 
            vjust = -6.5, 
            size = 3, 
            parse = TRUE)
dev.off()
