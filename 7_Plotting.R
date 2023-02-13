##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Plotting ouput of Switch tool per       ###
###            sample and per variant.                 ###
##########################################################

###########################
### Plotting per sample ###
###########################
#Plotting the Switch Error Rate for all phasing types you are comparing per sample

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



############################
### Plotting per variant ###
############################
#Plotting the Switch Error Rate of one phasing type along the chromosome


# Load libraries # 
library(SNPRelate) 
library(stats)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Load Data
df <- read.delim('pedigreephased.variant.switch.txt', header=F, sep='')
colnames(df) <- c('rm1', 'pos', 'rm2', 'rm3', 'trioSER')
df$trioSER <- ifelse(df$trioSER=='NaN', NA, df$trioSER)

# Remove rows with NAs for both readbase and trio
df_complete <- na.omit(df)
df_complete <- df_complete %>% select(pos, trioSER)
#for several columns:
#sel <- apply( df[,c("readbaseSER","trioSER")], 1, function(x) all(is.na(x)) )
#df_complete <- df[!sel,]

# Once we have the complete df without NAs, reorder it by ascending position of chr.
df_complete$pos <- as.numeric(df_complete$pos)
df_complete <-  df_complete %>% arrange(pos)

# We have 331,006 SNPs, so we can't plot per SNP and have to plot per window
# make a dataframe with windows in SNP indexes
win <- 1000
round.down <- nrow(df_complete) - (nrow(df_complete) %% win) # rounds down to previous 100th
change <- (nrow(df_complete) - (ncol(df_complete) %% win)) / win
# makes intervals of 100 SNPs 
int <- data.frame('start'=seq(1,round.down,win),
                  'end'=seq(1,round.down,win)+win-1)

# Define Var matrix 
var <- as.data.frame(matrix(NA,nrow=nrow(int), ncol=1))
pos <- c()

#Transform df_complete to a matrix
t <- df_complete %>% select(trioSER)
dat.t <- t(t)

# Sum 
for(i in 1:nrow(int)){
  # Sum number of mismatches
  var[i,] <- sum(dat.t[,seq(int$start[i],int$end[i])]) / win
  # Extract mean position of SNPs in the window
  pos[i] <- mean(df_complete$pos[int$start[i]:int$end[i]])
}

# Plot mismatches in snps windows:
#Function for transparency
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}


# Plot
jpeg(paste0('Single trio phasing vs Mendelian Inheritance along the Scaffold.jpg'), width=12,height=6,unit='in',quality = 1000,res=800)
par(mfrow=c(1,1))

plot(0,pch='',xlab="1000 snps window index (ss14)",ylab='switch error rate', xlim=c(0,max(df_complete$pos)),
     ylim=c(0,10), main = "Comparison of single trio pedigree phasing to Mendelian Inheritance along ss14")

lines(var[,1]~pos, col=add.alpha('black',.7)) # plot readbase

dev.off()
