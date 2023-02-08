##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: Plot SER along chromosome.              ###
##########################################################
# Load libraries # 
library(SNPRelate) 
library(stats)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Load Data
df <- read.delim('TRIO_Super-Scaffold_14.variant.switch.txt', header=F, sep='')
colnames(df) <- c('rm1', 'pos', 'rm2', 'rm3', 'trioSER')
df$trioSER <- ifelse(trio$df=='NaN', NA, df$trioSER)

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
