# Process overview



This rmd should guide us through the process of pruning each population
to its core unrelated individuals to use in the reference panel. The
process happens on a subset of the overall beta matrix. Specifically we
want to get the overall beta along with a data frame that has
individuals and population information. Then we will subset this matrix,
make each value relative to the population mean and then prune using
inefficient R function developed by yours truly.

## Loading necessary files and functions

The package we need is SNPRelate to get the beta matrix.

``` r
library(SNPRelate)
```

    ## Loading required package: gdsfmt

    ## SNPRelate -- supported by Streaming SIMD Extensions 2 (SSE2)

The files we will load are:  
- The interval gds  
- A depth dataframe with mean DP / individual  
- A renaming list to transform the names in the DP df to the names in
the vcf

``` r
gds <- SNPRelate::snpgdsOpen('int4_allpops_miss10_mac5_ld.gds') #gds
dp <- read.table('BamCoverageOldNames.txt') # DP df
names(dp) <- c('INDV','MEAN_DEPTH') # col names
renamer <- read.table('newnames.list') # renaming list
names(renamer) <- c('old','new') # col names 
head(dp)
```

    ##      INDV MEAN_DEPTH
    ## 1 861530H    36.3520
    ## 2 862783H    32.8768
    ## 3 871031H    33.7349
    ## 4 871960H    33.6248
    ## 5 875626H    29.3849
    ## 6 877004H    27.3487

``` r
head(renamer)
```

    ##           old    new
    ## 1     861530H 861530
    ## 2     862783H 862783
    ## 3 P3_868867B1 868867
    ## 4     871031H 871031
    ## 5     871960H 871960
    ## 6     875626H 875626

The functions we will need come from [the pedigree
repository](https://github.com/topalw/pedigree). They are pretty slow
and shit but get the job done. If I ever get around to improving them it
would be nice.

``` r
find.names <- function(beta.mat,cutoff){
  df <- data.frame('id1'=NULL,
                   'id2'=NULL,
                   'beta'=NULL)
  end <- nrow(beta.mat) - 1
  for(i in 1:end){
    start <- i + 1
    for(j in start:nrow(beta.mat)){
      if(beta.mat[i,j] >= cutoff){
        tmp.df <- data.frame('id1'=rownames(beta.mat)[i],
                             'id2'=rownames(beta.mat)[j],
                             'beta'=beta.mat[i,j])
        df <- rbind(df,tmp.df)
      }
    }
  }
  return(df)
}
make.ind.df <- function(pairs){
  tmp <- data.frame('id'=unique(c(pairs$id1,pairs$id2)))  # make individual df
  # make negative to match decreasing order !!!
  tmp$dp <- -dp$MEAN_DEPTH[match(renamer$old[match(tmp$id,renamer$new)],dp$INDV)] # add dp
  tmp$links <- 0 #  add link count
  for(i in 1:nrow(tmp)){
    tmp$links[i] <- sum(c(pairs$id1,pairs$id2)==tmp$id[i])
  }
  tmp <- tmp[order(tmp[,3],tmp[,2], decreasing = T),]
  return(tmp)
}

# total function that prunes a beta mat on a cutoff of relatedness 
find.unrelated <- function(b,cutoff){
  bu <- b
  # initialize pairs   
  pairs <- find.names(b,cutoff) # get pairs 
  while(nrow(pairs)!=0){
    # find and remove worst candidate
    tmp <- make.ind.df(pairs)
    rm.ind <- tmp[1,1]
    bu <- bu[-match(rm.ind,rownames(bu)),-match(rm.ind,rownames(bu))]
    # recalculate pairs
    pairs <- find.names(bu,cutoff)
  }
  return(bu)
}
```

## Data preparation

For this process we need a beta object from SNPRelate and the pop
identifier df.  
We also *remove M026452 from everything* !

``` r
b <- SNPRelate::snpgdsIndivBeta(gds, autosome.only = F)
```

    ## Individual Inbreeding and Relatedness (beta estimator):
    ## Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
    ##     # of samples: 504
    ##     # of SNPs: 91,874
    ##     using 1 thread
    ## Individual Beta:    the sum of all selected genotypes (0,1,2) = 79250138
    ## CPU capabilities: Double-Precision SSE2
    ## Tue Nov 15 11:19:47 2022    (internal increment: 64128)
    ## [..................................................]  0%, ETC: ---        [==================================================] 100%, completed, 4s
    ## Tue Nov 15 11:19:51 2022    Done.

``` r
pop.df <- data.frame(
  'id' = b$sample.id,
  'pop' = ifelse(nchar(b$sample.id)>=5,'CH',substr(b$sample.id,0,2))
)
head(pop.df)
```

    ##       id pop
    ## 1 861530  CH
    ## 2 862783  CH
    ## 3 871031  CH
    ## 4 871960  CH
    ## 5 875626  CH
    ## 6 877004  CH

``` r
bs <- 'M026452' 
pop.df <- pop.df[pop.df$id != bs,]
dp <- dp[dp$INDV != bs,]
renamer <- renamer[renamer$new != bs,]

rels <- 0.25 / 2^c(0,1,2,3,4,5,6) # kinship levels 
```

## Pruning

Now we can subset this matrix to 1 matrix per population (except Georgia
which will be manually pruned since the average beta will be biased by
the fact that its exactly 1 family).

``` r
pops <- unique(pop.df$pop)
pops <- pops[pops != 'GE'] # this one manually
unrelated <- c()
rownames(b$beta) <- b$sample.id
rel <- data.frame('pop'=character(),'b'=double())
rel2 <- data.frame('pop'=character(),'b'=double())
rel3 <- data.frame('pop'=character(),'b'=double())
for(pop in pops){ # SLOW
  samples <- pop.df$id[pop.df$pop == pop] # get sample names 
  rel.b <- snpgdsIndivBeta(gds,autosome.only = F,sample.id =samples)$beta
  rownames(rel.b) <- samples # function needs rownames on matrix
  write.table(rel.b, paste0('beta_tables/',pop,'_beta.table'))
  tmp.unrel <- find.unrelated(rel.b, rels[4])
  unrelated <- c(unrelated,rownames(tmp.unrel)) # write unrelated names 
  tmp.b.rel <- b$beta[b$sample.id %in% samples,
                      b$sample.id %in% samples]
  n <- nrow(rel.b)
  rel <- rbind(rel, data.frame('pop' = rep(paste0(pop,'_pop'), (n^2-n)/2),
                                 'b' = rel.b[upper.tri(rel.b)])
                )
  n <- nrow(tmp.unrel)
  rel2 <- rbind(rel2, data.frame('pop' = rep(paste0(pop,'_unrel'), (n^2-n)/2),
                                 'b' =  tmp.unrel[upper.tri( tmp.unrel)])
                )
    n <- nrow(tmp.b.rel)
rel3 <- rbind(rel3, data.frame('pop' = rep(paste0(pop,'_all'), (n^2-n)/2),
                                 'b' =  tmp.b.rel[upper.tri( tmp.b.rel)])
                )
  }
```

A masochist that went through the code would see that we created a few
objects that appear interesting. One is a dataframe of a population
identifier and the beta estimates when using the population specific
mean. The other is a similar df but using a b.rel matrix which we
created (in the previous code chunk) by making beta relative to the
overall mean value. Finally we also calculated the beta before and after
individual pruning. This leads to a nice plot.

``` r
rel3 <- rbind(rel,rel2,rel3)
pdf('plots/3df_comparing_betas.pdf',height=8,width=16)
boxplot(rel3$b~factor(rel3$pop,levels=unique(rel3$pop)),las=2,xlab='',ylab='beta')
abline(h=0)
tmp <- aggregate(rel3$b , by=list(rel3$pop), FUN=mean)
tmp$y <- match(tmp$Group.1, levels(factor(rel3$pop,levels=unique(rel3$pop))))
points(tmp$x~tmp$y,col='orange',cex=.8,pch=16)
dev.off()
```

    ## png 
    ##   2

``` r
write.table(rel3, 'plots/3df.data',quote=F)
```

Finally we calculate the unrelated in Georgia manually using the overall
beta matrix.

``` r
geo <- pop.df[pop.df$pop =='GE',1]
gb <- b$beta[b$sample.id %in% geo,b$sample.id %in% geo]
diag(gb) <- NA
rownames(gb) <- geo
write.table(gb,'beta_tables/GE_beta.table')
colnames(gb) <- geo
corrplot::corrplot.mixed(gb, lower = 'number',upper='shade',
                         number.digits=4,)
```

![](finding_unrelated_files/figure-gfm/Georgians-1.png)<!-- -->

This leaves GEO1 and GEO06 to be added since they are the only ones that
obey the \< 0.03125 rule.

``` r
unrelated <- c(unrelated,c('GE01','GE06'))
# write unrelated all 
write.table(unrelated,paste0('unrel_lists/unrelated_all_',rels[4],'.list'),quote=F,row.names = F, col.names = F)
# write unrelated lib names 
write.table(renamer$old[match(unrelated,renamer$new)],paste0('unrel_lists/unrelated_all_',rels[4],'_libnames.list'),quote=F,row.names = F, col.names = F)
# write 1 file / pop in unrel lists 
per.pop <- data.frame('id'=unrelated,
                      'pop'=pop.df$pop[match(unrelated, pop.df$id)])
for(pop in unique(per.pop$pop)){
  write.table(per.pop$id[per.pop$pop==pop],
              paste0('unrel_lists/',pop,'_unrelated_',rels[4],'.list'),
              quote=F,row.names = F,col.names = F)
}
```
