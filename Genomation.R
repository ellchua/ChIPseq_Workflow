######################
## Set working directory
######################
getwd()

#and if not desired location, change with setwd()

setwd("/Users/ellora/Desktop/chip_analysis")

######################
## Load libraries
######################
library(GenomicRanges)
library(genomation)
library(RColorBrewer)

######################
## Create the file of peaks for windows
######################
## readGeneric has option to remove chromosomes with unusual names
## resize is a function in GenomicRanges so make sure to run the library for that first

combined_peaks <- readGeneric("/Users/ellora/Desktop/nfkb/IDR/idrValues.txt", chr = 1, start = 2, end = 3)
combined_peaks <- resize(combined_peaks, width = 1000, fix = "center")

RelA_peaks <- readNarrowPeak("/Users/ellora/Desktop/nfkb/RelA_rep1/macs2/RelA_rep1_sorted_peaks.narrowPeak")
RelA_peaks <- resize(RelA_peaks, width = 1000, fix = "center")

RelA_RelB_peaks <- union(combined_peaks, RelA_peaks)
######################
## Create score matrix list
######################

sml <- ScoreMatrixList(targets = c("~/Desktop/nfkb/RelB_rep1/bwa/RelB_rep1_sorted.bam", "~/Desktop/nfkb/RelB_rep2/bwa/RelB_rep2_sorted.bam"), 
                       windows = combined_peaks, bin.num = 50)
sml_scaled <- scaleScoreMatrixList(sml)

remunus_sml <- ScoreMatrixList(targets = c("~/Desktop/nfkb/RelB_rep1/bwa/RelB_rep1_sorted.bam", "~/Desktop/nfkb/RelB_rep2/bwa/RelB_rep2_sorted.bam"), 
                               windows = remunus_peaks, bin.num = 50)
remunus_sml_scaled <- scaleScoreMatrixList(remunus_sml)

RelA_RelB_sml <- ScoreMatrixList(targets = c("~/Desktop/nfkb/RelB_rep1/bwa/RelB_rep1_sorted.bam", "~/Desktop/nfkb/RelB_rep2/bwa/RelB_rep2_sorted.bam",
                                             "~/Desktop/nfkb/RelA_rep1/bwa/aln_RelA_rep1_sorted.bam"), 
                                 windows = RelA_RelB_peaks, bin.num = 50)
RelA_RelB_sml_scaled <- scaleScoreMatrixList(RelA_RelB_sml)

RelAB_sml <- ScoreMatrixList(targets = c("~/Desktop/nfkb/RelB_rep1/bwa/RelB_rep1_sorted.bam", "~/Desktop/nfkb/RelB_rep2/bwa/RelB_rep2_sorted.bam",
                                             "~/Desktop/nfkb/RelA_rep1/bwa/aln_RelA_rep1_sorted.bam"), 
                                 windows = remunus_RelAB_peaks, bin.num = 50)
RelAB_sml_scaled <- scaleScoreMatrixList(RelAB_sml)

######################
## Plot heat map
######################
## use display.brewer.all() to display all brewer palettes

multiHeatMatrix(sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Blues"))
multiHeatMatrix(remunus_sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Blues"))
multiHeatMatrix(RelA_RelB_sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Blues"))
multiHeatMatrix(RelAB_sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Blues"))
multiHeatMatrix(RelAB_sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Purples"))      

######################
## Check for kmeans center
######################
foo = as(RelAB_sml_scaled[[1]], "matrix")
foo = cbind(as(RelAB_sml_scaled[[1]], "matrix"), as(RelAB_sml_scaled[[2]], "matrix"), 
            as(RelAB_sml_scaled[[3]], "matrix"))
bar = sapply(1:10, function(x){kmeans(foo, centers=x)$tot.withinss})
plot(bar, type="l")

######################
## Plot clustered heat map
######################
ham = multiHeatMatrix(RelAB_sml_scaled, xcoords=c(-500, 500), clustfun=function(x) kmeans(x, centers=5)$cluster, winsorize=c(0,95),  
                      common.scale=TRUE, col=brewer.pal(n = 9, name = "Reds"))
head(ham)
table(ham)
ham$cluster                         