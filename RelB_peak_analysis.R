######################
## Set working directory
######################
getwd()

#and if not desired location, change with setwd()

setwd("/Users/ellora/Desktop/chip_analysis")

######################
## Load libraries 
######################
library(GenomicFeatures)
library(biomaRt)

######################
## Load gene model
######################
# returns a GRanges object of genes i.e. every entry is a gene and its start and end location
# but all this tells you is theres a gene there with a specific gene identifier
# so it's much more useful to actually know what the gene name is. Is it a protein-coding gene? Is it a lincRNA?

hg.gtf.db <- makeTxDbFromGFF("/Users/ellora/Desktop/nfkb/Homo_sapiens.GRCh38.94.chr.gtf", format="gtf")
ensembl.genes = genes(hg.gtf.db)
ensembl.transcripts = transcripts(hg.gtf.db)

######################
## Create mart 
######################
# connects to ensembl and says we're using human (hsapiens) and version 94 of the annotation

human = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",host="asia.ensembl.org", 
                   dataset="hsapiens_gene_ensembl", version="94")

######################
## Query biomaRt
######################
# says go to biomaRt and return me a data.frame containing lots of attributes 
# (things we're interested in about the gene) using a set of gene identifiers (ensembl.genes$gene_id)

bm.annotations = getBM(attributes=c("ensembl_gene_id", "entrezgene", "gene_biotype", "hgnc_symbol", "description"), 
                       mart=human, filters="ensembl_gene_id", values=ensembl.genes$gene_id, 
                       uniqueRows=TRUE)

######################
## Add metadata
######################
# these few lines simply match the genes you have with their metadata and add it to the "ensembl.genes" object

ensembl.genes$hgnc_symbol = bm.annotations$hgnc_symbol[match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]
ensembl.genes$gene_biotype = bm.annotations$gene_biotype[match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]
ensembl.genes$description = bm.annotations$description[match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]
ensembl.genes$entrezgene = bm.annotations$entrezgene[match(ensembl.genes$gene_id, bm.annotations$ensembl_gene_id)]

######################
## Filter conservative IDR peaks
######################
# read the output file obtained from idr command into a table

idrValues <- read.delim("/Users/ellora/Desktop/nfkb/IDR/idrValues.txt", sep="\t", header=FALSE)

# IDR cutoff

idr_threshold <- 0.05

# scaled IDR score

scaled_idr_val <- round(min(-125*log2(idr_threshold), 1000))

# We filter out the peaks that pass the idr cutoff of 0.05
# Here I use the score int but Nathan uses global IDR

passPeaks <- idrValues[idrValues$V5 >= scaled_idr_val,] #remember we want every column, hence the comma at the end

######################
## Create GRanges object
######################
# For easy addition of metadata

idr.gr <- GRanges(passPeaks[,1], IRanges(passPeaks[,2], passPeaks[,3]))

#remove weird chromosomes (keepStandardChromosomes work on GRanges object)
#another way is to subset peaks -- peaks.gr <- peaks.gr[seqnames(peaks.gr) %in% paste("chr", c(1:22, "X", "Y", "M")), sep="")]

idr.gr <- keepStandardChromosomes(idr.gr, species = "Homo sapiens", pruning.mode = "coarse")

######################
## Write out file for GREAT analysis
######################

great_idr <- idr.gr
great_idr <- as.data.frame(great_idr)
write.table(great_idr, file="/Users/ellora/Desktop/chip_analysis/great_idr.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

######################
## Initialize metadata columns
######################

idr.gr$gene_id = ""
idr.gr$hgnc_symbol = ""
idr.gr$gene_biotype = ""
idr.gr$description = ""
idr.gr$overlap_promoter = FALSE
idr.gr$promoter = ""
idr.gr$nearest_gene = ""
idr.gr$nearest_gene_hgnc_symbol = ""
idr.gr$nearest_gene_biotype = ""
idr.gr$nearest_gene_description = ""
idr.gr$nearest_gene_distance = ""

######################
## Find overlaps and add metadata
######################
# Finding overlaps between peaks in idr.gr and entries in ensembl.genes

overlaps <- findOverlaps(idr.gr, ensembl.genes)

idr.gr[queryHits(overlaps)]$gene_id = ensembl.genes[subjectHits(overlaps)]$gene_id
idr.gr[queryHits(overlaps)]$hgnc_symbol = ensembl.genes[subjectHits(overlaps)]$hgnc_symbol
idr.gr[queryHits(overlaps)]$gene_biotype = ensembl.genes[subjectHits(overlaps)]$gene_biotype
idr.gr[queryHits(overlaps)]$description = ensembl.genes[subjectHits(overlaps)]$description

# Finding overlaps between peaks in idr.gr and promoters in ensemble.genes
# ?promoters

prom_overlaps <- findOverlaps(idr.gr, promoters(ensembl.genes))

idr.gr[queryHits(prom_overlaps)]$overlap_promoter = TRUE # if peak overlaps with promoter region, TRUE
idr.gr[queryHits(prom_overlaps)]$promoter = ensembl.genes[subjectHits(prom_overlaps)]$gene_id

# If not overlapping with gene or promoting, finding nearest gene
# nearest(idr.gr, ensembl.genes) - for each peak what is the nearest gene
# nearest(ensembl.genes, idr.gr) - for each gene what is the nearest peak

nearest.gene <- nearest(idr.gr, ensembl.genes)
idr.gr$nearest_gene = ensembl.genes[nearest.gene]$gene_id
idr.gr$nearest_gene_hgnc_symbol = ensembl.genes[nearest.gene]$hgnc_symbol
idr.gr$nearest_gene_biotype = ensembl.genes[nearest.gene]$gene_biotype
idr.gr$nearest_gene_description = ensembl.genes[nearest.gene]$description

# Finding distance to nearest gene

dist <- distanceToNearest(idr.gr,ensembl.genes)

idr.gr[queryHits(dist)]$nearest_gene_distance = mcols(dist)$distance

######################
## Write file out
######################

peaks_info <- as.data.frame(idr.gr)
write.csv(peaks_info, file="/Users/ellora/Desktop/chip_analysis/peaks_info.csv")

######################
## Genomation
######################
library(genomation)

## readGeneric has option to remove chromosomes with unusual names
combined_peaks <- readGeneric("/Users/ellora/Desktop/nfkb/IDR/idrValues.txt", chr = 1, start = 2, end = 3)
combined_peaks <- resize(combined_peaks, width = 1000, fix = "center") #resize is a function in GenomicRanges so make sure to run the library for that first
sml <- ScoreMatrixList(targets = c("~/Desktop/nfkb/RelB_rep1/bwa/RelB_rep1_sorted.bam", "~/Desktop/nfkb/RelB_rep2/bwa/RelB_rep2_sorted.bam"), windows = combined_peaks, bin.num = 50)
sml_scaled <- scaleScoreMatrixList(sml)


install.packages("RColorBrewer")
library("RColorBrewer")
multiHeatMatrix(sml_scaled, xcoords = c(-500, 500), col = brewer.pal(n = 9, name = "Blues"))


