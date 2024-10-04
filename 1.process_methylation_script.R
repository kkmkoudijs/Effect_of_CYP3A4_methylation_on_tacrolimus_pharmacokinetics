# From QC report: array type = 'MethylEPIC v1 manifest B5'
# Initial number of samples = 23
# Initial number of CpG probes = 865,918

library(methylationArrayAnalysis)
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# use hg38, b5 version: https://github.com/achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38
# Downloaded from:
# library(devtools)
# install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
# install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38") 
# instead of 'IlluminaHumanMethylation450kanno.ilmn12.hg19' and 'IlluminaHumanMethylation450kmanifest' in tutorial
# use these packages: 


library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)

# The manifest contains all of the annotation information for each of the CpG probes on the ann865k array.
ilm10b5.hg38 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
head(ilm10b5.hg38)

# read in the sample sheet for the experiment
targets <- read.metharray.sheet(base = "Input", pattern="MET2020-244-014.csv")
targets # Basename column is not filled
targets$Basename <- paste0(targets$Slide,"_",targets$Array)

# read in the raw data from the IDAT files
rgSet <- read.metharray.exp(base = "Input/idats", targets)
sampleNames(rgSet) <- targets$Sample_Name # give the samples descriptive names
rgSet # Detects as 'IlluminaHumanMethylationEPIC' (correct afaik) but with 'ilm10b4.hg19' annotation 

rgSet@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38") # Force this annotation
annotation(rgSet)

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

table(targets$Array) # 8 different slides
table(targets$Slide) # 3 different arrays

# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Slide)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Slide)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Slide)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Slide)), fill=pal, 
       bg="white")

qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Slide, pdf="Output/qcReport.pdf")
warnings() # Type conversion errors, does not seem serious

# remove poor quality samples
keep <- colMeans(detP) < 0.05 
rgSet <- rgSet[,keep]
rgSet # All 23 samples are OK

# Note: preprocessQuantile function (Touleimat and Tost 2012) is 
# more suited for datasets where you do not expect global differences between your samples, for example a single tissue.
# Note that after normalisation, the data is housed in a GenomicRatioSet object.

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# Figure S1: visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Slide,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Slide)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Slide,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Slide)), 
       text.col=brewer.pal(8,"Dark2"))

# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Slide)])
legend("top", legend=levels(factor(targets$Slide)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Slide)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Slide)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Slide)), text.col=pal, 
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Slide)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Slide)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Slide)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Slide)), text.col=pal,
       cex=0.7, bg="white")


# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep) # Kept 829468 out of 865859 probes (95.8%); 

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ilm10b5.hg38$Name[ilm10b5.hg38$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

mSetSqFlt
featureNames(mSetSqFlt)

included_probe_indexes <- 
  c(
    which(grepl("CYP3A4",ilm10b5.hg38[,"UCSC_RefGene_Name"]))
  )


Included_probes <- 
  as.data.frame(ilm10b5.hg38[included_probe_indexes,colnames(ilm10b5.hg38)])


Included_probes <-
  Included_probes[setdiff(x = 1:nrow(Included_probes), y = which(grepl("CYP3A43",Included_probes$UCSC_RefGene_Name))),]


Included_probes <-
  Included_probes[,c("UCSC_RefGene_Name","Name","CHR_hg38","Start_hg38","End_hg38","Strand_hg38")]

Start_CYP3A4 <- 99756960
Included_probes$Relative_start <- as.integer(Included_probes$Start_hg38) - Start_CYP3A4
Included_probes$Relative_end <- as.integer(Included_probes$End_hg38) - Start_CYP3A4


Included_probes <- 
  Included_probes[order(Included_probes$Start_hg38),]



All_beta_values <- getBeta(mSetSq)
All_beta_values <- t(All_beta_values[rownames(All_beta_values) %in% Included_probes$Name,])

All_M_values <- getM(mSetSq)
All_M_values <- t(All_M_values[rownames(All_M_values) %in% Included_probes$Name,])

rownames(Included_probes) <- 1:nrow(Included_probes)
colnames(Included_probes)[2] <- "Probe_name"
Included_probes <- Included_probes[order(Included_probes$Probe_name),]
All_beta_values <- All_beta_values[,order(colnames(All_beta_values))]

write.table(Included_probes,  file = paste0("Output/Table_S2_Included_probes.csv"), row.names = F, sep = ";")
write.table(cor(All_beta_values),  file = paste0("Output/Table_S3_cor_All_beta_values.csv"), row.names = T, sep = ";")

saveAsCsv <- function(object, name){
  
  CSV_output <- as.data.frame(object)
  CSV_output$sample_name <- rownames(CSV_output)
  CSV_output <- CSV_output[,c(ncol(CSV_output),1:(ncol(CSV_output)-1))]
  write.table(CSV_output,  file = paste0("Output/",name,".csv"), row.names = F, sep = ";")
  
}
saveAsCsv(All_beta_values, "All_beta_values")
saveAsCsv(All_M_values, "All_M_values")