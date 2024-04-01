
#                              ANALYSE MICROBOME DADA2 & PHYLOSEQ

# First install a recent version of Rstudio: https://www.rstudio.com/products/rstudio/download/#download

# Install last version of R : https://cran.r-project.org/bin/windows/base/
#Select it in Tools Global Options Rversion in Rstudio.

#Installation of packages for Microbiome analyses usig DADA2

#command to install using BioManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("Biostrings")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("dada2")
BiocManager::install("DESeq2")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("devtools")
BiocManager::install("Tax4Fun")
BiocManager::install("Rtools")
install.packages("microbial")
BiocManager::install("edgeR")

#for Tax4Fun, install manually (download archive from http://tax4fun.gobics.de/)
#Tax4Fun_0.3.1.tar.gz
#install qiime2R
devtools::install_github("jbisanz/qiime2R")

#Barplots & Venn diagram packages
install.packages("microbiome")

install.packages("git2r")
install.packages("remotes")
remotes::install_github("Russel88/MicEco")




#Load libraries

library(phyloseq)
library(DECIPHER)
library(dada2)
library(vegan)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(msa)
library(phangorn)
library(microbial)
#install Tax4Fun manually after downloading from http://tax4fun.gobics.de/
#(need to install first Rtools from https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html)
library(Tax4Fun)
library(qiime2R)
library("MicEco")

#starting the dada2 pipeline
setwd("C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_18S/18SALL")

path <- "C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_18S/18SALL/HN00171775"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)



# Forward and reverse fastq filenames have format: SAMPLENAME-18S_1.fastq and SAMPLENAME-18S_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

#Filter and trim: Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter: Macrogne primers used: 18S V4 region V4F (20 pb) - CCAGCAGCCGCGGTAATTCC; V4R (18 pb)- ACTTTCGTTCTTGATTAA

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(20, 18), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
# On Windows set multithread=FALSE; attention ici on trime les primers de 20 pb en amont de chaque sequence for et Rev (trimleft)


head(out) 
#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

#Remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted length. 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 349:433]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## [1]  20 232
sum(seqtab.nochim)/sum(seqtab2)
## [1] 0.964263

write.csv(seqtab.nochim, "seqtabnochimAll18S.csv")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Stat_assemblyDADA2_Kakada18SAll.csv")



#Assign taxonomy (download first silva_nr_v138_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v138.fa.gz from silva website)
taxa <- addSpecies(taxa, "C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/silva_species_assignment_v138.fa.gz")


#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(taxa, file="taxaAll18S_species.RDS")
saveRDS(seqtab.nochim, file="seqtabnochimAll18S.RDS")

#load metadata file
samdf <- read.csv2("C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_18S/metadata_all_corrected.csv")

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE to verify correspondance 

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Field", "Organ", "Disease", "OrgHealth", "FieldOrgHealth") # indicate column names from your metadata that you want to keep
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

#save files in RDS format (easier to reload later in R with loadRDS function)

saveRDS(seqtab.nochim, "seqtab.nochim.RDS")
saveRDS(taxa, "taxa.RDS")
write.csv(taxa, "taxa18S.csv")
write.csv(seqtab.nochim, "seqtab.nochim18S.csv")


# extract info in a single table (without metadata)
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa6<-as.data.frame(taxa)
taxa6$ID <- rownames(taxa)
all_data5  <- merge(taxa6,seqtab6,by='ID')
write.csv(all_data5, "Oeum2024_ALL18S.csv")


#phyloseq object without phylogenetic tree

ps_18S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_18S


#Add ASV numbering instead of fasta sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_18S))
names(dna) <- taxa_names(ps_18S)
ps_18S  <- merge_phyloseq(ps_18S, dna)
taxa_names(ps_18S) <- paste0("ASV", seq(ntaxa(ps_18S)))
ps_18S



#remove mitochondria and chloroplast and Phylum Phragmoplastophyta (plants and algae)
ps <- ps_18S %>% subset_taxa(Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" |is.na(Order)) 
ps <- ps %>% subset_taxa(Phylum!="Phragmoplastophyta"  |is.na(Phylum))

#check if removed correctrly by using plotbar
plotbar(ps, level = "Phylum", color = NULL,group = "FieldOrgHealth",
        top = 30,return = FALSE,fontsize.x = 5,  fontsize.y = 12)

#Remove sample from dataset (here control)
ps <- prune_samples(sample_names(ps) != "RVCTRL", ps) 


#remove a list of samples (removing leaf samples since they contain onyl rice reads)
#First, create a list of the samples that you want to remove
Samples_toRemove <- c("PSL1","PSL2","PSL3","PSL4","PSL5","PSL6","PSL7","PSL8","PSL9","PSL10","PSL11", "PSL12", "PSL13", "PSL14", "PSL15", "PSL16","PSL17","PSL18","PSL19", "PSL20", "RVL1", "RVL2", "RVL3","RVL4", "RVL5", "RVL6","RVL7", "RVL8", "RVL9","RVL10", "RVL11", "RVL12", "RVL13", "RVL14","RVL15", "RVL16", "RVL17","RVL18", "RVL19", "RVL20")

#To see what samples get removed, run the following; note, I have a column called sample_id

subset_samples(ps, sample_id %in% Samples_toRemove)
#This will return a ps object that contains the samples you want to remove

#To remove those from your phyloseq object
ps_noleaf <- subset_samples(ps, !(sample_id %in% Samples_toRemove))
#This will return a ps object with the samples removed
ps_noleaf

#plotbars at different taxonomic ranks and top

plotbar(ps_noleaf, level = "Phylum", color = NULL,group = "FieldOrgHealth",
        top = 30,return = FALSE,fontsize.x = 5,  fontsize.y = 12)

plotbar(ps_noleaf, level = "Genus", color = NULL,group = "FieldOrgHealth",
        top = 30,return = FALSE,fontsize.x = 5,  fontsize.y = 12)


# order
desired_order <- list("PSRo1","PSRo2","PSRo3","PSRo4","PSRo5","PSRo6","PSRo7","PSRo8","PSRo9","PSRo10","PSRo11","PSRo12","PSRo13","PSRo14","PSRo15","PSRo16","PSRo17","PSRo18","PSRo19","PSRo20","RVRo1","RVRo2","RVRo3","RVRo4","RVRo5","RVRo6","RVRo7","RVRo8","RVRo9","RVRo10","RVRo11","RVRo12","RVRo13","RVRo14","RVRo15","RVRo16","RVRo17","RVRo18","RVRo19","RVRo20","PSRh1","PSRh2","PSRh3","PSRh4","PSRh5","PSRh6","PSRh7","PSRh8","PSRh9","PSRh10","PSRh11","PSRh12","PSRh13","PSRh14","PSRh15","PSRh16","PSRh17","PSRh18","PSRh19","PSRh20","RVRh1","RVRh2","RVRh3","RVRh4","RVRh5","RVRh6","RVRh7","RVRh8","RVRh9","RVRh10","RVRh11","RVRh12","RVRh13","RVRh14","RVRh15","RVRh16","RVRh17","RVRh18","RVRh19","RVRh20")


p <- plot_bar(ps_noleaf, fill="Phylum") 
p$data$Sample <- factor(p$data$FieldOrgHealth, levels = desired_order)
print(p)



#Create files for each data field
ps_RV <- subset_samples(ps_16S, Field=="Rovieng")
ps_PS <- subset_samples(ps_16S, Field=="PreakSdei")

#Create files per organ per field
ps_RVL <- subset_samples(ps_RV, Organ=="Leaf")
ps_RVL

ps_RVRo <- subset_samples(ps_RV, Organ=="Root")
ps_RVRo

ps_RVRh <- subset_samples(ps_RV, Organ=="Rhizo")
ps_RVRh

ps_PSL <- subset_samples(ps_PS, Organ=="Leaf")
ps_PSL

ps_PSRo <- subset_samples(ps_PS, Organ=="Root")
ps_PSRo

ps_PSRh <- subset_samples(ps_PS, Organ=="Rhizo")
ps_PSRh

ps_LeavesALL <- subset_samples(ps_16S, Organ=="Leaf")
ps_LeavesALL

ps_RootALL <- subset_samples(ps_16S, Organ=="Root")
ps_RootALL

ps_RhizoALL <- subset_samples(ps_16S, Organ=="Rhizo")
ps_RhizoALL


#plot bars
plot_bar(ps_noleaf, fill="Phylum")


# order
desired_order <- list("PSL1","PSL2","PSL3","PSL4","PSL5","PSL6","PSL7","PSL8","PSL9","PSL10","PSL11","PSL12","PSL13","PSL14","PSL15","PSL16","PSL17","PSL18","PSL19","PSL20","RVL1","RVL2","RVL3","RVL4","RVL5","RVL6","RVL7","RVL8","RVL9","RVL10","RVL11","RVL12","RVL13","RVL14","RVL15","RVL16","RVL17","RVL18","RVL19","RVL20","PSRo1","PSRo2","PSRo3","PSRo4","PSRo5","PSRo6","PSRo7","PSRo8","PSRo9","PSRo10","PSRo11","PSRo12","PSRo13","PSRo14","PSRo15","PSRo16","PSRo17","PSRo18","PSRo19","PSRo20","RVRo1","RVRo2","RVRo3","RVRo4","RVRo5","RVRo6","RVRo7","RVRo8","RVRo9","RVRo10","RVRo11","RVRo12","RVRo13","RVRo14","RVRo15","RVRo16","RVRo17","RVRo18","RVRo19","RVRo20","PSRh1","PSRh2","PSRh3","PSRh4","PSRh5","PSRh6","PSRh7","PSRh8","PSRh9","PSRh10","PSRh11","PSRh12","PSRh13","PSRh14","PSRh15","PSRh16","PSRh17","PSRh18","PSRh19","PSRh20","RVRh1","RVRh2","RVRh3","RVRh4","RVRh5","RVRh6","RVRh7","RVRh8","RVRh9","RVRh10","RVRh11","RVRh12","RVRh13","RVRh14","RVRh15","RVRh16","RVRh17","RVRh18","RVRh19","RVRh20")


p <- plot_bar(ps, fill="Phylum") 
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)

desired_order2 <- list("PreakSdeiRootNormal","PreakSdeiRootDisease","RoviengRootNormal","RoviengRootDisease","PreakSdeiRhizoNormal","PreakSdeiRhizoDisease","RoviengRhizoNormal","RoviengRhizoDisease")

p1 <- plotbar(ps_noleaf, level = "Order", color = NULL,group = "FieldOrgHealth",
              top = 20,return = FALSE,fontsize.x = 5,  fontsize.y = 12)
p1$data$FieldOrgHealth <- factor(p1$data$FieldOrgHealth, levels = desired_order2)
print(p1)


p1 <- plotbar(ps_18S_1000F, level = "Genus", color = NULL,group = "FieldOrgHealth",
              top = 20,return = FALSE,fontsize.x = 5,  fontsize.y = 12)
p1$data$FieldOrgHealth <- factor(p1$data$FieldOrgHealth, levels = desired_order2)
print(p1)



#Rarefy at 1000
set.seed(1234)
ps_18S_1000F <- rarefy_even_depth(ps_noleaf, sample.size = 1000,
                                   rngseed = F, replace = TRUE) 
ps_18S_1000F    


#Reload 18S data from RDS files to create figures                          
seqtab.nochim <- readRDS("seqtabnochimAll18S.RDS")
taxa <- readRDS("taxaAll18S_species.RDS")
samdf <- read.csv2("C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_18S/18SALL/metadata_all_corrected.csv")

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE pour verifier que les noms concordent

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Field", "Organ", "Disease", "FieldHealth", "OrgHealth", "FieldOrgHealth","DiseasePercentage") # a changer avec le nom de tes colonnes de m?tadonn?es
samdf <- samdf[rownames(seqtab.nochim), keep.cols]
#phyloseq object without phylogenetic tree

ps_18S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_18S
ps <- ps_18S

#Add ASV numbering instead of fasta sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps  <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#Rarefaction curves

set.seed(42)

calculate_rarefaction_curves <- function(ps, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(ps, measures, depth) {
    if(max(sample_sums(ps)) < depth) return()
    ps <- prune_samples(sample_sums(ps) >= depth, ps)
    
    rarified_ps <- rarefy_even_depth(ps, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_ps, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, ps = ps, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}


rarefaction_curve_data <- calculate_rarefaction_curves(ps, c('Observed', 'Shannon'), rep(c(1, 10, 50, 100, 200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000), each = 10))

summary(rarefaction_curve_data)


rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(ps)), by.x = 'Sample', by.y = 'row.names')


ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = FieldOrgHealth,
    group = Sample
  )
) + geom_line(
) + geom_pointrange(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + theme_bw()

#Removing grey font and dots
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Organ,
    group = Sample
  )
) + geom_line(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + theme_bw()

#save curves to an object "p1"  for figure compilation
p1 <- ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Organ,
    group = Sample
  )
) + geom_line(
) + facet_wrap(
  facets = ~ Measure,
  scales = 'free_y'
) + theme_bw()



