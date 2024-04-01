
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
install.packages("metagMisc")

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

setwd("C:/your working directory")

#starting the dada2 pipeline
path <- "C:/"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)



# Forward and reverse fastq filenames have format: SAMPLENAME-16S_1.fastq and SAMPLENAME-16S_2.fastq
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

#Filter:  Macrogen 16S primers are 341F :  CCTACGGGNGGCWGCAG (17 nt) and  805R GACTACHVGGGTATCTAATCC (21 nt)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(17, 21), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
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
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 371:435]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## [1]  20 232
sum(seqtab.nochim)/sum(seqtab2)
## [1] 0.964263

write.csv(seqtab.nochim, "seqtabnochim.csv")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "trackfile.csv")


#Assign taxonomy (download first silva_nr_v132_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v132.fa.gz from silva website)
taxa <- addSpecies(taxa, "C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/silva_species_assignment_v132.fa.gz")

#Reload seqtab.nochim to test new taxonomi assignment with RDP and silva update v138
seqtab.nochim <- readRDS(file ="C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Microbiota_Rovieng_PreakSdei2021/fastQ_16S/All16S/seqtabnochimAll16S.RDS")
taxa <- readRDS(file="taxaAll16S_species.RDS")

#Assign taxonomy (different databases options given, most recent and updated is silva138)
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v138.fa.gz from silva website)
taxa <- addSpecies(taxa2, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_species_assignment_v138.fa.gz")

#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save seqtab and taxa as RDS files (easier to reload later in R)
saveRDS(taxa2, file="taxaAll16S_species_silva138.RDS")
saveRDS(seqtab.nochim, file="seqtabnochimAll16Strimleft17-21.RDS")

#load metadata file
samdf <- read.csv2("C:/metadata_all.csv") #path to your metadtaa in csv format

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE pour verifier que les noms concordent

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "Field", "Organ", "Disease", "OrgHealth", "FieldOrgHealth") # indicate which columun you want to keep from your metadata
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

#save in RDS format the seqtab 

saveRDS(seqtab.nochim, "seqtab.nochim.RDS")
saveRDS(taxa, "taxa.RDS")
write.csv(taxa, "taxa.csv")
write.csv(seqtab.nochim, "seqtab.nochim.csv")



# extract info in a single table (without metadata)foranaluysis in excel and formatting to Namco
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa6<-as.data.frame(taxa)
taxa6$ID <- rownames(taxa)
all_data5  <- merge(taxa6,seqtab6,by='ID')
write.csv(all_data5, "Oeum2024_16S.csv")


#phyloseq object without phylogenetic tree

ps_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_16S


#Add ASV numbering instead of fasta sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_16S))
names(dna) <- taxa_names(ps_16S)
ps_16S  <- merge_phyloseq(ps_16S, dna)
taxa_names(ps_16S) <- paste0("ASV", seq(ntaxa(ps_16S)))
ps_16S

write.table(taxa, "taxa16S.xls")


#remove mitochondria and chloroplast 
ps_16S <- ps_16S %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Order!="Chloroplast" | is.na(Order) ) 

#remove Eukaryota
ps_16S <- ps_16S %>% subset_taxa( Kingdom!= "Eukaryota" | is.na(Kingdom) & Class!="Chloroplast" | is.na(Class) ) 
ps_16S

#Remove sample from dataset (here control)
ps_16S <- prune_samples(sample_names(ps_16S) != "RVCTRL", ps_16S) 

#Create files for each data field
ps_RV <- subset_samples(ps_16S, Field=="Rovieng")
ps_PS <- subset_samples(ps_16S, Field=="PreakSdei")


#plot bars
plot_bar(ps_RV, fill="Phylum")


# order
desired_order <- list("S1", "S2", "S3", "S16", "S17", "S18", "S4", "S5", "S6", "S19", "S20", "S21", "S13", "S14", "S15", "S28", "S29", "S30", "S7", "S8", "S9", "S22", "S23", "S24", "S10", "S11", "S12", "S25", "S26", "S27", "S31", "S32", "S33",  "S34", "S35")


p <- plot_bar(ps_16S, fill="Phylum") 
p$data$Sample <- factor(p$data$Sample, levels = desired_order)
print(p)

#Rarefy at 10000
set.seed(1234)
ps_16S_10000F <- rarefy_even_depth(ps_16S, sample.size = 10000,
                                   rngseed = F, replace = TRUE) 
ps_16S_10000F    


#Rarefaction curves

set.seed(42)

calculate_rarefaction_curves <- function(ps_16S, measures, depths) {
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


rarefaction_curve_data <- calculate_rarefaction_curves(ps_16S, c('Observed', 'Shannon'), rep(c(1, 10, 50, 100, 200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000,17000, 20000), each = 10))

summary(rarefaction_curve_data)


rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(ps)), by.x = 'Sample', by.y = 'row.names')

#plotting rarefaction curves
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
#-----------------------------------------------------

#alpha diversity
plot_richness(ps_RVL, x="FieldOrgHealth", measures=c("Shannon", "Simpson"), color="FieldOrgHealth")
plot_richness(ps, x="condition", measures=c("Shannon", "Simpson"), color="medium")

plot_richness(ps, x="medium", measures=c("Shannon", "Simpson"), color="medium")

plot_richness(ps, x="medium", measures=c("condition", "Simpson"), color="medium") + geom_boxplot()

plot_richness(ps, x="medium", measures=c("condition", "Shannon"), color="medium") + geom_boxplot()

#Ordonate with x using sortby:
plot_richness(ps, x="FieldOrgHealth", measures=c("Field", "Shannon"), color="Field", sortby = "Shannon") + geom_boxplot()
plot_richness(ps, x="medium", measures=c("Shannon", "Simpson"), color="medium", sortby = "Shannon") + geom_boxplot()

#remove grey font
plot_richness(ps, x="medium", measures=c("Shannon", "Simpson"), color="medium", sortby = "Shannon") + geom_boxplot() + theme_bw()

#Export a table with richness scores
rich = estimate_richness(ps.rarefied)
rich

#NMDS of diversity
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps_16S, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS",distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="FieldOrgHealth", title="Bray NMDS") + theme_bw()

# PCoA plot using the unweighted UniFrac as distance (need tree in phyloseq object)
wunifrac_dist = phyloseq::distance(ps_16S, method="unifrac", weighted=F)
ordination = ordinate(ps, method="PCoA", distance=wunifrac_dist)
plot_ordination(ps_RVL, ordination, color="Field") + theme(aspect.ratio=1)
#Test if conditions differ significantly from each other using the permutational ANOVA (PERMANOVA) analysis:

adonis(wunifrac_dist ~ sample_data(ps_16S_5000F)$condition)

#Plot Bars (library "microbial")
library("microbial")

plotbar(ps_16S_10000F_filter10, level = "Phylum", color = NULL,group = "FieldOrgHealth",
        top = 10,return = FALSE,fontsize.x = 10,  fontsize.y = 12)

plotbar(ps_16S_10000F_filter10, level = "Class", color = NULL,group = "FieldOrgHealth",
        top = 50,return = FALSE,fontsize.x = 10,  fontsize.y = 12)


desired_order <- list("TSA10", "TSA50", "FR", "NGN", "NFB", "Root")

plotbar(ps_16S_filter2, level = "Phylum", color = NULL,group = "medium",
        top = 10,return = FALSE,fontsize.x = 5,  fontsize.y = 12)


#Top 30 phylum
p1 <- plotbar(ps_PS, level = "Phylum", color = NULL,group = "FieldOrgHealth",
              top = 30,return = FALSE,fontsize.x = 5,  fontsize.y = 12)
p1$data$medium <- factor(p1$data$medium, levels = desired_order)
print(p1)

#Top25 class
p2 <- plotbar(ps_16S_10000F_filter10, level = "Class", color = NULL,group = "FieldOrgHealth",
              top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p2$data$medium <- factor(p2$data$medium, levels = desired_order)
print(p2)

#Top 25 Order
p3 <- plotbar(ps_16S_10000F_filter10, level = "Order", color = NULL,group = "FieldOrgHealth",
              top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p3$data$medium <- factor(p3$data$medium, levels = desired_order)
print(p3)


#top25 genus
p4 <- plotbar(ps_16S_10000F_filter10, level = "Genus", color = NULL,group = "FieldOrgHealth",
             top = 25,return = FALSE,fontsize.x = 10,  fontsize.y = 12)
p4$data$medium <- factor(p4$data$medium, levels = desired_order)
print(p4)

#figure plot compilation
cowplot::plot_grid(p1, p2, p3, p4,nrow = 2, align = 'v', axis = 'lr', rel_heights = c(0.1, 0.1))




