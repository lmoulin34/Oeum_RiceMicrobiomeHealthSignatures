#Rscript toproduce box plots with p values of kurskal wallis + post-hoc pairwise Wilcoxon tests with results as letters 

install.packages("rstatix")
install.packages("ggpubr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("multcompView")
install.packages("emmeans")
install.packages("multcomp")

library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(multcompView)
library(emmeans)
library(multcomp)

#Import Table
setwd('C:/Users/moulinl/Documents/Thèse_Kakada_OEUM/Plant_tests/Xoo_symptoms_2024')

data<-read.csv2('BiocontrolPhkar2024.csv')


# Indiquer quel est la colonne facteur (ici c'est la colonne Treatment)
data$treatment <- as.factor(data$treatment)
data$Symptom <- as.numeric(data$Symptom)
data$Height <- as.numeric(data$Height)
data$SDW <- as.numeric(data$SDW)


#view distribution of data
hist(data$Symptom)


#Reorder data 
#placette1
data <- data %>% reorder_levels(treatment, order = c("NoInoc", "3533", "3560", "3561", "3562"))

#Statistiques descriptives
data %>% group_by(treatment) %>% get_summary_stats (Symptom, type = "mean_sd")
ggboxplot(data, x = "treatment", y = "Symptom") +  geom_jitter(fill="black")


#non-parametric test

res.kruskal <- data %>% kruskal_test(Symptom ~ treatment)
res.kruskal

pwc <- data %>% 
  dunn_test(Symptom ~ treatment, p.adjust.method = "bonferroni") 
pwc


pwc <- pwc %>% add_xy_position(x = "treatment")


ggboxplot(data, x = "treatment", y = "Symptom") +
  theme_classic()+
  theme(axis.title=element_text(size=10),
        plot.title = element_text(size=12),
        plot.subtitle = element_text(size=10),
        legend.text=element_text(size=8),
        legend.title = element_text(size = 10),
        axis.text.y = element_text(size=8,angle=0),
        axis.text.x = element_text(size=8,angle=90,hjust=1),
        legend.position="none") +  geom_jitter(fill="black") +stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
       caption = get_pwc_label(pwc))

----------------------------------------------------------------
  #Visualiser stat sous forme de lettre
  
  #need to create a function first:
  
  tri.to.squ<-function(x)
  {
    rn <- row.names(x)
    cn <- colnames(x)
    an <- unique(c(cn,rn))
    myval <-  x[!is.na(x)]
    mymat <-  matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
    for(ext in 1:length(cn))
    {
      for(int in 1:length(rn))
      {
        if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
        mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
        mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      }
      
    }
    return(mymat)
  } 

#Statistical tests

res.kruskal <- kruskal.test(Symptom ~treatment, data=data)
res.kruskal

pp <- pairwise.wilcox.test(data$Symptom, data$treatment,p.adjust.method ="bonferroni" )
pp

#then
mymat <-tri.to.squ(pp$p.value)
mymat


myletters <- multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
myletters

myletters_df <- data.frame(treatment=names(myletters$Letters),letter = myletters$Letters )
myletters_df

ggplot(data, aes(x=treatment, y=Symptom, colour=treatment, fill=treatment))+
  geom_boxplot(outlier.alpha = 0, alpha=0.25)+
  geom_jitter(width=0.25)+  
  stat_summary(fun.y = 'mean', colour="black", geom="point", 
               shape=18, size=2) +
  theme_classic()+
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))+ 
  geom_text(data = myletters_df, aes(label = letter, y = 40 ), colour="black", size=4)


