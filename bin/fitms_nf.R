#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
organ = "Breast"
genome.v  ="hg38"
args = commandArgs(trailingOnly=TRUE)


sample <- args[1]
all <- args[2]
clonal_any <- args[3]
clonal_NA <- args[4]
clonal_early <- args[5]
clonal_late <- args[6]
subclonal <- args[7]

clonality_list <- c('sample', 'all', 'clonal_any', 'clonal_NA', 'clonal_early', 'clonal_late', 'subclonal')

for (state in 2:7){
  filenameinput <- args[state]
  tab <- read.table(filenameinput, sep='\t')
  if(nrow(tab) >3){
    names(tab) <- tab[1,]
    tab <- tab[-1,]
    rownames(tab) <-NULL
    tab$position <- as.numeric(tab$position)
    res <- tabToSNVcatalogue(tab, genome.v)
    if(state == 2){
       df <- data.frame(res$catalogue)
       names(df)[names(df) == 'catalogue'] <- paste0(sample, '_',clonality_list[state])
        }
    else{
      df$'placeholder' <- res$catalogue$catalogue
      names(df)[names(df) == 'placeholder'] <- paste0(sample, '_',clonality_list[state])
      }
    }
}

write.csv(df, paste0(sample, '_clonality_state_catalogue.csv'))


plotSubsSignatures(signature_data_matrix = df,output_file = paste0(sample, "_SNV_catalogues.pdf"))

res <-FitMS(catalogues = df,
           organ =organ, 
           #27/10/2022 run used gini, now 09/01/2022 trying without
           #exposureFilterType="giniScaledThreshold",
           useBootstrap = TRUE, 
           nboot = 200)

plotFitMS(res, 'results/')
