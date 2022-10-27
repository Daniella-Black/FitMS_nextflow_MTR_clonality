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

clonality_in <- c(all, clonal_any, clonal_NA, clonal_early, clonal_late, subclonal)
clonality_print <- list()
clonality <- list()

for(state in clonality_in){
  if(state != 'False'){
    clonality <- append(clonality, state)
    state_print <- [ i for i, a in locals().items() if a == state][0]
    clonality_print <- append(clonality_print, state_print)
  }
}

counter = -1
for (state in clonality){
  counter = counter +1
  filenameinput = state
  tab <- read.table(filenameinput, sep='\t')
  names(tab) <- tab[1,]
  tab <- tab[-1,]
  rownames(tab) <-NULL
  tab$position <- as.numeric(tab$position)
  res <- tabToSNVcatalogue(tab, genome.v)
  if(counter == 0){
    df <- data.frame(res$catalogue)
    df = df.rename(columns={'catalogue': paste0(sample, '_', clonality_print[counter)})
  }
  else
    df[paste0(sample, '_', clonality_print[counter))] <- res$catalogue
}

write.csv(df, paste0(sample, '_clonality_state_catalogue.csv'))


plotSubsSignatures(signature_data_matrix = df,output_file = paste0(sample, "_SNV_catalogues.pdf"))

res <-FitMS(catalogues = df,
            organ =organ, 
            exposureFilterType="giniScaledThreshold",
            useBootstrap = TRUE, 
            nboot = 200)

plotFitMS(res, 'results/')
