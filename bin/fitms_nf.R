#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
organ = "Breast"
genome.v  ="hg38"
args = commandArgs(trailingOnly=TRUE)

print_var <- function(v1) {
  deparse(substitute(v1))
}

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
    state_print <- print_var(state)
    clonality_print <- append(clonality_print, state_print)
  }
}

counter =0
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
    names(df)[names(df) == 'catalogue'] <- paste0(sample, '_',clonality_print[counter])
  }
  else{
    df$'placeholder' <- res$catalogue$catalogue
    names(df)[names(df) == 'placeholder'] <- paste0(sample, '_',clonality_print[counter])
    }
}

write.csv(df, paste0(sample, '_clonality_state_catalogue.csv'))


plotSubsSignatures(signature_data_matrix = df,output_file = paste0(sample, "_SNV_catalogues.pdf"))

res <-FitMS(catalogues = df,
            organ =organ, 
            exposureFilterType="giniScaledThreshold",
            useBootstrap = TRUE, 
            nboot = 200)

plotFitMS(res, 'results/')
