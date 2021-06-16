library(readr)

# Loads the published COSMIC mutational signatures in a format that can be
# used by MutationalPatterns
load_cosmic_matrix<- function(mut.order="supplied_results/mut_sig.order.txt") {
  
  sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                  "signatures_probabilities.txt", sep = "")
  
  cancer_signatures <-  read.table(sp_url, sep = "\t", header = TRUE)
  # Match the order of the mutation types to MutationalPatterns standard
  mut_mat <-  read_tsv(mut.order)
  new_order <-  match(mut_mat$ext.context, cancer_signatures$Somatic.Mutation.Type) 
  # Reorder cancer signatures dataframe
  cancer_signatures <-  cancer_signatures[as.vector(new_order),]
  # Add trinucletiode changes names as row.names
  row.names(cancer_signatures) <-  cancer_signatures$Somatic.Mutation.Type
  # Keep only 96 contributions of the signatures in matrix
  cancer_signatures <-  as.matrix(cancer_signatures[,4:33])
  return(cancer_signatures)
}

# Prop = [0.1, 0.9] vector
# Signatures = [21, 7] vector

create_test_sample <- function(signatures, prop, num) {
  cosmic_signatures <- load_cosmic_matrix()
  process <- cosmic_signatures[, signatures]
  
  exposures <- matrix(prop, nrow = length(signatures), ncol = num)
  
  return ((process %*% exposures) * 10^5)
}

cosmic_signatures[, "Signature.21"]
