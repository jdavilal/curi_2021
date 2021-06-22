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

#Creating sample to match signature distribution (vector of all mutations)

create_signature_sample_vector <- function(signature.matrix, number.of.mutations){
  #creates a vector of mutation types
  mutations <- rownames(signature.matrix)
  #creates a vector of probabilities of each mutation type
  probabilities = as.vector(signature.matrix[,1])
  
  #creates sample of desired number of total mutations to match the distribution of mutations in the signature
  signature.sample = sample(mutations, size = number.of.mutations, replace = TRUE, prob = probabilities)
  
  #returns a vector with length = number.of.mutations, each element is a mutation type (string)
  return(signature.sample)
}

#Creating sample to match signature distribution (matrix of probabilities)
create_signature_sample_matrix <- function(signature.matrix, number.of.mutations){
  #creates a vector of mutation types
  mutations <- rownames(signature.matrix)
  #creates a vector of probabilities of each mutation type
  probabilities = as.vector(signature.matrix[,1])
  #creates sample of desired number of total mutations to match the distribution of mutations in the signature (vector)
  signature.sample = sample(mutations, size = number.of.mutations, replace = TRUE, prob = probabilities)
  
  #converts sample vector into a table of probabilities
  tab <- table(signature.sample)/sum(table(signature.sample))
  #converts table of probabilities to a data frame and defines column names
  sample.probs.df <- data.frame(tab)
  colnames(sample.probs.df) <- c("mutations", "frequencies")
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #creates a data.frame of mutation types to join with sample probabilities
  df <- data.frame(mutations)
  
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  sample.df <- left_join_NA(df, sample.probs.df, by = "mutations")
  
  #creates a matrix of sample frequencies and defines row names as mutation types
  sample.matrix<- as.matrix(sample.df$frequencies)
  rownames(sample.matrix)= mutations
  #returns matrix
  return(sample.matrix)
}

#Calculating the cosine similarity between a vector of mutations and a mutational signature

signature_cosine_similarity <- function(mutations.vector, signature.matrix){
  #create a table of probabilities for the vector of mutations
  tab <- table(mutations.vector)/sum(table(mutations.vector))
  #transform table to a data frame and assign column names to mutations and frequencies
  mutations.df <- data.frame(tab)
  colnames(mutations.df) <- c("mutations", "frequencies")
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #creates a vector of the 96 mutation types
  mutations <- rownames(signature.matrix)
  #converts vector of 96 mutation types to a data frame
  mutation.types.df <- data.frame(mutations)
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  complete.mutations.df <- left_join_NA(mutation.types.df, mutations.df, by = "mutations")
  
  #converts data frame of mutation types and frequencies to a data frame and assigns rownames
  mutations.matrix <- as.matrix(complete.mutations.df$frequencies)
  rownames(mutations.matrix) = mutations
  
  #calculates and returns cosine similarity to the matrix for the signature of interest 
  cosine.similarity.matrix <- cos_sim_matrix(signature.matrix, mutations.matrix)
  return(cosine.similarity.matrix[,1])
}
