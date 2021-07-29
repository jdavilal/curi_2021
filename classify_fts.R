get_classification_df <- function(mutations_vector, signature_labels, signatures_matrix) {
  
  #converts mutations_vector into a table of probabilities
  tab <- table(mutations_vector)/sum(table(mutations_vector))
  
  #converts table of probabilities to a data frame and defines column names
  mutation.probs.df <- data.frame(tab)
  colnames(mutation.probs.df) <- c("mutations", "frequencies")
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #creates a vector of all 96 snv mutation types
  mutations <- rownames(signatures[[1]])
  
  #creates a data.frame of mutation types to join with sample probabilities
  df <- data.frame(mutations)
  
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  mut.profile.df <- left_join_NA(df, mutation.probs.df, by = "mutations")
  
  #Transform mut.profile.df to mutational profile matrix, add pseudocount to column of probabilities 
  mut.profile<- as.matrix(mut.profile.df$frequencies) + 0.0001
  rownames(mut.profile)= mutations
  
  fts.res <- fit_to_signatures(mut.profile, signatures.matrix)
  
}