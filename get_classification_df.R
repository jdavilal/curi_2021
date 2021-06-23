

get_classification_df <- function(samples, sig.names, signatures) {
  
  names <- vector(mode = "character")
  signature.sample <- vector(mode = "character")
  
  for (i in 1:length(samples)) {
    names <- append(names, rep(sig.names[i], length(samples[[i]])))
    signature.sample <- append(signature.sample, samples[[i]])
  }
  
  df <- data.frame(signature.sample, names)
  
  #converts sample vector into a table of probabilities
  tab <- table(signature.sample)/sum(table(signature.sample))

  #converts table of probabilities to a data frame and defines column names
  sample.probs.df <- data.frame(tab)
  colnames(sample.probs.df) <- c("mutations", "frequencies")
  sample.probs.df
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #creates a vector of mutation types
  mutations <- rownames(signatures[[1]])
  
  #creates a data.frame of mutation types to join with sample probabilities
  df <- data.frame(mutations)
  df
  
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  sample.df <- left_join_NA(df, sample.probs.df, by = "mutations")
  sample.df
  
  sample.matrix<- as.matrix(sample.df$frequencies) + 0.0001
  rownames(sample.matrix)= mutations
  
  nmf_res <- fit_to_signatures(sample.matrix, do.call(cbind, signatures))
  nmf_res$signatures <-do.call(cbind, signatures)
  
  colnames(nmf_res$signatures) <- sig.names
  rownames(nmf_res$contribution) <- sig.names
  
  #prob.list <- vector(mode = "list", length = length(samples))
  
  prob.list <- list()
  
  for (i in length(samples)){
    prob.list[[i]]<-vector()
  }

  for (x in mutations) {
    probs <- extract_all_prob(x, 1, nmf_res)
    for (i in 1:length(samples)){
    prob.list[i] <- append(prob.list[i], probs[[i]])
    print(prob.list[i])
    }
  }
    
  
  return (prob.list)
}
