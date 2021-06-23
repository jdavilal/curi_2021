

get_classification_df <- function(samples, sig.names, signatures) {
  
  names <- vector(mode = "character")
  signature.sample <- vector(mode = "character")
  
  for (i in 1:length(samples)) {
    names <- append(names, rep(sig.names[i], length(samples[[i]])))
    signature.sample <- append(signature.sample, samples[[i]])
  }
  
  df.final <- data.frame(signature.sample, names)
  colnames(df.final) <- c("mutations", "truth")
  
  
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
  
  
  prob.list <- list()
  
  for (i in length(samples)){
    prob.list[[i]]<-vector()
  }

  for (m in mutations) {
    probs <- extract_all_prob(m, 1, nmf_res)
    for (n in 1:length(samples)){
    prob.list[[n]][m] <- probs[[n]]
    }
  }
  
  prob.list.df = data.frame(prob.list)
  rownames(prob.list.df) <- NULL
  colnames(prob.list.df) <- sig.names
  
  prob.list.df$mutations <- mutations
  
  df.final <- left_join(df.final, prob.list.df, by = "mutations")
  
  classify <- vector(mode = "character")
  for (row in df.final){
    x <- colnames(df.final[which.max(df.final[row,])])
    append(classify, x)
  }
  
 # df.final <- df.final %>%
    #group_by(mutations) %>%
    #mutate(classify = colnames(df.final)[which.max()]) %>%
    #ungroup()
  
  #colnames(df.final)[which.max(df.final[1,])]
  
  #classify.df <- df.final %>%
    #pivot_longer(sig.names, names_to = "process", values_to = "probability") %>%
    #group_by(mutations) %>%
    #filter(probability == max(probability)) %>%
    #mutate(classify = process) %>%
    #select(mutations, classify)
  
  #df.final <- df.final %>%
    #left_join(classify.df, by = "mutations")
  
  return (classify)

}
