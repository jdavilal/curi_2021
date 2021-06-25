
#Arguments:
#1) list where each element is a sample generated from a certain signature and
#   each subelement is a trinucleotide mutation type (as a string)
#2) vector of strings with names of signatures in the same order as the elements
#   in the list from argument 1
#3) list where each element is the true mutational signature vector


get_classification_df <- function(samples, sig.names, signatures) {
  
  #Empty vectors for names of mutation sources (signatures) and randomly generated samples
  names <- vector(mode = "character")
  signature.sample <- vector(mode = "character")
  
  for (i in 1:length(samples)) {
    #Adds each signature name for every element/profile in the sample to empty vector
    names <- append(names, rep(sig.names[i], length(samples[[i]])))
    #Adds each element/profile in the sample to empty vector
    signature.sample <- append(signature.sample, samples[[i]])
  }
  
  #Puts the long generated mutations and their source into a data frame and names columns
  df.final <- data.frame(signature.sample, names)
  colnames(df.final) <- c("mutations", "truth")
  
  
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
  
  #creates a vector of mutation types
  mutations <- rownames(signatures[[1]])
  
  #creates a data.frame of mutation types to join with sample probabilities
  df <- data.frame(mutations)

  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  sample.df <- left_join_NA(df, sample.probs.df, by = "mutations")

  #Add pseudocount to column of probabilities 
  sample.matrix<- as.matrix(sample.df$frequencies) + 0.0001
  rownames(sample.matrix)= mutations
  
  #Computes exposure matrix using sample matrix and true signature matrix
  nmf_res <- fit_to_signatures(sample.matrix, do.call(cbind, signatures))
  #Add element with original signatures so that Bayes function accepts this object
  nmf_res$signatures <-do.call(cbind, signatures)
  
  #name the rows and columns that represent mutational signatures for compatibility with Bayes function
  colnames(nmf_res$signatures) <- sig.names
  rownames(nmf_res$contribution) <- sig.names
  
  
  #empty list to hold probability of each signature given randomly generated mutation
  prob.list <- list()
  
  for (i in length(samples)){
    #each element of empty list is vector of probabilities for every element/profile from mutational sample
    prob.list[[i]]<-vector()
  }

  for (m in mutations) {
    #each mutation type (string) is sent to Bayes function along with "nmf" results
    probs <- extract_all_prob(m, 1, nmf_res)
    for (n in 1:length(samples)){
      #probability for given mutation is sent to all corresponding signatures then proceeds to next mutation
      prob.list[[n]][m] <- probs[[n]]
    }
  }
  
  #convert posterior probability list into a dataframe
  prob.list.df = data.frame(prob.list)
  rownames(prob.list.df) <- NULL
  #each column represents a signature
  colnames(prob.list.df) <- sig.names
  #add column with mutation names
  prob.list.df$mutations <- mutations
  
  #join posterior probabilities with randomized mutations+source (keep only probabilities that
  #correspond to an existing randomized mutation)
  df.final <- left_join(df.final, prob.list.df, by = "mutations")
  
  #empty vector to hold Bayesian classifier
  classify <- vector(mode = "character")
  
  for (row in 1:nrow(df.final)){
    #store posterior probabilities in each row for all signatures
    probabilities <- df.final[row,3:(2+length(samples))]
    #finds the index for the maximum probability and retrieves name of signature that is the max
    classifier <- colnames(probabilities)[which.max(probabilities[1,])]
    #each element is signature name (string) that has the maximum posterior probability
    classify[row]<- classifier
  }
  
  #add column to the data frame that represents Bayesian classifier
  df.final$classify <- classify
  
  #add column to the data frame that determines if classifier is different from source
  df.final <- df.final %>%
    #misclassification when the source string does not match the Bayesian classifier
    mutate(misclassification = ifelse(classify != truth, 1, 0))
  

  return (df.final)

}
