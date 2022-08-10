

extract_all_prob <- function(mutation, sample, nmf_res) {
  
  if (class(colnames(nmf_res$signature)) == "NULL") {
    stop("Column names for nmf_res$signature is empty")
  } else if (class(rownames(nmf_res$contribution)) == "NULL") {
    stop("Row names for nmf_res$contribution is empty")
  }
  
  # Creates empty list of length equal to amount of signatures
  prob = vector(mode = "list", length=length(colnames(nmf_res$signatures)))
  
  # initailize accumulator variable and denominator of Bayes Theorem
  i <- 1
  denominator = 0
  
  # Fills prob with prior and likelihood for each signature
  for (signature in colnames(nmf_res$signatures)) {
    # Naming each element of prob
    names(prob)[i] <- signature
    
    # Pr(Mutation|Signature) - likelihood
    a <- nmf_res$signatures[mutation,signature]/sum(nmf_res$signatures[,signature])
    
    # Pr(Signature) - prior
    b <- nmf_res$contribution[signature,sample]/sum(nmf_res$contribution[,sample])
    
    # Appending likelihood and prior into prob
    prob[[i]] <- c(a, b)
    
    # Adding product of likelihood and prior to denominator
    denominator = denominator + a * b
    
    # Accumulating
    i <- i + 1
  }
  
  # Intialize acummulator and vector to store posteriors
  i <- 1
  posteriors <- vector(mode = "list", length = length(colnames(nmf_res$signatures)))
  
  # Fills posteriors
  for (x in colnames(nmf_res$signatures)) {
    names(posteriors)[i] <- x
    # Calculates numerator of Bayes formula for each signature
    numerator <- prob[[x]][1] * prob[[x]][2]
    # Finishes Bayes calculation and inserts value into the vector posteriors
    posteriors[[i]] <- numerator / denominator
    
    # Accumulating
    i <- i + 1
  }
  
  return (posteriors)
}