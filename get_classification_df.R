

get_classification_df <- function(samples, sig.names, signatures) {
  
  names <- vector(mode = "character")
  mutation <- vector(mode = "character")
  
  for (i in 1:length(samples)) {
    names <- append(names, rep(sig.names[i], length(samples[[i]])))
    mutation <- append(mutation, samples[[i]])
  }
  
  df <- data.frame(mutation, names)
  
  return (df)
}