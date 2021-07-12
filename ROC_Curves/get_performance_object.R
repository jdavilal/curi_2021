#ffpe.signature and cosmic.v3 must be loaded ahead of time
get_performance_object <- function(cosmic_name, total_mut, prop_ffpe){
  cosmic.sig.mat <- as.matrix(cosmic.v3[,cosmic_name])
  rownames(cosmic.sig.mat) <- mutations
  cosmic.mut.no = total_mut-(total_mut*prop_ffpe)
  ffpe.mut.no = total_mut*prop_ffpe
  
  cosmic.samp <- create_signature_sample_vector_seed(cosmic.sig.mat, cosmic.mut.no)
  ffpe.samp <- create_signature_sample_vector_seed(ffpe.signature, ffpe.mut.no)
  
  class.df.1 <- get_classification_df(list(cosmic.samp, ffpe.samp), c(cosmic_name, "FFPE"),
                                      list(cosmic.sig.mat, ffpe.signature))
  
  class.df.2 <- class.df.1 %>%
    mutate(indicator = ifelse(truth == cosmic_name, 1, 0))
  
  pred <- prediction(class.df.2[3], class.df.2[7])
  perf <- performance(pred, "tpr", "fpr")
  auc.perf <- performance(pred, measure = "auc")
  auc.val <- auc.perf@y.values[[1]][1]
  acc.perf <- performance(pred, measure = "acc")
  #find accuracy for 0.5 threshold/cutoff
  #get line segment using threshold immediately before and after 0.5
  low.ind = max(which(slot(acc.perf, "x.values")[[1]] >=0.5))
  high.ind = min(which(slot(acc.perf, "x.values")[[1]] <=0.5))
  #point slop formula of a line
  slope = (slot(acc.perf, "y.values")[[1]][high.ind] - 
             slot(acc.perf, "y.values")[[1]][low.ind])/
    (slot(acc.perf, "x.values")[[1]][high.ind] - slot(acc.perf, "x.values")[[1]][low.ind])
  acc.50 = slope*(0.5 - slot(acc.perf, "x.values")[[1]][low.ind]) + 
    slot(acc.perf, "y.values")[[1]][low.ind]
  #find the cutoff that makes the max accuracy
  ind = which.max( slot(acc.perf, "y.values")[[1]] )
  acc.max = slot(acc.perf, "y.values")[[1]][ind]
  cutoff = slot(acc.perf, "x.values")[[1]][ind]
  
  return(list(rocplot = perf, auc = auc.val, accplot = acc.perf, 
              acc50 = acc.50, acc.max = acc.max, cutoff = cutoff))
}

vector_to_matrix <- function(mutations.vector){
  tab <- table(mutations.vector)/sum(table(mutations.vector))
  #transform table to a data frame and assign column names to mutations and frequencies
  mutations.df <- data.frame(tab)
  colnames(mutations.df) <- c("mutations", "frequencies")
  
  #defines a function that joins two data frames and replaces NA with 0
  left_join_NA <- function(x, y, by) {
    left_join(x = x, y = y, by = by) %>% 
      mutate_each(funs(replace(., which(is.na(.)), 0)))
  }
  
  #converts vector of 96 mutation types to a data frame
  mutation.types.df <- data.frame(mutations)
  #joins together sample probabilities with data frame of mutation types so that all 96 mutation types are included
  complete.mutations.df <- left_join_NA(mutation.types.df, mutations.df, by = "mutations")
  
  #converts data frame of mutation types and frequencies to a data frame and assigns rownames
  mutations.matrix <- as.matrix(complete.mutations.df$frequencies)
  rownames(mutations.matrix) = mutations
  
  return(mutations.matrix)
}

get_profile_matrices <- function(cosmic_name, total_mut, prop_ffpe) {
  cosmic.sig.mat <- as.matrix(cosmic.v3[,cosmic_name])
  rownames(cosmic.sig.mat) <- mutations
  cosmic.mut.no = total_mut-(total_mut*prop_ffpe)
  ffpe.mut.no = total_mut*prop_ffpe
  
  cosmic.samp <- create_signature_sample_vector_seed(cosmic.sig.mat, cosmic.mut.no)
  ffpe.samp <- create_signature_sample_vector_seed(ffpe.signature, ffpe.mut.no)
  
  cosmic.matrix <- vector_to_matrix(cosmic.samp)
  colnames(cosmic.matrix) <- c("Original")
  cosmic.ffpe.matrix <- vector_to_matrix(c(cosmic.samp, ffpe.samp))
  colnames(cosmic.ffpe.matrix) <- c("FFPE Added")
  
  class.df.1 <- get_classification_df(list(cosmic.samp, ffpe.samp), c(cosmic_name, "FFPE"),
                                      list(cosmic.sig.mat, ffpe.signature))
  
  class.df.2 <- class.df.1 %>%
    filter(classify != "FFPE")
  
  reconstructed.vector <- class.df.2$mutations
  reconstructed.matrix <- vector_to_matrix(reconstructed.vector)
  colnames(reconstructed.matrix) <- c("Reconstructed")
  
  cos.sim.mat <- cos_sim_matrix(cosmic.matrix, reconstructed.matrix)
  cos.sim <- cos.sim.mat[1,1]
  
  return(list(cosmic.sig.mat, cosmic.matrix, cosmic.ffpe.matrix, reconstructed.matrix, cos.sim))
}




