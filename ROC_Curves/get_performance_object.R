#ffpe.signature and cosmic.v3 must be loaded ahead of time
get_performance_object <- function(cosmic_name, total_mut, prop_ffpe){
  cosmic.sig.mat <- as.matrix(cosmic.v3[,cosmic_name])
  rownames(cosmic.sig.mat) <- mutations
  cosmic.mut.no = total_mut-(total_mut*prop_ffpe)
  ffpe.mut.no = total_mut*prop_ffpe
  
  cosmic.samp <- create_signature_sample_vector(cosmic.sig.mat, cosmic.mut.no)
  ffpe.samp <- create_signature_sample_vector(ffpe.signature, ffpe.mut.no)
  
  class.df.1 <- get_classification_df(list(cosmic.samp, ffpe.samp), c(cosmic_name, "FFPE"),
                                      list(cosmic.sig.mat, ffpe.signature))
  
  class.df.2 <- class.df.1 %>%
    mutate(indicator = ifelse(truth == cosmic_name, 1, 0))
  
  pred <- prediction(class.df.2[3], class.df.2[7])
  perf <- performance(pred, "tpr", "fpr")
  
  return(perf)
}

