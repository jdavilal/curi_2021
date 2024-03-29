# Script for converting mutation samples into vcf format

create_gr_from_sample <- function(sample, specie, chromosome) {
  
  # Function arguments validation
  if (type(chromosome) != "character") {
    stop("argument chromosome is not of type character")
  } else if (str_sub(chromosome, 1, 3) != "chr" || suppressWarnings(is.na(as.numeric(str_sub(chromosome, 4))))) {
    stop("argument chromosome is not of format: chr***. e.g chr17")
  } else if (!is.data.frame(sample)) {
    stop("argument sample is not a dataframe")
  } else if (class(specie) != "BSgenome") {
    stop("argument specie is not a Bsgenome object")
  }
  
  # Instantiate vector objects for function output
  seq <- getSeq(specie, chromosome)
  position <- vector(mode = "double")
  id <- vector(mode = "character")
  info <- vector(mode = "character")
  sample1 <- vector(mode = "character")
  sample2 <- vector(mode = "character")
  # Original nucleotide
  ref <- vector(mode = "character")
  # New nucleotide
  alt <- vector(mode = "character")
  
  print("Determining position of mutations...")
  
  # Iteration begins at 10000 to avoid telomeres
  i <- 10000
  
  # Loops through every mutation in sample
  for (m in 1:length(sample$mutations)) {
    
    # Reformat mutation string and stores substitutions into ref and alt
    mut.reformatted <- str_replace_all(sample$mutations[m], "[[:punct:]>]", "")
    ref <- c(ref, str_sub(mut.reformatted, 2, 2))
    alt <- c(alt, str_sub(mut.reformatted, 3, 3))
    str_sub(mut.reformatted, 3, 3) <- ""
    
    info <- c(info, paste("SOMATIC", ";TRUTH=", sample$truth[m], ";PROB", "=", round(sample[m, 3], 3), ",", round(sample[m, 4], 3), sep = ""))
    
    # Searches chromosome sequence for mutation match
    while (mut.reformatted != toString(seq[i:(i+2)])) {
      i = i + 1
    }
    
    # Once matched, store position
    position <- c(position, i+1)
    id <- c(id, paste(str_sub(chromosome, 4, 4), ":", i+1, "_", str_sub(mut.reformatted, 1, 1), "/", str_sub(mut.reformatted, 3, 3), sep=""))
    # if (sample$truth[m] == "FFPE") {
    #   sample1 <- c(sample1, paste("0/1:", sample(10:50, 1), sep = ""))
    #   sample2 <- c(sample2, paste("0/1:", sample(50:100, 1), sep = ""))
    # } else {
    #   sample1 <- c(sample1, paste("0/1:", sample(51:100, 1), sep = ""))
    #   sample2 <- c(sample2, paste("0/1:", sample(10:50, 1, replace=T), sep = ""))
    # }
    
    sample1 <- c(sample1, paste("0/1:", sample(30:60, 1), sep = ""))
    sample2 <- c(sample2, paste("0/1:", sample(30:60, 1), sep = ""))
  }
  
  print("...Done!")
  
  # Creates GRange object with position, ref, and alt
  gr <- GRanges(Rle(c("chr1"), c(1)), IRanges(start = position, end = position))
  
  # Additional metacolumns added
  gr$ID <- "."
  gr$REF <- ref
  gr$ALT <- alt
  gr$QUAL <- 100
  gr$FILTER <- "."
  gr$INFO <- info
  gr$FORMAT <- "GT:GQ"
  gr$SAMPLE1 <- sample1
  gr$SAMPLE2 <- sample2
  
  return (gr)
}

write_grange_to_vcf <- function(grange, file.name) {
  
  # Create empty file
  file.data <- ""
  write(file.data, file.name)
  
  # Loads in template for a vcf file and writes to new file
  file.data <- read_file("../vcf-template.vcf")
  file.data <- str_remove_all(file.data, "\r")
  
  # Appends to file.data meta column names
  for (name in colnames(mcols(grange))) {
    file.data <- paste(file.data,"\t", name, sep="")
  }
  
  write(file.data, file.name)
  
  # Loops through rows in grange
  for (x in 1:length(grange)) {
    
    # Append CHROM number and position to file.data
    file.data <- paste(file.data, "\n","1\t", start(ranges(grange)[x]), sep="")
  
    # Loops through each metacolumn and appends meta data to file.data
    for (col in 1:length(colnames(mcols(grange)))) {
      meta.data <- mcols(grange)[colnames(mcols(grange))[col]][[1]][x]
      file.data <- paste(file.data,"\t", meta.data, sep="")
    }
  }
  
  # Loads file.data to new file
  write(file.data, file.name)
}

# Shortcut function for turning sample to vcf file
write_sample_to_vcf <- function(sample, file.name, specie, chromosome) {
  
  # Creates grange object from sample
  gr <- create_gr_from_sample(sample, specie, chromosome)
  print(paste("Writing to", file.name))
  
  # Writes grange to a vcf file
  write_grange_to_vcf(gr, file.name)
  
  return (gr)
}

# For converting a data frame to a grange object
convert_df_to_gr <- function(df, specie, chromosome) {

  # Creates a grange object with the mutations from the data frame
  gr <- create_gr_from_sample(df$mutations, specie, chromosome)
  
  # Get amount of signature present in sample
  # sig.amount <- df %>%
  #   group_by(truth) %>%
  #   summarise() %>%
  #   count()
  # 
  # # Add metacolumns to grange object
  # for (i in 3:(sig.amount$n + 2)) {
  #   values(gr)[colnames(df)[i]] <- df[i]
  # }
  
  for (x in 1:length(values(gr))) {
    
  }
  
  return (gr)
}


write_df_to_vcf <- function(df, specie, chromosome, file.name) {
  gr <- convert_df_to_gr(df, specie, chromosome)
  
  write_grange_to_vcf(gr, file.name)
}

get_mut_pos <- function(mutation) {
  
}

create_classify_df <- function(vcf.file, sequence) {
  
  file.name <- str_remove(vcf.file, ".vcf")
  files <- c(paste(file.name, "_Refined.vcf", sep = ""), paste(file.name, "_Artifact.vcf", sep = ""))
  
  mutations.vector <- vector(mode = "character")
  truth.vector <- vector(mode = "character")
  classify.vector <- vector(mode = "character")
  misclass.vector <- vector(mode = "double")
  
  for (file in 1:2) {
    sample <- ReadVCF(files[file])
    
    positions <- sample[1]$data$POS
    truths <- str_remove_all(sample[1]$data$INFO, "SOMATIC;TRUTH=")
    
    for (i in 1:length(positions)) {
      mutation <- paste(sequence[positions[i] - 1], "[", sequence[positions[i]], ">", sample[1]$data[i]$ALT, "]", sequence[positions[i] + 1], sep = "")
      mutations.vector <- c(mutations.vector, mutation)
      
      if (str_split(truths, ";")[[i]][1] != "FFPE") {
        other.sig.name <- str_split(truths, ";")[[i]][1]
      }
      
      truth.vector <- c(truth.vector, str_split(truths, ";")[[i]][1])
      
      if (file == 1) {
        classify.vector <- c(classify.vector, other.sig.name)
      } else {
        classify.vector <- c(classify.vector, "FFPE")
      }
    }
    
    
  }
  
  df <- data.frame(mutations = mutations.vector, truth = truth.vector, classify = classify.vector)
  
  df <- df %>%
    mutate(misclassification = ifelse(truth != classify, 1, 0))

  return (df)
}





