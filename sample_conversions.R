# Script for converting mutation samples into vcf format

create_gr_from_sample <- function(sample, specie, chromosome) {
  
  # Loads chromosome sequence
  seq <- getSeq(specie, chromosome)
  
  # Initiate accummulator
  i <- 10000
  position <- vector(mode = "double")
  
  # ID
  id <- vector(mode = "character")
  
  # Original nucleotide
  ref <- vector(mode = "character")
  # New nucleotide
  alt <- vector(mode = "character")
  
  print("Determining position of mutations...")
  
  # Loops through every mutation in sample
  for (mutation in sample) {
    
    # Reformat mutation string and stores substitutions into ref and alt
    mut.reformatted <- str_replace_all(mutation, "[[:punct:]>]", "")
    ref <- c(ref, str_sub(mut.reformatted, 2, 2))
    alt <- c(alt, str_sub(mut.reformatted, 3, 3))
    str_sub(mut.reformatted, 3, 3) <- ""
    
    # Searches chromosome sequence for mutation match
    while (mut.reformatted != toString(seq[i:(i+2)])) {
      i = i + 1
    }
    
    # Once matched, store position
    position <- c(position, i+1)
    id <- c(id, paste(str_sub(chromosome, 4, 4), ":", i+1, "_", str_sub(mut.reformatted, 2, 2), "/", str_sub(mut.reformatted, 3, 3), sep=""))
  }
  
  print("...Done!")
  
  # Creates GRange object with position, ref, and alt
  gr <- GRanges(Rle(c("chr1"), c(1)), IRanges(start = position, end = position))
  
  gr$ID <- id
  gr$REF <- ref
  gr$ALT <- alt
  gr$QUAL <- 0
  gr$FILTER <- "PASS"
  gr$INFO <- "."
  
  return (gr)
}

write_grange_to_vcf <- function(grange, file.name) {
  
  # Create empty string
  file.data <- ""
  
  # Create empty file
  write(file.data, file.name)
  
  # Loads in template for a vcf file and writes to new file
  file.data <- read_file("vcf-template.vcf")
  
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
