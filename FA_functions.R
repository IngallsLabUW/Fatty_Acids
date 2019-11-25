## Function definitions ##

# TODO --------------------------------------------------------
# Organize directory referencing a little easier? Set value ahead of time?
# Nested loop or list comprehension for renaming the uploads. Functions?
# Order function descriptions neatly at beginning of script.
# Cleanly name items in the upload, rather than manually changing it to SN, RZ, etc.
# Figure out a better way to choose columns out of MSDIAL
# Add function documentation and comments to clarify the process.
# Make function for always downloading the most up to date Ingalls lab standards.
# Figure out a way to preserve the QC parameter values.
# Fix the StandardizeVariables function

library(ggplot2)
library(rlist)
library(stringr)
library(tidyverse)
library(tidyr)
options(scipen=999)

SetHeader <- function(df) {
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

FilterUnknowns <- function(df) {
  df <- df %>%  
    filter(Metabolite.name != 'Unknown') %>%
    select(-c(Average.Rt.min., Formula, Ontology, INCHIKEY, SMILES, Isotope.tracking.parent.ID, Isotope.tracking.weight.number, 
              MS1.isotopic.spectrum, MS.MS.spectrum, Average.Mz, Post.curation.result, Fill.., Annotation.tag..VS1.0., RT.matched, 
              m.z.matched, MS.MS.matched, Manually.modified, Total.score:Fragment.presence..))
}

RemoveCsv <- function(full.filepaths) {
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}

ChangeClasses <- function(df) {
  for (i in c(10:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the Skyline output.
  #
  # Args
  #   msdial.file: Raw output file from Skyline.
  #
  # Returns
  #   run.type: list of labels identifying the run types, isolated from Replicate.Name.
  #   Options conssist of samples (smp), pooled (poo), standards (std), and blanks (blk).
  #
  run.type <- tolower(str_extract(msdial.file$Replicate.Name, "(?<=_)[^_]+(?=_)"))
  print(paste("Your runtypes are:", toString(unique(run.type))))
}

TrimWhitespace <- function (x) gsub("^\\s+|\\s+$", "", x)

StandardizeVariables <- function(df) {
  if (c("ReplicateName", "AreaValue", "MZValue", "RTValue", "SNValue") %in% colnames(df))
  {
    df <- df %>%
      rename(Replicate.Name = ReplicateName) %>%
      rename(Area.Value = AreaValue) %>%
      rename(MZ.Value = MZValue) %>%
      rename(RT.Value = RTValue) %>%
      rename(SN.Value = SNValue)
  }
  return(df)
}

IdentifyDuplicates <- function(df) {
  test <- which(duplicated(df$Compound.Name))
  duplicates <- as.data.frame(df$Compound.Name[test]) %>%
    rename(Compound.Name = 1) %>%
    arrange(Compound.Name)
  
  return(duplicates)
}


# Do we need this function?
CheckBlankMatcher <- function(blank.matcher) {
  # Takes a blank matcher file and separates any multi-value variable
  # columns into their own row.
  #
  # Args:
  #   blank.matcher: CSV entered by user to match samples with
  #   appropriate blanks.
  #
  # Returns:
  #   blank.matcher: new CSV with any duplicate values separated
  #   into their own rows.
  #
  blank.matcher <- do.call("rbind", Map(data.frame,
                                        Blank.Name = strsplit(as.character(blank.matcher$Blank.Name), ","),
                                        Replicate.Name = (blank.matcher$Replicate.Name))
  )
  blank.matcher <- blank.matcher[c(2, 1)]
  
  return(blank.matcher)
}