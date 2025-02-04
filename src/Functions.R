## Function definitions ##

SetHeader <- function(df) {
  # Remove empty or unnecessary lines from machine output, and make column names headers.
  #
  # Args
  #   df: Raw output file from MSDial.
  #
  # Returns
  #   df: modified dataframe with correct headers and no empty lines.
  #
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given list of filepaths.
  #
  # Args
  #   Character strings of filepaths in a directory.
  #
  # Returns
  #   Character strings of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*", "", no.path)
  
  return(no.path)
}

ChangeClasses <- function(df) {
  for (i in c(12:ncol(df))) {
    df[, i] <- as.numeric(as.character(df[, i]))
  }
  return(df)
}

IdentifyRunTypes <- function(msdial.file) {
  # Identify run typfes and return each unique value present in the machine output.
  #
  # Args
  #   msdial.file: Raw output file from MSDial.
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

RearrangeDatasets <- function(df, parameter) {
  df <- df %>%
    tidyr::gather(
      key = "Replicate.Name",
      value = "parameter",
      starts_with("X")) %>%
    select(Replicate.Name, parameter, everything())
  
  names(df)[2] <- parameter
  
  return(df)
}

FindStdDev <- function(df) {
  df.first <- df %>%
    group_by(Metabolite.Name) %>%
    group_split()
  df.midframe <- lapply(df.first, function(x) mutate(x, Std.dev = sd(RT.Value, na.rm = TRUE)))
  df.final <- bind_rows(df.midframe)
  
  return(df.final)
}

StandardizeMetabolites <- function(df) {
  df.standardized <- df %>%
    mutate(Metabolite.Name = ifelse(str_detect(Metabolite.Name, "Ingalls_"), sapply(strsplit(Metabolite.Name, "_"), `[`, 2), Metabolite.Name)) 
  
  df.standardized$Replicate.Name <- gsub("^.{0,1}", "", df.standardized$Replicate.Name)
  
  return(df.standardized)
}
