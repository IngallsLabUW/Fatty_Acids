## Function definitions ##

# Import original output scripts from MSDial.
# Unify variable names and classes.

# TODO --------------------------------------------------------
# Organize directory referencing a little easier? Set value ahead of time?
# Nested loop or list comprehension for renaming the uploads. Functions?
# Order function descriptions neatly at beginning of script.
# Cleanly name items in the upload, rather than manually changing it to SN, RZ, etc.

library(tidyverse)
library(tidyr)
# options(scipen=999)

#FA_env <- new.env()


SetHeader <- function(df) {
  df <- df[!(is.na(df[1]) | df[1]==""), ]
  colnames(df) <- make.names(as.character(unlist(df[1,])))
  df <- df[-1, ]
  
  return(df)
}

FilterUnknowns <- function(df) {
  df <- df %>%  
    filter(Metabolite.name != 'Unknown') %>%
    select(-c(Average.Rt.min., Formula, Ontology, INCHIKEY, SMILES, Isotope.tracking.parent.ID, Isotope.tracking.weight.number, MS1.isotopic.spectrum, MS.MS.spectrum,
              Average.Mz))
  }

RemoveCsv <- function(full.filepaths) {
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}


# Processing begins here --------------------------------------------------

filenames <- RemoveCsv(list.files(path = 'data', pattern = '*.csv'))

for(i in filenames) {
  filepath <- file.path('data', paste(i,".csv", sep = ""))
  assign(i, read.csv(filepath))
}

my.dataframes <- lapply(mget(ls(pattern = "_")), dim)

## Trim, set header, filter unknowns
# TODO (rlionheart): Make a loop or function
SN.pos <- header.true(SN_HILICPos_EddyTransect)
SN.pos <- filter.unknowns(SN.pos)

RT.pos <- header.true(RT_HILICPos_EddyTransect)
RT.pos <- filter.unknowns(RT.pos)

Area.pos <- header.true(Area_HILICPos_EddyTransect)
Area.pos <- filter.unknowns(Area.pos)

MZ.pos <- header.true(Mz_HILICPos_EddyTransect)
MZ.pos <- filter.unknowns(MZ.pos)










