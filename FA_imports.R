source("FA_Functions.R")

# Fill out user data --------------------------------------------------
# The matching variable should be included in every MSDial file in the data_raw folder
matching.variable <- "FA"

# Ensure correct columns are dropped from imports.
columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')


# Import all MSDial files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

# Set header, filter unknowns ---------------------------------------
runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------

# Positive
SN <- SN_KM1906_FA_DepthProfiles %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "SN.Value",
    starts_with("X")) %>%
  select(Replicate.Name, SN.Value, everything())

RT <- RT_KM1906_FA_DepthProfiles %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "RT.Value",
    starts_with("X")) %>%
  select(Replicate.Name, RT.Value, everything())

Area <- Area_KM1906_FA_DepthProfiles %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "Area.Value",
    starts_with("X")) %>%
  select(Replicate.Name, Area.Value, everything())

MZ <- Mz_KM1906_FA_DepthProfiles %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "MZ.Value",
    starts_with("X")) %>%
  select(Replicate.Name, MZ.Value, everything())


# Combine to one dataset --------------------------------------------------
combined <- Area %>%
  left_join(MZ) %>%
  left_join(SN) %>%
  left_join(RT) %>%
  select(Replicate.Name, Area.Value, MZ.Value, RT.Value, SN.Value, everything())

combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)

rm(list = ls()[!ls() %in% c("combined", lsf.str())])