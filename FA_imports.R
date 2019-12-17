source("FA_Functions.R")

# Fill out user data --------------------------------------------------
# The matching variable should be included in every MSDial file in the data_raw folder
matching.variable <- "FA"

# Ensure correct columns are dropped from imports.
columns.to.drop <- c('Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')
# 'Average.Rt.min.', 'RT.similarity'

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
  headers.set[[df]] <- headers.set[[df]] %>% rename(Metabolite.Name = Metabolite.name) 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!str_detect(Metabolite.Name, "w/o MS2|Unknown"))
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------
Area <- RearrangeDatasets(Area_KM1906_FA_DepthProfiles, parameter = "Area.Value")
Mz   <- RearrangeDatasets(Mz_KM1906_FA_DepthProfiles, parameter = "Mz.Value")
RT   <- RearrangeDatasets(RT_KM1906_FA_DepthProfiles, parameter = "RT.Value")
SN   <- RearrangeDatasets(SN_KM1906_FA_DepthProfiles, parameter = "SN.Value")


# Combine to one dataset --------------------------------------------------
combined <- Area %>%
  left_join(Mz) %>%
  left_join(SN) %>%
  left_join(RT) %>%
  select(Replicate.Name, Area.Value, Mz.Value, RT.Value, SN.Value, everything())

# Standardize dataset --------------------------------------------------
combined.final <- StandardizeMetabolites(combined) %>%
  filter(!str_detect(Replicate.Name, "Blk_Blk_20"))

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_FA_combined_", currentDate, ".csv", sep = "")

write.csv(combined.final, csvFileName, row.names = FALSE)

rm(list = ls())