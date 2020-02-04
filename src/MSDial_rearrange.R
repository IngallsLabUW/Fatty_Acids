# Set header, filter unknowns ---------------------------------------
columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 'SMILES', 
                     'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 'RT.similarity', 
                     'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

runs <- grep(matching.pattern, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))

headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) {
  headers.set[[df]] <- headers.set[[df]] %>% filter(!str_detect(Metabolite.name, "Unknown|w/o"))
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
  headers.set[[df]] <- headers.set[[df]] %>% rename(Metabolite.Name = Metabolite.name)
}

# Change variable classes ---------------------------------------------------------------------
classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs

list2env(classes.changed, globalenv())


# Rearrange data and combine to one dataframe -------------------------------------------------
  
# Cyano
Area <- RearrangeDatasets(Area.FA.data, parameter = "Area.Value")
Mz   <- RearrangeDatasets(Mz.FA.data, parameter = "Mz.Value")
RT   <- RearrangeDatasets(RT.FA.data, parameter = "RT.Value")
SN   <- RearrangeDatasets(SN.FA.data, parameter = "SN.Value")

# Combine to one dataset
combined.final <- Area %>%
  left_join(Mz) %>%
  left_join(SN) %>%
  left_join(RT) %>%
  select(Replicate.Name, Area.Value, Mz.Value, RT.Value, SN.Value, everything())


# Standardize dataset --------------------------------------------------
combined.final <- StandardizeMetabolites(combined.final)

currentDate <- Sys.Date()
csvFileName <- paste("data_intermediate/MSDial_combined_", file.pattern, "_", currentDate, ".csv", sep = "")

write.csv(combined.final, csvFileName, row.names = FALSE)

rm(list = setdiff(ls()[!ls() %in% c("file.pattern")], lsf.str()))