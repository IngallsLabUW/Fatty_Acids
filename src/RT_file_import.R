# Retention time prediction 

pattern = "combined"

# Import QC'd files ----------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_intermediate/", pattern = pattern))
filepath <- file.path("data_intermediate", paste(filename, ".csv", sep = ""))

fatty.acids <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

# Create file isolating runs that have comments.
file.comments <- fatty.acids %>%
  select(Replicate.Name, Metabolite.Name, Comment) %>%
  filter(!nchar(Comment) < 5)

# Join expected fatty acid file with data.  
FA.expected <- fatty.acids %>%
  left_join(FA.expected) %>%
  mutate(Cmpd.with.Std = ifelse(str_detect(Metabolite.Name, "_Std"), "Standard.Compound", 
                                ifelse(str_detect(Metabolite.Name, "IS"), "Internal.Standard", "NonStandard.Compound"))) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE),
         Min.RT.Value = min(RT.Value, na.rm = TRUE),
         Max.RT.Value = max(RT.Value, na.rm = TRUE)) %>%
  select(Replicate.Name, Metabolite.Name, RT.Value, RT.Expected, Mean.RT.Value:Max.RT.Value, Cmpd.with.Std) 

rm(list = setdiff(ls()[!ls() %in% c("FA.expected")], lsf.str()))