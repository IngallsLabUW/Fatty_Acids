# Parameter assignment ----------------------------------------------------
area.min   <- 1000
RT.flex    <- 0.4 # This will need to go into standards for predictive RTs. 
blk.thresh <- 0.3
SN.min     <- 4


################################################################################

# Retention Time Table ----------------------------------------------------
RT.table <- combined %>%
  filter(Run.Type == "std") %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(RT.min = min(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.max = max(RT.Value, na.rm = TRUE)) %>%
  select(Metabolite.Name, RT.Value, RT.Expected, RT.min, RT.max, Run.Type) %>%
  unique()

# Blank Table ------------------------------------------------------------
blank.table <- combined %>%
  filter(Run.Type == "blk") %>%
  mutate(Blk.Area = Area.Value) %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(Blk.min = min(Area.Value)) %>%
  mutate(Blk.max = max(Area.Value)) %>%
  select(Metabolite.Name:Blk.max) %>%
  select(-Blk.Area) %>%
  unique()







# Create datasets for different flag types ------------------------------
SN.Area.Flags <- combined %>%
  arrange(Metabolite.Name) %>%
  mutate(SN.Flag       = ifelse(((SN.Value) < SN.min), "SN.Flag", NA)) %>%
  mutate(Area.Min.Flag = ifelse((Area.Value < area.min), "Area.Min.Flag", NA))

# Joining datasets---------------------------------------
add.RT.Flag <- SN.Area.Flags %>%
  group_by(Metabolite.Name) %>%
  left_join(RT.table, by = c("Metabolite.Name", "RT.Value", "RT.Expected", "Run.Type")) %>%
  mutate(RT.Flag = ifelse((RT.Value >= (RT.max + RT.flex) | RT.Value <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(-c("RT.max", "RT.min"))

add.blk.Flag <- add.RT.Flag %>%
  left_join(blank.table, by = c("Metabolite.Name", "Run.Type", "RT.Expected")) %>%
  mutate(Blank.Flag = ifelse((Area.Value / Blk.max) < blk.thresh, "Blank.Flag", NA)) %>%
  select(Replicate.Name, Metabolite.Name, Area.Value:RT.Value, RT.Expected, SN.Value, Run.Type, SN.Flag:Blank.Flag) %>%
  select(-c("Blk.min", "Blk.max"))


# Finally, combine all the flags. -----------------------------------------
final.table <- add.blk.Flag %>%
  mutate(all.Flags      = paste(SN.Flag, Area.Min.Flag, RT.Flag, Blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      = as.character(all.Flags %>% str_remove_all("NA, ") %>% str_remove_all("NA"))) %>%
  mutate(all.Flags      = ifelse(all.Flags == "", NA, all.Flags)) %>%
  mutate(Area.with.QC   = ifelse(is.na(Area.Min.Flag), Area.Value, NA)) %>%
  select(Replicate.Name:Area.Value, Area.with.QC, MZ.Value, SN.Value, RT.Value, RT.Expected, everything()) %>%
  ungroup(Metabolite.Name) %>%
  mutate(Metabolite.Name = as.character(Metabolite.Name)) 

RT.comparisons <- final.table %>%
  select(Replicate.Name, Metabolite.Name, RT.Value, RT.Expected) %>%
  mutate(RT.Difference = (abs(RT.Value - RT.Expected)))

# Print to file with comments and new name! -----------------------------
Description <- c("Hello! Welcome to the world of MSDIAL QE Quality Control! ",
                 "Minimum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: " ,
                 "Processed on: ")
Value <- c(NA, area.min, RT.flex, blk.thresh, SN.min, Sys.time())

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)


rm(list = ls()[!ls() %in% c("final.table", lsf.str())])

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/QC_Output_", currentDate, ".csv", sep = "")

write.csv(final.table, csvFileName, row.names = FALSE)


