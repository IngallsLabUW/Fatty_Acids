# Quality control script
source("FA_Functions.R")

# Use those FAs that have _Std in the Metabolite.Name. compare to RT expected, not value. know that 1.7 difference between expected and value is pretty significant in the RP ? space. 
# Check out weird blank flags. 
# do different RT predictions/QCs for size fractionation. 

pattern = "combined"

# Import QC'd files ----------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

combined <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

# Create file isolating runs that have comments.
file.comments <- combined %>%
  select(Replicate.Name, Comment) %>%
  filter(!nchar(Comment) < 5)

# Import expected retention time files
FA.expected <- read.csv("data_extras/FA_Expected_RT.csv", stringsAsFactors = FALSE) %>%
  rename(Metabolite.Name = Name) %>%
  rename(RT.Expected = RT) %>%
  rename(MZ.Value = m.z) %>%
  rename(Adduct.Type = Charge)
  
combined.expected <- combined %>%
  left_join(FA.expected %>% select(Metabolite.Name, RT.Expected)) %>%
  select(Replicate.Name, Metabolite.Name, Area.Value:SN.Value, RT.Expected, Reference.RT, Run.Type)


# Separate replicates by size fractionation ---------------------------------------------------
size.fraction_0.2 <- combined.expected %>%
  filter(!str_detect(Replicate.Name, "0.3")) %>%
  mutate(Cmpd.with.Std = ifelse(str_detect(Metabolite.Name, "_Std"), "Standard.Compound", 
                                ifelse(str_detect(Metabolite.Name, "IS"), "Internal.Standard","NonStandard.Compound")))

# First plot of retention time ----------------------------------------------------
first.plot <- ggplot(size.fraction_0.2, aes(x = Metabolite.Name, y = RT.Value, fill = Cmpd.with.Std)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
ggtitle("Fatty Acids: 0.2 Size Fraction")
print(first.plot)


################################################################################
# Retention Time Table ----------------------------------------------------
RT.Table <- size.fraction_0.2 %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE),
         Min.RT.Value = min(RT.Value, na.rm = TRUE),
         Max.RT.Value = max(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.Diff = RT.Value - RT.Expected) %>% 
  mutate(RT.Diff.abs = abs(RT.Value - RT.Expected)) %>%
  mutate(Min.RT.Diff = min(RT.Diff),
         Max.RT.Diff = max(RT.Diff)) %>%
  mutate(Midrange.RT.Diff = (min(RT.Diff) + max(RT.Diff)) / 2) %>%
  mutate(High.Low = ifelse(RT.Diff > Midrange.RT.Diff, "High", "Low")) %>%
  select(Replicate.Name, Metabolite.Name, RT.Expected, RT.Value, Mean.RT.Value:High.Low)
  

# RT differences plot
RT.Table.Plot <- RT.Table %>%
  group_by(Metabolite.Name) %>%
  # TESTING AREA #
  select(Metabolite.Name, Replicate.Name, RT.Diff.abs, Midrange.RT.Diff, High.Low) %>%
  unique() 

RT.Plot <- ggplot(RT.Table.Plot, aes(x = Replicate.Name, y = RT.Diff.abs, fill = High.Low)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Metabolite.Name, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 5),
        legend.position = "top",
        strip.text = element_text(size = 10))
  ggtitle("Retention Time Differences")
print(RT.Plot)


Tolerance.Table.High <- RT.Table %>%
  filter(High.Low == "High") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.High = mean(RT.Value)) %>%
  mutate(Ave.High.Diff = abs(RT.Expected - Ave.High)) %>%
  select(-RT.Value) %>%
  unique()
  
Tolerance.Table.Low <- RT.Table %>%
  filter(High.Low == "Low") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.Low= mean(RT.Value)) %>%
  mutate(Ave.Low.Diff = abs(RT.Expected - Ave.Low)) %>%
  select(-RT.Value) %>%
  unique()



# Tolerance production ----------------------------------------------------
Full.Tolerance.Table <- size.fraction_0.2 %>%
  select(Metabolite.Name, Replicate.Name, RT.Expected, RT.Value) %>%
  filter(Metabolite.Name == "FA 16:0_Std") %>% 
  filter(str_detect(Replicate.Name, "Smp|Blk")) %>%
  # TESTING AREA #
  left_join(Tolerance.Table.High) %>%
  left_join(Tolerance.Table.Low) %>%
  select(Metabolite.Name, RT.Expected, Ave.High:Ave.Low.Diff) %>%
  unique()

FA16_HighTolerance = unique(Full.Tolerance.Table$Ave.High.Diff)
FA16_LowTolerance = unique(Full.Tolerance.Table$Ave.Low.Diff)




################################################################################

  filter(str_detect(Replicate.Name, regex("std", ignore_case = TRUE))) %>%
  filter(!str_detect(Replicate.Name, "IS")) %>%
  mutate(Replicate.Name = gsub("(.*)_.*", "\\1", Replicate.Name)) %>%
  group_by(Metabolite.Name, Replicate.Name) %>%
  mutate(RT.Value.ave = mean(RT.Value)) %>%
  select(-RT.Value) %>%
  unique() %>%
  mutate(RT.Value = na_if(RT.Value, 0))


RT.table2 <- RT.table %>%
  mutate(Std.Type = ifelse(str_detect(Metabolite.Name, "IS"), "Internal_std", "Standard")) %>%
  mutate(Names = gsub("(.*)_.*", "\\1", Metabolite.Name)) %>%
  group_by(Names) %>%
  group_split()


RT.midframe <- lapply(RT.table2, function(x) group_by(x, Replicate.Name))

testing <- RT.table2[[1]]
testing <- testing %>% 
  select(Replicate.Name, Metabolite.Name, RT.Difference, Std.Type) %>%
  unique() %>%
  group_by(Replicate.Name, Metabolite.Name) %>%
  mutate(IS.Std.Diff = (abs(RT.Difference[Std.Type == "Internal_std"] - RT.Difference[Std.Type == "Standard"])))


RT.midframe2 <- lapply(RT.midframe, function(x) mutate(x, IS.Std.Diff = (abs(RT.Difference[Std.Type == "Standard"] - RT.Difference[Std.Type == "Internal_std"]))))



Wei.IS.smp.data.transect <- do.call(rbind, Wei.IS.mid_frame2) %>%
  filter(!str_detect(Compound.Name, ",")) %>%
  rename(Sample.Name = ReplicateName) %>%
  select(Sample.Name:Area.with.QC, Concentration_nM, umol.in.vial_IS)

  
  #filter(str_detect(Replicate.Name, regex("std", ignore_case = TRUE))) %>%

  

