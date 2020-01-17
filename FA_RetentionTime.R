# Quality control script
source("FA_Functions.R")

# Use those FAs that have _Std in the Metabolite.Name. Compare to RT expected, not value.  

pattern = "combined"

# Import QC'd files ----------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

fatty.acids <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

# Create file isolating runs that have comments.
file.comments <- fatty.acids %>%
  select(Replicate.Name, Metabolite.Name, Comment) %>%
  filter(!nchar(Comment) < 5)

# Import expected retention time files ----------------------------------------------------
FA.expected <- read.csv("data_extras/FA_Expected_RT.csv", stringsAsFactors = FALSE) %>%
  rename(Metabolite.Name = Name) %>%
  rename(RT.Expected = RT) %>%
  select(Metabolite.Name, RT.Expected)
  
FA.expected <- fatty.acids %>%
  left_join(FA.expected) %>%
  mutate(Cmpd.with.Std = ifelse(str_detect(Metabolite.Name, "_Std"), "Standard.Compound", 
                                ifelse(str_detect(Metabolite.Name, "IS"), "Internal.Standard", "NonStandard.Compound"))) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE),
         Min.RT.Value = min(RT.Value, na.rm = TRUE),
         Max.RT.Value = max(RT.Value, na.rm = TRUE)) %>%
  select(Replicate.Name, Metabolite.Name, RT.Value, RT.Expected, Mean.RT.Value:Max.RT.Value, Cmpd.with.Std) 

# Separate replicates by size fractionation ---------------------------------------------------
# Currently we are just doing the 0.2 size fraction.
size.fraction_0.2 <- FA.expected %>%
  filter(!str_detect(Replicate.Name, "0.3")) 

size.fraction_0.3 <- FA.expected %>%
  filter(!str_detect(Replicate.Name, "0.2")) 

# Create random values
dummy.data <- size.fraction_0.2 %>%
  group_by(Metabolite.Name, Replicate.Name) %>%
  mutate(Random.RT = runif(1, RT.Value, RT.Value + 1)) %>%
  group_by(Metabolite.Name) %>%
  mutate(Random.Mean = mean(Random.RT)) %>%
  select(Replicate.Name:RT.Expected, Mean.RT.Value, Random.RT, Random.Mean, Cmpd.with.Std)


# First plot of retention times ----------------------------------------------------
real.RT.plot <- ggplot(size.fraction_0.2, aes(x = Metabolite.Name, y = Mean.RT.Value, fill = Cmpd.with.Std)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Mean.RT.Value, digits = 2)), 
            position=position_dodge(width=0.9), 
            vjust=-0.25, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
ggtitle("Fatty Acids: 0.2 Size Fraction")

# And with dummy data
dummy.plot <- ggplot(dummy.data, aes(x = Metabolite.Name, y = Random.Mean, fill = Cmpd.with.Std)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Random.Mean, digits = 2)), 
            position=position_dodge(width=0.9), 
            vjust=-0.25, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Random RT Values")

require(gridExtra)
grid.arrange(real.RT.plot, dummy.plot, nrow=2)

################################################################################
# Retention Time Tables ----------------------------------------------------
real.standards <- size.fraction_0.2 %>%
  filter(!str_detect(Replicate.Name, "IS")) %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  mutate(RT.Diff = RT.Value - RT.Expected) %>% 
  select(Replicate.Name, Metabolite.Name, RT.Expected, RT.Value, Mean.RT.Value:RT.Diff)

dummy.standards <- dummy.data %>%
  filter(!str_detect(Replicate.Name, "IS")) %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  group_by(Metabolite.Name) %>%
  mutate(Random.RT.Mean = mean(Random.RT)) 


# Real + Dummy standards plotted ----------------------------------------------------
std.RT.plot <- ggplot(real.standards, aes(x = Metabolite.Name, y = Mean.RT.Value, label = Mean.RT.Value)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Mean.RT.Value, digits = 2)), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Fatty Acids: Standard Compounds")

dummy.std.RT.plot <- ggplot(dummy.standards, aes(x = Metabolite.Name, y = Random.Mean)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Random.Mean, digits = 2)), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Random Standards")

stds.together <- dummy.standards %>%
  select(Metabolite.Name, Mean.RT.Value, Random.RT.Mean) %>%
  rename(Real.RT.Value = Mean.RT.Value) %>%
  rename(Random.RT.Value = Random.RT.Mean) %>%
  unique()

stds.together <- melt(stds.together) %>%
  rename(RT.Type = variable) %>%
  rename(Retention.Time = value)
stds.together.plot <- ggplot(stds.together, aes(Metabolite.Name, Retention.Time, fill = RT.Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Real vs Random Standards Retention Times")
stds.together.plot <- ggplotly(stds.together.plot)
stds.together.plot

# Real + Dummy standards plotted ----------------------------------------------------
stds.differences <- stds.together %>%
  group_by(Metabolite.Name) %>%
  mutate(Difference = abs(Retention.Time[RT.Type == "Real.RT.Value"] - Retention.Time[RT.Type == "Random.RT.Value"])) %>%
  select(Metabolite.Name, Difference) %>%
  unique() %>%
  ungroup() %>%
  mutate(Min.diff = min(Difference)) %>%
  mutate(Max.diff = max(Difference)) %>%
  select(Min.diff, Max.diff) %>%
  unique()

high.estimate <- stds.differences$Max.diff
low.estimate <- stds.differences$Min.diff

# High and Low estimate ----------------------------------------------------
high.low.predictions <- size.fraction_0.2 %>%
  select(Replicate.Name, Metabolite.Name, RT.Value) %>%
  mutate(High.RT = RT.Value + high.estimate) %>%
  mutate(Low.RT = RT.Value + low.estimate) %>%
  filter(!str_detect(Metabolite.Name, "Std|IS"))

high.low.plot.table <- melt(high.low) %>%
  rename(RT.Value.Prediction = variable) %>%
  rename(RT.Value = value)


prediction.plot <- ggplot(high.low.plot.table, aes(fill = RT.Value.Prediction, y = RT.Value, x = Metabolite.Name)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("firebrick3", "navyblue", "skyblue")) +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) + 
  ggtitle("Predicted Fatty Acid Retention Times")
prediction.plot <- ggplotly(prediction.plot)
prediction.plot
