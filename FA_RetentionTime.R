# Quality control script
source("FA_Functions.R")

# Use those FAs that have _Std in the Metabolite.Name. Compare to RT expected, not value.  

pattern = "combined"

# Import QC'd files ----------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

combined <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) 

# Create file isolating runs that have comments.
file.comments <- combined %>%
  select(Replicate.Name, Metabolite.Name, Comment) %>%
  filter(!nchar(Comment) < 5)

# Import expected retention time files
FA.expected <- read.csv("data_extras/FA_Expected_RT.csv", stringsAsFactors = FALSE) %>%
  rename(Metabolite.Name = Name) %>%
  rename(RT.Expected = RT) %>%
  rename(MZ.Value = m.z) %>%
  rename(Adduct.Type = Charge)
  
combined.expected <- combined %>%
  left_join(FA.expected %>% select(Metabolite.Name, RT.Expected)) %>%
  select(Replicate.Name, Metabolite.Name, Area.Value:SN.Value, RT.Expected, Reference.RT, Run.Type) %>%
  mutate(Cmpd.with.Std = ifelse(str_detect(Metabolite.Name, "_Std"), "Standard.Compound", 
                                ifelse(str_detect(Metabolite.Name, "IS"), "Internal.Standard", "NonStandard.Compound")))


# Separate replicates by size fractionation ---------------------------------------------------
# Currently we are just doing the 0.2 size fraction.
size.fraction_0.2 <- combined.expected %>%
  filter(!str_detect(Replicate.Name, "0.3")) 

size.fraction_0.3 <- combined.expected %>%
  filter(!str_detect(Replicate.Name, "0.2")) 

dummy.data <- size.fraction_0.2 %>%
  select(Replicate.Name, Metabolite.Name, RT.Value, RT.Expected, Cmpd.with.Std) %>%
  mutate(Random = RT.Value + sample(1, size = nrow(dummy.data), replace = TRUE))


# First plot of retention time ----------------------------------------------------
all.RT.plot <- ggplot(size.fraction_0.2, aes(x = Metabolite.Name, y = RT.Value, fill = Cmpd.with.Std)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
ggtitle("Fatty Acids: 0.2 Size Fraction")
print(all.RT.plot)

# And with dummy data
dummy.plot <- ggplot(dummy.data, aes(x = Metabolite.Name, y = Random, fill = Cmpd.with.Std)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Random RT Values")
print(dummy.plot)

require(gridExtra)
grid.arrange(all.RT.plot, dummy.plot, nrow=2)

################################################################################
# Retention Time Table ----------------------------------------------------
RT.Table <- size.fraction_0.2 %>%
  filter(!str_detect(Replicate.Name, "IS")) %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE),
         Min.RT.Value = min(RT.Value, na.rm = TRUE),
         Max.RT.Value = max(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.Diff = RT.Value - RT.Expected) %>% 
  select(Replicate.Name, Metabolite.Name, RT.Expected, RT.Value, Mean.RT.Value:RT.Diff)

std.RT.plot <- ggplot(RT.Table, aes(x = Metabolite.Name, y = RT.Value)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Fatty Acids: Standard Compounds")
print(std.RT.plot)

################################################################################
# FA 18 linear model test ----------------------------------------------------
FA.18 <- RT.Table %>%
  filter(str_detect(Metabolite.Name, "18")) %>%
  ungroup()

FA.18 <- FindStdDev(FA.18) 
FA.18 <- FA.18 %>% 
  mutate(Metabolite.Name = recode(Metabolite.Name,
                                  "FA 18:0_Std" = "0",
                                  "FA 18:1 cis-6_Std"= "1",
                                  "FA 18:2_Std" = "2",
                                  "FA 18:3_Std" = "3")) 
FA.18$Metabolite.Name <- as.numeric(as.character(FA.18$Metabolite.Name))

# Separate Replicates
ggplot(FA.18, aes(fill=Replicate.Name, y=RT.Value, x=Metabolite.Name)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=RT.Value - Std.dev, ymax=RT.Value + Std.dev), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Fatty Acid 18: 0-4")

# Grouped Replicates with Regression
FA.18.plot <- ggplot(FA.18, aes(x = Metabolite.Name, y = Mean.RT.Value)) +
  geom_point(aes(color = "Averaged RT Value")) +
  geom_point(aes(x = Metabolite.Name, y = RT.Expected, color = "Expected RT Value")) +
  geom_smooth(method = "lm", fullrange = TRUE, color = "black", formula = y ~ x) +
  ylim(0, 15) +
  xlim(0, 5) +
  theme(strip.text = element_text(size = 10)) +
  ggtitle("Fatty Acid 18: 0-4")
print(FA.18.plot)

FA18.model <- lm(formula = Mean.RT.Value ~ Metabolite.Name, data = FA.18)
FA18.fitted <- predict(FA18.model)
FA18.resid <- residuals(FA18.model)

myY <- (FA18.model$coefficients[2]*4) + FA18.model$coefficients[1]
realY <- size.fraction_0.2 %>%
  filter(Metabolite.Name == "FA 18:4") %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(!str_detect(Replicate.Name, "IS")) %>%
  summarize(mean(RT.Value))
  


################################################################################
# FA 20 test ----------------------------------------------------
FA.20 <- size.fraction_0.2 %>%
  filter(str_detect(Metabolite.Name, "20")) %>%
  filter(str_detect(Replicate.Name, "_Std_")) %>%
  filter(!str_detect(Metabolite.Name, "IS")) %>%
  group_by(Metabolite.Name) %>%
  mutate(Mean.RT.Value = mean(RT.Value, na.rm = TRUE))

FA.20 <- FindStdDev(FA.20) 

FA.20 <- FA.20 %>% # 4 + 5 are standards
  mutate(Metabolite.Name = recode(Metabolite.Name,
                                  "FA 20:0" = "0",
                                  "FA 20:1" = "1",
                                  "FA 20:2" = "2",
                                  "FA 20:3" = "3",
                                  "FA 20:4_Std" = "4",
                                  "FA 20:5_Std" = "5")) 
FA.20$Metabolite.Name <- as.numeric(as.character(FA.20$Metabolite.Name))

# Replicates included 
ggplot(FA.20, aes(fill=Replicate.Name, y=RT.Value, x=Metabolite.Name)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=RT.Value - Std.dev, ymax=RT.Value + Std.dev), width=.2,
                position=position_dodge(.9)) +
  ggtitle("Fatty Acid 20: 0-5")

# Replicates grouped
ggplot(FA.20, aes(x = Metabolite.Name, y = Mean.RT.Value)) +
  geom_point(aes(color = "Averaged RT Value")) +
  geom_point(aes(x = Metabolite.Name, y = RT.Expected, color = "Expected RT Value")) +
  geom_smooth(method = "lm", color = "black", formula = y ~ x) +
  ylim(0, 15) +
  xlim(0, 5) +
  theme(strip.text = element_text(size = 10)) +
  ggtitle("Fatty Acid 20: 0-4")

FA20.model <- lm(formula = Mean.RT.Value ~ Metabolite.Name, data = FA.20)
FA20.fitted <- predict(FA20.model)
FA20.resid <- residuals(FA20.model)
################################################################################
# K means clustering ----------------------------------------------------
ggplot(RT.Table, aes(RT.Value, Replicate.Name, color = Metabolite.Name)) + 
  geom_point() +
  ggtitle("Standard Retention Times")

ggplot() +
  geom_point(data=RT.Table, aes(x=RT.Value, y=Replicate.Name, color = Metabolite.Name, size = 1)) + 
  geom_point(data=RT.Table, aes(x=RT.Expected, y=Replicate.Name)) +
  ggtitle("Expected vs Real Retention Times")

cluster.test <- RT.Table %>%
  arrange(Metabolite.Name)

set.seed(20)
RTCluster <- kmeans(cluster.test[, 8], 2, nstart = 20)
RTCluster

RTCluster$cluster <- as.factor(RTCluster$cluster)
ggplot(cluster.test, aes(RT.Diff, Replicate.Name, color = RTCluster$cluster)) + 
  geom_point() +
  ggtitle("K-means clustering: RT Value Differences")

cluster.test$cluster <- RTCluster$cluster

# RT differences plot
RT.Table.clustered <- cluster.test %>%
  group_by(Metabolite.Name) %>%
  mutate(High.Low = as.character(ifelse(cluster == 1, "Low", "High"))) %>%
  select(-cluster) %>%
  # TESTING AREA #
  unique() 

RT.Plot <- ggplot(RT.Table.clustered, aes(x = Replicate.Name, y = RT.Diff, fill = High.Low)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Metabolite.Name, scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        axis.text.y = element_text(size = 5),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Retention Time Differences")
print(RT.Plot)

################################################################################
## Constructing Tolerances

Tolerance.Table.High <- RT.Table.clustered %>%
  filter(High.Low == "High") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.High = mean(RT.Value)) %>%
  mutate(Ave.High.Diff = abs(RT.Expected - Ave.High)) %>%
  select(-RT.Value) %>%
  unique()

Tolerance.Table.Low <- RT.Table.clustered %>%
  filter(High.Low == "Low") %>%
  unique() %>%
  ## TESTING
  select(Metabolite.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>%
  group_by(Metabolite.Name) %>%
  mutate(Ave.Low= mean(RT.Value)) %>%
  mutate(Ave.Low.Diff = abs(RT.Expected - Ave.Low)) %>%
  select(-RT.Value) %>%
  unique()

# Tolerance production ----------------------------------------------------
## in progress
Full.Tolerance.Table <- size.fraction_0.2 %>%
  select(Metabolite.Name, Replicate.Name, RT.Expected, RT.Value) %>%
  #filter(Metabolite.Name == "FA 16:0_Std") %>% 
  filter(str_detect(Metabolite.Name, "_Std")) %>%
  # TESTING AREA #
  left_join(Tolerance.Table.High) %>%
  left_join(Tolerance.Table.Low) %>%
  select(Metabolite.Name, RT.Expected, Ave.High:Ave.Low.Diff) %>%
  unique()

FA16_HighTolerance = unique(Full.Tolerance.Table$Ave.High.Diff)
FA16_LowTolerance = unique(Full.Tolerance.Table$Ave.Low.Diff)



  

