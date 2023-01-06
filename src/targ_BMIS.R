# B-MIS 

# Match QC'd data with Internal Standards list ----------------------------------------------------------
data_withIS <- QCd.data %>%
  filter(Precursor.Ion.Name %in% Internal.Standards$Compound.Name)

data_noIS <- QCd.data %>%
  filter(!Precursor.Ion.Name %in% Internal.Standards$Compound.Name)

# Create Internal Standard data -----------------------------------------------------------------------
IS_data <- data_withIS %>%
  select(Replicate.Name, Precursor.Ion.Name, Area.with.QC) %>%
  mutate(Mass.Feature = Precursor.Ion.Name) %>%
  select(-Precursor.Ion.Name)

# Add injection volume from sampkey here. See original QC.
# Rbind sampkey to create IS data to identify problematic compounds/replicates. See original QC.
 
# Identify internal standards without an Area, i.e. any NA values.
IS.Issues <- IS_data[is.na(IS_data$Area.with.QC), ]
write.csv(IS.Issues, paste("mesocenter_pos_output/MSDial_InternalStdIssues_", currentDate, ".csv", sep = ""))

# Visualize raw areas of Internal Standards -------------------------------------------------------------
IS_raw_area_plot <- ggplot(IS_data, aes(x = Replicate.Name, y = Area.with.QC)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~Mass.Feature, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10)) +
  ggtitle("Internal Standard Raw Areas")

plotFileName <- paste("mesocenter_pos_output/IS.Raw.Areas_", currentDate, ".png", sep = "")

ggsave(file = plotFileName, dpi = 600, width = 8, height = 6, units = "in")
print(IS_raw_area_plot)

# Edit data so names match, test that names are equal across sample sets---------------------------------
data_long <- data_noIS %>%
  mutate(Mass.Feature = Precursor.Ion.Name) %>%
  select(Replicate.Name, Mass.Feature, Area.with.QC) %>%
  arrange(Replicate.Name)

## Test here for replicate names are equal across sets. See original QC.

# Caluclate mean values for each Internal Standard--------------------------------------------------------
IS_means <- IS_data %>%
  mutate(Mass.Feature = as.factor(Mass.Feature)) %>%
  group_by(Mass.Feature) %>%
  summarise(Average.Area = mean(as.numeric(Area.with.QC), na.rm = TRUE)) %>%
  mutate(Mass.Feature = as.character(Mass.Feature))

# Normalize to each internal Standard--------------------------------------------------------------------
data_bound <- rbind(IS_data, data_long) %>%
  arrange(Mass.Feature)

Split_Dat <- list()

for (i in 1:length(unique(IS_data$Mass.Feature))) {
  Split_Dat[[i]] <- data_bound %>%
    mutate(MIS = unique(IS_data$Mass.Feature)[i]) %>%
    left_join(IS_data %>%
                dplyr::rename(MIS = Mass.Feature, IS_Area = Area.with.QC) %>%
                select(MIS, Replicate.Name, IS_Area), by = c("Replicate.Name", "MIS")) %>%
    left_join(IS_means %>%
                dplyr::rename(MIS = Mass.Feature), by = "MIS") %>%
    mutate(Adjusted.Area = Area.with.QC/IS_Area*Average.Area)
}

data_area_normed <- do.call(rbind, Split_Dat) %>%
  select(-IS_Area, -Average.Area)

# Standardize name structure to: Date_type_ID_replicate_anythingextra -----------------------------------
data_standardized <- data_area_normed %>%
  separate(Replicate.Name, c("runDate", "type", "SampID", "replicate"), "_") %>%
  mutate(Run.Cmpd = paste(data_area_normed$Replicate.Name, data_area_normed$Mass.Feature))

# Find the B-MIS for each MassFeature-------------------------------------------------------------------

# Look only at the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo),
# then choose which IS reduces the RSD the most (Poo.Picked.IS)
Poodata <- data_standardized %>%
  filter(type == "Poo") %>%
  group_by(SampID, Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted.Area, na.rm = TRUE) / mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(Mass.Feature, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE)) %>%
  mutate(RSD_ofPoo = ifelse(RSD_ofPoo == "NaN", NA, RSD_ofPoo))

Poodata2 <- Poodata %>%
  left_join(Poodata %>% group_by(Mass.Feature) %>%
              summarise(Poo.Picked.IS = unique(MIS)[which.min(RSD_ofPoo)] [1]))


# Get the original RSD, calculate RSD change, decide if MIS is acceptable -------------------------------
Poodata3 <- left_join(Poodata, Poodata2 %>%
                       #filter(MIS == "Inj_vol" ) %>% ## Filtered just for this dumb run
                       mutate(Orig_RSD = RSD_ofPoo) %>%
                       select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percent.Change = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percent.Change > cut.off & Orig_RSD > cut.off2))

# Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----------------------------------------

# Adds a column that has the BMIS, not just Poo.Picked.IS
# Changes the FinalBMIS to inject_volume if it's no good
Fixed.poodata <- Poodata3 %>%
  filter(MIS == Poo.Picked.IS) %>%
  mutate(FinalBMIS = ifelse(accept_MIS == "FALSE", "Inj_vol", Poo.Picked.IS)) %>%
  mutate(FinalRSD = RSD_ofPoo)

New.poodata <- Poodata %>%
  left_join(Fixed.poodata) %>% #select(Mass.Feature, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)

Try <- New.poodata %>%
  filter(FinalBMIS != "Inj_vol")

QuickReport <- print(paste("Percent of Mass Features that picked a BMIS:",
                           length(Try$Mass.Feature) / length(New.poodata$Mass.Feature), "|",
                           "RSD improvement cutoff", cut.off, "|",
                           "RSD minimum cutoff", cut.off2,
                           sep = " "))

# Save some quick reports here. See original QC.

# Evaluate and visualize the results of your BMIS cutoff. See original QC.


# Return data that is normalized via BMIS----------------------------------------------------------------
BMIS_normalized <- New.poodata %>% select(Mass.Feature, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(data_standardized, by = "Mass.Feature") %>%
  filter(MIS == FinalBMIS)

currentDate <- Sys.Date()
csvFileName <- paste("mesocenter_pos_output/MSDial_BMIS_Output_", column, "_", currentDate, ".csv", sep = "")


write.csv(BMIS_normalized, csvFileName, row.names = FALSE)
