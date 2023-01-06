# This code retrieves mol/L from peak areas of targeted compounds.

# Get response factors for transect compounds ----------------------------------
Full_data_RF <- Full_data %>%
  mutate(RF = Area.with.QC/Conc..uM) %>%
  filter(!Compound.Type == "Internal Standard") %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2))

# In HILIC compounds, filter mixes. See original QC.

# Calculate RF max and min using only standards in water.
Full_data_RF_dimensions <- Full_data_RF %>%
  filter(Type == "Standards_Water") %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(RF.max = max(RF),
         RF.min = min(RF))

Full_data_RF_dimensions$RF.max[is.infinite(Full_data_RF_dimensions$RF.max) | is.nan(Full_data_RF_dimensions$RF.max) ] <- NA
Full_data_RF_dimensions$RF.min[is.infinite(Full_data_RF_dimensions$RF.min) | is.nan(Full_data_RF_dimensions$RF.min) ] <- NA

Full_data_RF_dimensions <- Full_data_RF_dimensions %>%
  mutate(RF.diff = RF.max/RF.min) %>%
  unique()

# Calculate response factor ratios ----------------------------------------
# Calculate the response factor ratios using (Standards in Matrix - Water in Matrix) / (Standards in water) for each replicate.
temp_RF_ratios <- Full_data_RF %>%
  group_by(Precursor.Ion.Name, Type) %>%
  mutate(RF.mean.per_sampleID = mean(RF, na.rm = TRUE)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Type, RF.mean.per_sampleID) %>%
  unique()

print(paste("NAs or NaNs in the calculated response factor ratios:", TRUE %in% is.na(temp_RF_ratios)))
metabolite.issues <- temp_RF_ratios[is.nan(temp_RF_ratios$RF.mean.per_sampleID),]
print(unique(metabolite.issues$Precursor.Ion.Name))

## Here we calculate Full_data_RF_ratios, see original QC.
# Missing correct standard runtypes so skipping this step.
Full_data_RF_ratios <- temp_RF_ratios %>%
  mutate(RF.ratio = 1) %>%
  select(Precursor.Ion.Name, RF.ratio) %>%
  unique()

# If applicable, supplement data with information from calculated Ingalls QE.RF ratios. See original QC.
currentDate = Sys.Date()
write.csv(Full_data_RF_ratios, paste("mesocenter_pos_output/MSDIAL_ResponseFactorRatios_", currentDate, ".csv", sep = ""))

# Quantify samples for the BMIS'd dataset ---------------------------------
BMISd_data_filtered <- BMISd_data %>%
  separate(Run.Cmpd, c("Sample.Name"), extra = "drop", fill = "right") %>%
  mutate(Precursor.Ion.Name = Mass.Feature) %>%
  filter(Precursor.Ion.Name %in% Full_data_RF_ratios$Precursor.Ion.Name) %>%
  left_join(Full_data_RF_ratios) %>%
  left_join(Full_data_RF_dimensions %>% select(Precursor.Ion.Name, RF.max, RF.min) %>% 
              unique(), by = "Precursor.Ion.Name") %>%
  select(Precursor.Ion.Name, FinalBMIS, Sample.Name, Adjusted.Area, everything())

# Calculate umol/vial for compounds without an internal standard ----------
Quantitative_data <- BMISd_data_filtered %>%
  mutate(RF.ave = as.numeric(rowMeans(BMISd_data_filtered[, c("RF.min", "RF.max")]))) %>%
  mutate(umol.in.vial.ave = Adjusted.Area/RF.ave/RF.ratio,
         umol.in.vial.max = Adjusted.Area/RF.min/RF.ratio,
         umol.in.vial.min = Adjusted.Area/RF.max/RF.ratio) %>%
  select(Precursor.Ion.Name:Adjusted.Area, everything())

# Pull out data for matched internal standards ----------------------------
IS_key <- BMISd_data_filtered %>%
  select(FinalBMIS, Precursor.Ion.Name) %>%
  unique() %>%
  left_join(OG_IS_key %>% select(FinalBMIS, Concentration_nM))# %>%
 # filter(str_detect(FinalBMIS, Precursor.Ion.Name))

# Calculate umol/vial for compounds with matched internal standards -----------------
IS_data <- Full_data %>%
  filter(Precursor.Ion.Name %in% IS_key$FinalBMIS) %>%
  mutate(IS_Area = Area.with.QC,
         FinalBMIS = Precursor.Ion.Name) %>%
  select(IS_Area, FinalBMIS, Replicate.Name) %>%
  left_join(IS_key %>% select(FinalBMIS, Precursor.Ion.Name, Concentration_nM))

matched.IS.compounds <- data.frame(Compounds = c(IS_key[ ,"FinalBMIS"], as.character(IS_key[ ,"Precursor.Ion.Name"])))

IS_sample_data <- QCd_data %>%
  left_join(IS_data %>% select(FinalBMIS, Precursor.Ion.Name, Concentration_nM), by = "Precursor.Ion.Name") %>%
  unique() %>%
  filter(Precursor.Ion.Name %in% matched.IS.compounds$Compounds) %>%
  filter(!str_detect(Replicate.Name, "Std")) %>%
  mutate(Std.Type = ifelse(str_detect(Precursor.Ion.Name, ","), "Internal_std", "Standard")) %>%
  mutate(testcol1 = ifelse(str_detect(Precursor.Ion.Name, ","), 
                           sapply(strsplit(Precursor.Ion.Name, ","), `[`, 1), Precursor.Ion.Name)) %>%
  mutate(Names = ifelse(str_detect(testcol1, "-"), sapply(strsplit(testcol1, "-"), `[`, 2), testcol1)) %>%
  mutate(Pairs = ifelse(!str_detect(Precursor.Ion.Name, ","), 
                        Precursor.Ion.Name, paste(Names, "IS", sep = "_"))) %>%
  select(-c("Pairs", "testcol1")) %>%
  arrange(Replicate.Name) %>%
  group_by(Names) %>%
  group_split()

IS.mid_frame <- lapply(IS_sample_data, function(x) group_by(x, Replicate.Name))

IS.mid_frame2 <- lapply(IS.mid_frame,
                        function(x) mutate(x,
                        #umol.in.vial_IS = (Area.with.QC[Std.Type == "Standard"] / Area.with.QC[Std.Type == "Internal_std"]) * (Concentration_nM[Std.Type == "Standard"]/1000)))
                        umol.in.vial_IS = 100))

IS_sample_data <- do.call(rbind, IS.mid_frame2) %>%
  filter(!str_detect(Precursor.Ion.Name, ",")) %>%
  dplyr::rename(Sample.Name = Replicate.Name) %>%
  select(Sample.Name:Area.with.QC, Concentration_nM, umol.in.vial_IS)

rm(list = c("matched.IS.compounds", "QC_data", "IS.mid_frame", "IS.mid_frame2"))


# Add matched IS_smp info back into main frame ------------------------------------------------
All.Info <- Quantitative_data %>%
  select(Precursor.Ion.Name, runDate:replicate, Adjusted.Area, Area.with.QC, RF.ratio:umol.in.vial.min) %>%
  unite(Sample.Name, c("runDate", "type", "SampID", "replicate"), remove = FALSE) %>%
  left_join(IS_sample_data %>% select(Sample.Name, Precursor.Ion.Name, umol.in.vial_IS)) %>%
  mutate(umol.in.vial.ave = ifelse(is.na(umol.in.vial_IS), umol.in.vial.ave, umol.in.vial_IS),
         umol.in.vial.max = ifelse(is.na(umol.in.vial_IS), umol.in.vial.max, NA),
         umol.in.vial.min = ifelse(is.na(umol.in.vial_IS), umol.in.vial.min, NA)) %>%
  dplyr::rename(Replicate.Name = Sample.Name) %>%
  select(-runDate, -type, -replicate) # %>%
  #filter(!str_detect(Replicate.Name, "DDA"))

# Add in dilution factor and filtered volume --------------------------------------------------
All.Info.Quantitative <- All.Info %>%
  mutate(nmol.in.Enviro.ave = (umol.in.vial.ave*10^-6*Reconstitution.Volume/Volume.Filtered*1000*Dilution.Factor)) %>%
  left_join(Full_data %>% select(Precursor.Ion.Name, Emperical.Formula)) %>%
  unique()

# Get molecules of carbon and nitrogen ------------------------------------
All.Info.Molecules <- All.Info.Quantitative  %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"),
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1,
                    str_extract(Emperical.Formula, "N\\d"))) %>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(nmol.C.ave = nmol.in.Enviro.ave*C,
         nmol.N.ave = nmol.in.Enviro.ave*N ) %>%
  select(Precursor.Ion.Name, SampID, Replicate.Name, everything())

# Summarize for each metabolite ------------------------------------
All.Info.Summed <- All.Info.Molecules %>%
  group_by(Precursor.Ion.Name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T)) %>%
  arrange(desc(nmol.Enviro.med))

# Summarize total carbon and nitrogen for each compound ------------------------------------
Final.All.perSampID <- All.Info.Molecules %>%
  select(SampID, nmol.C.ave, nmol.N.ave) %>%
  group_by(SampID) %>%
  summarise(totalCmeasured_nM_perID = sum(as.numeric(nmol.C.ave), na.rm = TRUE),
            totalNmeasured_nM_perID = sum(as.numeric(nmol.N.ave), na.rm = TRUE))


# Calculate mole fractions of each compound ------------------------------------
Final.Quantitative <- All.Info.Molecules %>%
  unique() %>%
  left_join(Final.All.perSampID) %>%
  mutate(ratioCN = totalCmeasured_nM_perID / totalNmeasured_nM_perID) %>%
  mutate(molFractionC = nmol.C.ave/totalCmeasured_nM_perID,
         molFractionN = nmol.N.ave/totalNmeasured_nM_perID) %>%
  select(Precursor.Ion.Name, Replicate.Name, Adjusted.Area, Area.with.QC, RF.ratio:molFractionN) %>%
  unique()

Final.Quantitative.Summed <- Final.Quantitative %>%
  group_by(Precursor.Ion.Name) %>%
  summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
            nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
            nmol.C.med = median(nmol.C.ave, na.rm  = T),
            nmol.C.min = min(nmol.C.ave, na.rm  = T),
            nmol.C.max = max(nmol.C.ave, na.rm  = T),
            mol.C.Fraction.med = median(molFractionC, na.rm = T),
            mol.C.Fraction.min = min(molFractionC, na.rm = T),
            mol.C.Fraction.max = max(molFractionC, na.rm = T)) %>%
  arrange(desc(Precursor.Ion.Name))
