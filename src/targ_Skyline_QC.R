# Skyline TQS + QE Quality Control

# Import datafiles and accompanying master files --------------------------------------------------------------
if (instrument_pattern == "TQS") {
  filenames <- RemoveCsv(list.files(path = "data_extras", pattern = "master", ignore.case = TRUE))
  filepath <- file.path("data_extras", paste(filenames, ".csv", sep = ""))
  master_file <- assign(make.names(filenames), read.csv(filepath, stringsAsFactors = FALSE)) %>%
    dplyr::rename(Second.Trace = X2nd.trace)
}
skyline_output <- combined_skyline

# Sanity check for runtypes  ---------------------------------------------------------------------
# Stop program if this run has more or fewer runtypes than the normal std, blk, poo, and smp.
skyline_runtypes <- IdentifyRunTypes(skyline_output)

# Filter out redundant standard mixes in HILIC runs here, returns skyline_output -----------------

# Depending on instrument_pattern, create comparison tables --------------------------------------

if (instrument_pattern == "TQS") {
  # Fragment check here, producing fragments_checked. See original QC.
  fragments_checked <- skyline_output
  
  # Ion Ratio check here, producing ion_ratio_table. See original QC.
  # Blank Table check here, producing blank_table. See original QC.
  # Height Table check here, producing height_table. See original QC.
  # Signal to Noise check here, producing SN_table. See original QC.

} else{
  
  print(paste("This is a", instrument_pattern, "run. No fragmentation check necessary."))
  # Blank Table and Heigh Table checks for QE runs here. See original QC.
}


# Retention Times 
# Find the minimum and maximum Retention Times and take the average.
# Use this as a reference table for acceptable Retention Times.
RT_table <- skyline_output %>%
  filter(str_detect(Replicate.Name, regex("Std", ignore_case = TRUE))) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(RT.min = min(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.max = max(Retention.Time, na.rm = TRUE)) %>%
  mutate(RT.Reference = mean(as.numeric(Retention.Time), na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, RT.min, RT.max) %>%
  unique()

# Area  
# Isolate all pooled and sample Areas.
area_table <- skyline_output %>%
  select(Replicate.Name, Precursor.Ion.Name, Area) %>%
  filter(str_detect(Replicate.Name, regex("Smp|Poo", ignore_case = TRUE)))

# Signal to Noise calculations here, returns SN_table. See original QC.
# May want to improve this step.

# Construct final comparative table ---------------------------------------

if (instrument_pattern == "TQS") {
  # CheckFragments of standards and samples, returns all_standards and all_samples. See original QC.
  # Ion Ratio Flags. See original QC.
  all_samples <- skyline_output
  } else {
  all_samples <- skyline_output
}

# Retention Time Flags, returns RT_flags_added. See original QC.
RT_flags_added <- all_samples 

# Blank Flags, returns Blank_flags_added. See original QC.
Blank_flags_added <- RT_flags_added 

# Height Flags, returns Height_flags_added. See original QC.
Height_flags_added <- Blank_flags_added

# Area Flags  ---------------------------------------
# If the Area is less than the area.min value, add a flag.
Area_flags_added <- Height_flags_added %>%
  mutate(area.min.Flag = ifelse((Area < area.min), "area.min.Flag", NA)) %>%
  mutate(Area.with.QC   = ifelse(is.na(area.min.Flag), Area, NA)) %>%
  select(Replicate.Name:Area, Area.with.QC, everything())

# Signal to Noise Flags, returns SN_flags_added. See original QC.
SN_flags_added <- Area_flags_added

# All Flags  ---------------------------------------
# Add a column with all flags from the previous steps. 
semifinal_table <- SN_flags_added %>%
  unite(all.Flags, contains("Flag"), sep = ", ", remove = FALSE) %>%
  mutate(all.Flags = as.character(all.Flags %>% str_remove_all("NA,|NA") %>% trimws()))
semifinal_table$all.Flags <- gsub('^\\,|\\,$', '', semifinal_table$all.Flags)

final_table <- semifinal_table %>%
  select(Replicate.Name:Column, contains("Flag"))
final_table[final_table==""]<-NA

# Remove Secondary trace from TQS here, returns final_table.

# Standards & blank addition. See original QC.

# Print to file with comments and a new name ------------------------------
if (instrument_pattern == "TQS") {
  Description <- c(as.character(anydate(Sys.Date())),
                   "Hello! Welcome to the world of Skyline TQS Quality Control! ",
                   "Maximum height for a real peak: ",
                   "Minimum height for a real peak: ",
                   "Maximum area for a real peak: ",
                   "RT flexibility: ",
                   "Blank can be this fraction of a sample: ",
                   "S/N ratio: " ,
                   "Ion ratio flexibility", 
                   "Processed on: ")
  
  Value <- c(NA, NA, height.max, height.min, area.min, RT.flex, blk.thresh, SN.min, IR.flex, Sys.time())
  
} # else do something different here.

df <- data.frame(Description, Value)
final_table <- bind_rows(df, final_table)

# Remove unnecessary values here if you like.