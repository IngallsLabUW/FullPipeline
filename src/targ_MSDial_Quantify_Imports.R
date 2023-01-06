# Imports for the quantification step


# Import standards and filter NAs ---------------------------------------------------------------
# filename <- RemoveCsv(list.files(path = 'data_extras/', pattern = standards.pattern))
# filepath <- file.path('data_extras', paste(filename, ".csv", sep = ""))

Ingalls_Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  filter(Column == column) %>%
  dplyr::rename(Precursor.Ion.Name = Compound.Name) %>%
  select(Precursor.Ion.Name,Compound.Type, QE.RF.ratio, Conc..uM, HILICMix, Emperical.Formula) %>%
  filter(!is.na(Conc..uM)) 
Ingalls_Standards$Precursor.Ion.Name <- TrimWhitespace(Ingalls_Standards$Precursor.Ion.Name)

# Import BMIS'd sample file ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'mesocenter_pos_output/', pattern = BMIS_pattern))
filepath <- file.path('mesocenter_pos_output/', paste(filename, ".csv", sep = ""))
BMISd_data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

# Import QC'd files and remove parameter data ------------------------------
filename <- RemoveCsv(list.files(path = 'mesocenter_pos_output/', pattern = QC_pattern))
filepath <- file.path('mesocenter_pos_output/', paste(filename, ".csv", sep = ""))

QCd_data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-10) %>%
  select(-Description, -Value)

# Import Internal standards key ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_extras/', pattern = names_pattern))
filepath <- file.path('data_extras', paste(filename, ".csv", sep = ""))

OG_IS_key <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  dplyr::rename(FinalBMIS = Internal_Standards)

# Apply appropriate filters and isolate standards ---------------------------------------------------------------
Full_data <- QCd_data %>%
  filter(Precursor.Ion.Name %in% Ingalls_Standards$Precursor.Ion.Name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls_Standards, by = "Precursor.Ion.Name") %>%
  select(Replicate.Name, Precursor.Ion.Name, Compound.Type, everything()) %>%
  unique()
