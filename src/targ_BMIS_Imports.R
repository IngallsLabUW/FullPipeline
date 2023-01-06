# Imports for BMIS

# Imports -----------------------------------------------------------------
# Sample Key
SampKey.all <- read.csv("data_extras/Sample.Key.HILIC.csv", stringsAsFactors = FALSE, skip = 1,
                        header = TRUE) %>%
  mutate(Replicate.Name = Sample.Name) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) 

# Internal Standards
Internal.Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  filter(Column == column) %>%
  filter(Compound.Type == "Internal Standard")
Internal.Standards$Compound.Name <- TrimWhitespace(Internal.Standards$Compound.Name)

# QC'd output
filename <- RemoveCsv(list.files(path = "mesocenter_pos_output", pattern = QC_pattern))
filepath <- file.path("mesocenter_pos_output", paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) 

