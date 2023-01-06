## Skyline Rearrange and Compound Name Check

# Remove any illegal values
ReplaceNonvalues <- function(x) (gsub("#N/A", NA, x))

if (polarity == "pos") {
  
  skyline_HILIC_pos <- skyline_HILIC_pos %>%
    mutate(Column = "HILICpos")
  skyline_HILIC_neg <- skyline_HILIC_neg %>%
    mutate(Column = "HILICneg")
  
  combined_skyline <- skyline_HILIC_pos %>%
    rbind(skyline_HILIC_neg) %>%
    select(Replicate.Name, Precursor.Ion.Name, Retention.Time, Area, Background, Height, Mass.Error.PPM, Column) %>%
    mutate_all(ReplaceNonvalues) 
  
  # Change variable classes
  skyline_classes_changed <- ChangeClasses(combined_skyline, start.column = 3, end.column = 7) 
  
} 

