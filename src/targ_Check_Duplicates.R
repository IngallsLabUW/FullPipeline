# Duplicates testing

HILICS.duplicates <- IdentifyDuplicates(QCd.data)

if ("Column" %in% colnames(QCd.data)) {
  duplicates.testing <- QCd.data %>%
    filter(Precursor.Ion.Name %in% HILICS.duplicates$Precursor.Ion.Name) %>%
    group_by(Precursor.Ion.Name, Column) %>%
    mutate(Means = mean(Area.with.QC, na.rm = TRUE)) %>%
    mutate(Std.Devs = sd(Area.with.QC, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Precursor.Ion.Name, Column, Means, Std.Devs) %>%
    unique()

} else {
  print("Non-HILIC data: no duplicates to detect.")
}
