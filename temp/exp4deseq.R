library(readr)
V.HMC.map <- read_delim(
  "C:/Users/mark/Desktop/map.txt",
  "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

V.HMC.map <- V.HMC.map[, c("#SampleID",
                           "treatment", "day")]

names(V.HMC.map) <- c("SampleID", "treatment", "days.string")

temp <- strsplit(V.HMC.map$days.string, "\\+|h")
temp <- plyr::ldply(temp, rbind)
temp <- as.data.frame(temp)
temp[] <- lapply(temp[], as.character)
temp[] <- lapply(temp[], as.numeric)
names(temp) <- c("days", "hours")
temp$hours[is.na(temp$hours)] <- 0
V.HMC.map <- cbind.data.frame(V.HMC.map, temp)
V.HMC.map$days.hours <- temp$days + (temp$hours / 24)

biom_table_summary <-
  read_delim(
    "C:/Users/mark/Desktop/biom_table_summary.txt",
    " ",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE,
    skip = 15
  )

biom_table_summary$X1 <-
  sub(pattern = "\\:$", "", biom_table_summary$X1)
biom_table_summary[] <- lapply(biom_table_summary[], as.character)
biom_table_summary[] <- lapply(biom_table_summary[], as.numeric)

combo <- merge(
  x = V.HMC.map,
  y = biom_table_summary,
  by.x  = "SampleID",
  by.y = "X1",
  all = TRUE
)

table(V.HMC.map$treatment, V.HMC.map$days.hours)
