############################
####
#### R Script to Filter Data
####
############################

cat("\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering Individuals %%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")

cat(paste0("Age Range: ", stringDate, "\n\n"))

## Set-up
invisible(library(tidyverse))

## Load data set
anno <- read.delim("Data/anno.txt")
dim(anno)

## Remove samples in Genetic ID containing strings like ".SG", ".DG", "_contam", "_d", and other weird labels
remove.list <- paste(c(".SG", ".DG", "contam", "_d"), collapse = "|")
anno <- anno %>% filter(!str_detect(Genetic.ID, remove.list))
dim(anno)

## Remove samples that have familial status in the Group_ID
remove.list <- paste(c("brother", "sister", "father", "1d_degree_rel", "contam", "_o", "mother", "1d.rel", "son", "sibling"), collapse = "|") #, "daughter"
anno <- anno %>% filter(!str_detect(Group.ID, remove.list))
dim(anno)

## Remove rows with ".." in lat./long.
removeIndices <- which(anno$Lat. == "..")
anno <- anno %>% filter(!row_number() %in% removeIndices)
dim(anno)

## Ensure samples are processed on the 1240k library
removeIndices <- which(anno[, 64] == "..")
anno <- anno %>% filter(!row_number() %in% removeIndices)
dim(anno)

## Pick Countries
keep.list <- paste(c("Hungary", "Croatia", "Serbia", "Romania", "North Macedonia", "Germany", "Bulgaria", "Ukraine", "United Kingdom", "Greece", "Turkey", "France", "Italy", "Austria", "Switzerland", "Slovakia", "Czech Republic", "Belgium", "Netherlands", "Spain", "Portugal", "Denmark", "Norway", "Sweden", "Poland", "Belarus", "Moldova", "Lithuania, Latvia", "Estonia", "Albania", "Montenegro", "Kosovo", "Slovenia"), collapse = "|")
anno <- anno %>% filter(str_detect(Country, keep.list))
dim(anno)

## Pick Dates
names(anno)[11] <- "Date.BP"
range <- endDate:startDate
anno <- anno %>% filter(Date.BP %in% range)

## Write file
write.csv(anno, "Data/anno_R.csv")
cat("\nFiltered individuals.")
cat("\n\nNumber of Populations: ", length(unique(anno$Group.ID)))
cat("\n\nNumber of Individuals: ", length(unique(anno$Genetic.ID)))