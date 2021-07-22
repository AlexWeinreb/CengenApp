library(Seurat)

full_data <- load("../CengenApp/Dataset_6July_2021.rda")

rm(list = c("allCells", "allNeurons"))

save(list = ls(), file = "DataLite_6Jult_2021.rda")
