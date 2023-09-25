library(Seurat)
library(tidyverse)



merged_file <- readRDS("/Users/hailey/Downloads/merfish/Merged_SVZ_Seurat_dims14_Labelled_GOOD.rds")

metadata = merged_file@meta.data


#for cell type composition only
metadata = merged_file@meta.data

cell_counts <- metadata %>%
  group_by(orig.ident, CellType, Condition) %>%
  summarise(cell_count = n()) %>%
  ungroup()


cell_counts <- cell_counts %>%
  pivot_wider(names_from = CellType, values_from = cell_count, 
              id_cols = c(orig.ident, Condition))


write.csv(cell_counts, "/Users/hailey/Downloads/merfish/celltype_counts.csv", row.names = TRUE)




path = "/Users/hailey/Downloads/merfish/"

#####Get metadata for each sample#####
get_sample_metadata <- function(metadata, sample, path) {
  
  meta <- subset(metadata, orig.ident == sample)
  
  meta$library_id <- sample
  meta$patient <- sample
  
  meta_name = paste(sample, "_meta.csv", sep = "")
  
  full_path <- file.path(path, meta_name)
  
  write.csv(meta, full_path, row.names = TRUE)
  
  return(meta)
}

NX1_meta <- get_sample_metadata(metadata, "NX1", path)
NX2_meta <- get_sample_metadata(metadata, "NX2", path)
NX3_meta <- get_sample_metadata(metadata, "NX3", path)
HX1_meta <- get_sample_metadata(metadata, "NX3", path)
HX2_meta <- get_sample_metadata(metadata, "NX3", path)
HX3_meta <- get_sample_metadata(metadata, "NX3", path)



#####Get meta data for each sample#####


NX1_meta <- subset(metadata, orig.ident == "NX1")
NX1_meta$library_id <- "NX1"
NX1_meta$patient <- "m1"

NX1_meta_cell <- rownames(NX1_meta)
cell_type <- NX1_meta$CellType
new_df <- data.frame(Column1 = NX1_meta_cell, Column2 = cell_type)
names(new_df) <- c("X", "x")

NX2_meta <- subset(metadata, orig.ident == "NX2")
NX2_meta$library_id <- "NX2"
NX2_meta$patient <- "m2"


NX3_meta <- subset(metadata, orig.ident == "NX3")
NX3_meta$library_id <- "NX3"
NX3_meta$patient <- "m3"


HX1_meta <- subset(metadata, orig.ident == "HX1")
HX1_meta$library_id <- "HX1"
HX1_meta$patient <- "m1"


HX2_meta <- subset(metadata, orig.ident == "HX2")
HX2_meta$library_id <- "HX2"
HX2_meta$patient <- "m2"


HX3_meta <- subset(metadata, orig.ident == "HX3")
HX3_meta$library_id <- "HX3"
HX3_meta$patient <- "m3"


write.csv(new_df, "/Users/hailey/Downloads/merfish/NX1_celltype.csv", row.names = TRUE)


write.csv(HX1_meta, "/Users/hailey/Downloads/merfish/HX1_meta.csv", row.names = TRUE)
write.csv(HX2_meta, "/Users/hailey/Downloads/merfish/HX2_meta.csv", row.names = TRUE)
write.csv(HX3_meta, "/Users/hailey/Downloads/merfish/HX3_meta.csv", row.names = TRUE)
write.csv(NX1_meta, "/Users/hailey/Downloads/merfish/NX1_meta.csv", row.names = TRUE)
write.csv(NX2_meta, "/Users/hailey/Downloads/merfish/NX2_meta.csv", row.names = TRUE)
write.csv(NX3_meta, "/Users/hailey/Downloads/merfish/NX3_meta.csv", row.names = TRUE)

#combined_meta <- rbind(NX1_meta, NX2_meta, NX3_meta)

#subset counts
count_sparse <-  merged_file@assays[["SCT"]]@counts
counts_df <- as.data.frame(count_sparse)

NX1_counts <- counts_df[, grep("NX1", colnames(counts_df))]
dim(NX1_counts)
NX1_meta_cell <- rownames(NX1_meta)
NX1_count_final <- NX1_counts[, colnames(NX1_counts) %in% NX1_meta_cell]
NX1_transposed <- t(NX1_count_final)


NX2_counts <- counts_df[, grep("NX2", colnames(counts_df))]
NX2_meta_cell <- rownames(NX2_meta)
NX2_count_final <- NX2_counts[, colnames(NX2_counts) %in% NX2_meta_cell]
NX2_transposed <- t(NX2_count_final)


NX3_counts <- counts_df[, grep("NX3", colnames(counts_df))]
NX3_meta_cell <- rownames(NX3_meta)
NX3_count_final <- NX3_counts[, colnames(NX3_counts) %in% NX3_meta_cell]
NX3_transposed <- t(NX3_count_final)


HX1_counts <- counts_df[, grep("HX1", colnames(counts_df))]
HX1_meta_cell <- rownames(HX1_meta)
HX1_count_final <- HX1_counts[, colnames(HX1_counts) %in% HX1_meta_cell]
HX1_transposed <- t(HX1_count_final)


HX2_counts <- counts_df[, grep("HX2", colnames(counts_df))]
HX2_meta_cell <- rownames(HX2_meta)
HX2_count_final <- HX2_counts[, colnames(HX2_counts) %in% HX2_meta_cell]
HX2_transposed <- t(HX2_count_final)


HX3_counts <- counts_df[, grep("HX3", colnames(counts_df))]
HX3_meta_cell <- rownames(HX3_meta)
HX3_count_final <- HX3_counts[, colnames(HX3_counts) %in% HX3_meta_cell]
HX3_transposed <- t(HX3_count_final)


write.csv(NX1_count_final, "/Users/hailey/Downloads/merfish/count_NX1_untrans.csv", row.names = TRUE)


write.csv(HX1_transposed, "/Users/hailey/Downloads/merfish/count_HX1.csv", row.names = TRUE)
write.csv(HX2_transposed, "/Users/hailey/Downloads/merfish/count_HX2.csv", row.names = TRUE)
write.csv(HX3_transposed, "/Users/hailey/Downloads/merfish/count_HX3.csv", row.names = TRUE)
write.csv(NX1_transposed, "/Users/hailey/Downloads/merfish/count_NX1.csv", row.names = TRUE)
write.csv(NX2_transposed, "/Users/hailey/Downloads/merfish/count_NX2.csv", row.names = TRUE)
write.csv(NX3_transposed, "/Users/hailey/Downloads/merfish/count_NX3.csv", row.names = TRUE)


#subset coordinates

coors_all_NX1 <- merged_file@images[["NX1_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_NX1 <- as.data.frame(coors_all_NX1)
cell_names_coor_NX1 = merged_file@images[["NX1_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_NX1) <- cell_names_coor_NX1

coors_all_NX2 <- merged_file@images[["NX2_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_NX2 <- as.data.frame(coors_all_NX2)
cell_names_coor_NX2 = merged_file@images[["NX2_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_NX2) <- cell_names_coor_NX2

coors_all_NX3 <- merged_file@images[["NX3_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_NX3 <- as.data.frame(coors_all_NX3)
cell_names_coor_NX3 = merged_file@images[["NX3_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_NX3) <- cell_names_coor_NX3

coors_all_HX3 <- merged_file@images[["HX3_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_HX3 <- as.data.frame(coors_all_HX3)
cell_names_coor_HX3 = merged_file@images[["HX3_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_HX3) <- cell_names_coor_HX3

coors_all_HX2 <- merged_file@images[["HX2_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_HX2 <- as.data.frame(coors_all_HX2)
cell_names_coor_HX2 = merged_file@images[["HX2_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_HX2) <- cell_names_coor_HX2

coors_all_HX1 <- merged_file@images[["HX1_ALL"]]@boundaries[["centroids"]]@coords
coordinates_df_HX1 <- as.data.frame(coors_all_HX1)

cell_names_coor_HX1 = merged_file@images[["HX1_ALL"]]@boundaries[["centroids"]]@cells
rownames(coordinates_df_HX1) <- cell_names_coor_HX1

#combined_coor <- rbind(coordinates_df_NX1, coordinates_df_NX2, coordinates_df_NX3)


write.csv(coordinates_df_HX1, "/Users/hailey/Downloads/merfish/coor_HX1.csv", row.names = TRUE)
write.csv(coordinates_df_HX2, "/Users/hailey/Downloads/merfish/coor_HX2.csv", row.names = TRUE)
write.csv(coordinates_df_HX3, "/Users/hailey/Downloads/merfish/coor_HX3.csv", row.names = TRUE)
write.csv(coordinates_df_NX1, "/Users/hailey/Downloads/merfish/coor_NX1.csv", row.names = TRUE)
write.csv(coordinates_df_NX2, "/Users/hailey/Downloads/merfish/coor_NX2.csv", row.names = TRUE)
write.csv(coordinates_df_NX3, "/Users/hailey/Downloads/merfish/coor_NX3.csv", row.names = TRUE)


