require(Seurat)
rna <- readRDS('data/20251119_3922.rds')


# get spatial data
require(rhdf5)
s1 <- h5read("data/xenium/output-XETG00402__0054800__Region_1__20250822__221946/cell_feature_matrix.h5","/matrix")
# 13431 155086


require(DropletUtils)
read10xCounts("data/xenium/output-XETG00402__0054800__Region_1__20250822__221946/cell_feature_matrix")


# metric ideas
# (1) for all pairs of genes that are never co-expressed in scRNA, how often are they co-expressed in ST?
# (2) continuous measure of co-expression: log(obs% / exp%), pearson correlation between scRNA and ST
# (3) cosine distance + KL divergence?
# (4) integrate (StabMap/WNN), then batch integration metrics