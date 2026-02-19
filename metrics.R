require(Seurat)
rna <- readRDS('data/20251119_3922.rds')
require(SingleCellExperiment)
rna <- as.SingleCellExperiment(rna)

# get Ensembl IDs
require(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
gene_symbols <- rownames(rna)
results <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                 filters = "external_gene_name",
                 values = gene_symbols,
                 mart = ensembl)
rowData(rna)$ensemblID <- results$ensembl_gene_id[match(rownames(rna), results$external_gene_name)]
rowData(rna)$ensemblID[grep('ENSMUSG', rownames(rna))] <- rownames(rna)[grep('ENSMUSG', rownames(rna))]
rm(ensembl, gene_symbols, results)
# switch rownames to Ensembl IDs
rowData(rna)$Symbol <- rownames(rna)
rownames(rna) <- rowData(rna)$ensemblID


# read spatial data
require(DropletUtils)
# doesn't actually include locations
s1 <- read10xCounts("data/xenium/output-XETG00402__0054800__Region_1__20250822__221946/cell_feature_matrix")
###
# 6 types of gene: Gene Expression, Negative Control Probe, Genomic Control, Negative Control Codeword, Unassigned Codeword, Deprecated Codeword
###
s1 <- s1[which(rowData(s1)$Type == 'Gene Expression'), ]


###############################
# metric ideas
# (1) for all pairs of genes that are never co-expressed in scRNA, how often are they co-expressed in ST?
# (2) continuous measure of co-expression: log(obs% / exp%), pearson correlation between scRNA and ST
# (3) cosine distance + KL divergence?
# (4) integrate (StabMap/WNN), then batch integration metrics
###############################


###############################
# (1) for all pairs of genes that are never co-expressed in scRNA, how often are they co-expressed in ST?
###############################



