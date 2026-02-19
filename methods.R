###############################
# (1) for all pairs of genes that are never co-expressed in scRNA, how often are
# they co-expressed in ST?
# "spurious co-expressions per cell"
###############################

method1 <- function(refRNA, ST){
    # restrict gene set to intersection
    genes <- rownames(ST)[which(rownames(ST) %in% rownames(refRNA))]
    refRNA <- refRNA[genes, ]
    ST <- ST[genes, ]
    
    # minimum expression threshold?
    
    # scRNA co-expression
    nz <- assay(refRNA,'counts') > 0
    coexRNA <- tcrossprod(nz)
    coexRNA <- coexRNA > 0 # binary, co-expressed or not
    diag(coexRNA) <- 1
    rm(nz)
    
    # ST co-expression
    nz <- assay(ST,'counts') > 0
    coexST <- tcrossprod(nz) # how often genes are co-expressed
    rm(nz)
    
    stopifnot(all(rownames(coexRNA) == rownames(coexST)))
    stopifnot(all(colnames(coexRNA) == colnames(coexST)))
    
    # count spurious co-expression events
    # normalize by cell count (put over-segmentation and under-segmentation on equal footing. it's not perfect, but it's potentially informative)
    return(sum(coexST[which(coexRNA==0)]) / ncol(ST))
}







###############################
# (2) continuous measure of co-expression: log(obs% / exp%), pearson correlation
# between scRNA and ST
# "co-expression correlation"?
###############################

# refRNA <- rna
# ST <- s1

method2 <- function(refRNA, ST){
    # restrict gene set to intersection
    genes <- rownames(ST)[which(rownames(ST) %in% rownames(refRNA))]
    refRNA <- refRNA[genes, ]
    
    # only keep genes expressed in at least 20 cells in RNA
    nzRNA <- assay(refRNA,'counts') > 0
    keep <- which(rowSums(nzRNA) > 20)
    genes <- genes[keep]
    
    refRNA <- refRNA[genes, ]
    ST <- ST[genes, ]
    
    # scRNA co-expression
    nzRNA <- assay(refRNA,'counts') > 0
    obsRNA <- tcrossprod(nzRNA) / ncol(refRNA) # how often genes are co-expressed (% of cells)
    obsRNA[which(obsRNA < 0.01/ncol(refRNA))] <- 0.01/ncol(refRNA) # floor (can't take log of 0)    
    expRNA <- outer(rowMeans(nzRNA), rowMeans(nzRNA))
    expRNA[which(expRNA < 0.01/ncol(refRNA))] <- 0.01/ncol(refRNA) # floor (can't divide by 0)
    loRNA <- log(obsRNA) - log(expRNA)
    diag(loRNA) <- NA
    rm(nzRNA,obsRNA,expRNA)
    
    # ST co-expression
    nzST <- assay(ST,'counts') > 0
    obsST <- tcrossprod(nzST) / ncol(ST) # how often genes are co-expressed (% of cells)
    obsST[which(obsST < 0.01/ncol(ST))] <- 0.01/ncol(ST) # floor (can't take log of 0)    
    expST <- outer(rowMeans(nzST), rowMeans(nzST))
    expST[which(expST < 0.01/ncol(ST))] <- 0.01/ncol(ST) # floor (can't divide by 0)
    loST <- log(obsST) - log(expST)
    diag(loST) <- NA
    rm(nzST,obsST,expST)
    
    return(cor(as.numeric(loRNA),as.numeric(loST), use = 'complete.obs'))
}


sparse.cor <- function(x){
    n <- nrow(x)
    m <- ncol(x)
    ii <- unique(x@i)+1 # rows with a non-zero element
    
    Ex <- colMeans(x)
    nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans
    
    covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                    crossprod(t(Ex))*(n-length(ii))
    )/(n-1)
    sdvec <- sqrt(diag(covmat))
    covmat/crossprod(t(sdvec))
}
method2.5 <- function(refRNA, ST){
    # normalize
    refRNA <- Seurat::LogNormalize(data = assay(refRNA,'counts'))
    ST <- Seurat::LogNormalize(data = assay(ST,'counts'))
    
    # restrict gene set to intersection
    genes <- rownames(ST)[which(rownames(ST) %in% rownames(refRNA))]
    refRNA <- refRNA[genes, ]
    
    # only keep genes expressed in at least 100 cells in RNA
    nzRNA <- refRNA > 0
    keep <- which(rowSums(nzRNA) > 100)
    genes <- genes[keep]
    rm(nzRNA)
    
    refRNA <- refRNA[genes, ]
    ST <- ST[genes, ]
    
    # scRNA co-expression
    # corRNA <- cor(t(assay(refRNA,'logcounts'))) # can't do, no cor() for sparse matrix
    refRNA <- t(refRNA)
    corRNA <- sparse.cor(refRNA)
    corRNA[lower.tri(corRNA)] <- NA
    diag(lower.tri(corRNA)) <- NA

    # ST co-expression
    ST <- t(ST)
    corST <- sparse.cor(ST)
    corST[lower.tri(corST)] <- NA
    diag(lower.tri(corST)) <- NA
    
    
    return(cor(as.numeric(loRNA),as.numeric(loST), use = 'complete.obs'))
}


###############################
# # (3) project to circle (cosine distance?) + KL divergence?
# ""
###############################

method3 <- function(refRNA, ST){
    # restrict gene set to intersection
    genes <- rownames(ST)[which(rownames(ST) %in% rownames(refRNA))]
    refRNA <- refRNA[genes, ]
    ST <- ST[genes, ]
    
    # project onto circle
    sizes <- sqrt(colSums(assay(refRNA,'logcounts')^2))
    circRNA <- t(t(assay(refRNA,'logcounts')) / sizes)
    sizes <- sqrt(colSums(assay(refST,'logcounts')^2))
    circST <- t(t(assay(refRNA,'logcounts')) / sizes)
    
    
    

    return()
}





###############################
# # (4) integrate (StabMap/WNN), then batch integration metrics
# ""
###############################


m4_integrate <- function(refRNA, ST){
    require(StabMap)
    require(Seurat)
    require(matrixStats)
    rv <- rowVars(assay(refRNA,'logcounts'))
    keep <- which(rv >= sort(rv, decreasing = TRUE)[2000] |
                      rownames(refRNA) %in% rownames(ST))
    rm(rv)
    
    assay_list <- list(rna = assay(refRNA,'logcounts')[keep,],
                       st = assay(ST,'logcounts'))
    # get joint embedding
    stab <- stabMap(assay_list, reference_list = c("rna"),
                    suppressMessages = FALSE, maxFeatures = nrow(assay_list$rna),
                    plot = FALSE)
    mod <- factor(rep(c('rna','st'),
                      times = c(ncol(assay_list$rna), ncol(assay_list$st))))
    
    return(list(x = stab, label = mod))
}






