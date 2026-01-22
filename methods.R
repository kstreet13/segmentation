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
# ""
###############################

# refRNA <- rna
# ST <- s1

method2 <- function(refRNA, ST){
    # restrict gene set to intersection
    genes <- rownames(ST)[which(rownames(ST) %in% rownames(refRNA))]
    refRNA <- refRNA[genes, ]
    ST <- ST[genes, ]
    
    # scRNA co-expression
    nz <- assay(refRNA,'counts') > 0
    obsRNA <- tcrossprod(nz) / ncol(refRNA) # how often genes are co-expressed (% of cells)
    obsRNA[which(obsRNA < 0.01/ncol(refRNA))] <- 0.01/ncol(refRNA) # floor (can't take log of 0)    
    expRNA <- outer(rowMeans(nz), rowMeans(nz))
    expRNA[which(expRNA < 0.01/ncol(refRNA))] <- 0.01/ncol(refRNA) # floor (can't divide by 0)
    loRNA <- log(obsRNA) - log(expRNA)
    rm(nz,obsRNA,expRNA)
    
    # ST co-expression
    nz <- assay(ST,'counts') > 0
    obsST <- tcrossprod(nz) / ncol(ST) # how often genes are co-expressed (% of cells)
    obsST[which(obsST < 0.01/ncol(ST))] <- 0.01/ncol(ST) # floor (can't take log of 0)    
    expST <- outer(rowMeans(nz), rowMeans(nz))
    expST[which(expST < 0.01/ncol(ST))] <- 0.01/ncol(ST) # floor (can't divide by 0)
    loST <- log(obsST) - log(expST)
    rm(nz,obsST,expST)
    
    return(cor(as.numeric(loRNA),as.numeric(loST)))
}



x <- as.numeric(loRNA)
y <- as.numeric(loST)
ind <- which(abs(loRNA)>2 | abs(loST)>2)
plot(x[ind],y[ind], col=rgb(0,0,0,.1), asp=1)


