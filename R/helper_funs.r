test_all_genes <- function(gexp_A, gexp_B, test.alt='greater'){
    # rows = genes, columns = cells 
    # nrows(gexp_A) should be equal to nrows(gexp_B)
    n_genes <- nrow(gexp_A)
    pvals <- sapply(seq(n_genes), function(x) wilcox.test(gexp_A[x,], gexp_B[x,], alternative=test.alt)[['p.value']])
    names(pvals) <- rownames(gexp_A)
#     pvals <- pvals[order(pvals)]
    return(pvals)
}
                    

find.cluster.de.genes <- function(gexp, cell.com, test.alt='two.sided'){
    ## one vs all DE genes via wilcox
    #gexp: rows=genes, columns=cells
    #com: named character vector of cluster memberships
    coms <- unique(cell.com)
    
    com.pvals <- lapply(coms, function(x) {
        curr.com.cells <- names(cell.com)[cell.com == x]
        curr.other.cells <- names(cell.com)[! cell.com == x]
        
        curr.com.gexp <- gexp[, curr.com.cells]
        curr.other.gexp <- gexp[, curr.other.cells]
        
        pvals <- test_all_genes(curr.com.gexp, curr.other.gexp, test.alt=test.alt)
        
        if (length(pvals)==0){
            pvals <- rep(NA, length(nrow(gexp)))
        }
        
        return(pvals)
    })
    
    com.pvals <- Reduce(cbind, com.pvals)
    rownames(com.pvals) <- rownames(gexp)
    colnames(com.pvals) <- coms
    
    return(com.pvals)
}


test_all_genes_stats <- function(gexp_A, gexp_B, test.alt='greater'){
    # rows = genes, columns = cells 
    # nrows(gexp_A) should be equal to nrows(gexp_B)
    n_genes <- nrow(gexp_A)
    stats <- sapply(seq(n_genes), function(x) wilcox.test(gexp_A[x,], gexp_B[x,], alternative=test.alt)[['statistic']])
    names(stats) <- rownames(gexp_A)
#     pvals <- pvals[order(pvals)]
    return(stats)
}

find.cluster.test.stats <- function(gexp, cell.com, test.alt='two.sided'){
    ## one vs all DE genes via wilcox
    #gexp: rows=genes, columns=cells
    #com: named character vector of cluster memberships
    coms <- unique(cell.com)
    
    com.stats <- lapply(coms, function(x) {
        curr.com.cells <- names(cell.com)[cell.com == x]
        curr.other.cells <- names(cell.com)[! cell.com == x]
        
        curr.com.gexp <- gexp[, curr.com.cells]
        curr.other.gexp <- gexp[, curr.other.cells]
        
        stats <- test_all_genes_stats(curr.com.gexp, curr.other.gexp, test.alt=test.alt)
        
        if (length(stats)==0){
            stats <- rep(NA, length(nrow(gexp)))
        }
        
        return(stats)
    })
    
    com.stats <- Reduce(cbind, com.stats)
    rownames(com.stats) <- rownames(gexp)
    colnames(com.stats) <- coms
    
    return(com.stats)
}


find.log.fold.change <- function(gexp, cell.com, log.fn){
    ## one vs all log fold change 
    #gexp: rows=genes, columns=cells
    #com: named character vector of cluster membershps 
    coms <- unique(cell.com)
    
    logfcs <- lapply(coms, function(x){
        curr.com.gexp <- data.frame(gexp[,cell.com==x])
        curr.other.gexp <- data.frame(gexp[,(!cell.com==x)])
        
        curr.com.means <- rowMeans(curr.com.gexp, na.rm = T)
        curr.other.means <- rowMeans(curr.other.gexp, na.rm = T)
        
#         curr.com.means <- sapply(seq(nrow(curr.com.gexp)), function(i) median(as.numeric(curr.com.gexp[i,])))
#         curr.other.means <- sapply(seq(nrow(curr.other.gexp)), function(i) median(as.numeric(curr.other.gexp[i,])))
        
#         print(head(curr.com.means))
#         print(head(curr.other.means))
        curr.lfc <- log.fn(curr.com.means) - log.fn(curr.other.means)
#         print(head(curr.lfc))
        
        curr.lfc
    })
    
    logfcs <- Reduce(cbind, logfcs)
    rownames(logfcs) <- rownames(logfcs)
    colnames(logfcs) <- coms
    
    return(logfcs)
}                         

## jaccard similarity index
jacc.sim <- function(setA, setB){
    ## setA, setB: two sets (character lists)
    
    jacc.sim <- length(intersect(setA, setB))/length(union(setA, setB))
    
    return(jacc.sim)
    
}
                    
find.all.jacc.sims <- function(listA, listB, named=TRUE){
    ### computes jaccard similarity for all pairs  
    ## listA, listB: named list of named lists
    ## named: use names of sublists? 
    
    all.j.sims <- matrix(NA, nrow = length(listA), ncol = length(listB))
    rownames(all.j.sims) <- names(listA)
    colnames(all.j.sims) <- names(listB)
    
    for (i in names(listA)){
        for (j in names(listB)){
            if (named){
                setA <- names(listA[[i]])
                setB <- names(listB[[j]])
            } else if (!named){
                setA <- listA[[i]]
                setB <- listB[[j]]
            }
            
            curr.j.sim <- jacc.sim(setA, setB)
            all.j.sims[i, j] <- curr.j.sim
        }
    }
    
    return(all.j.sims)
}                    
                    

cluster.pseudobulk <- function(gexp, cell.com){
    ## sum single cell expression for each cluster
    #gexp: rows=genes, columns=cells
    #com: named character vector of cluster membershps - order matched to gexp 
    
    all.clusters <- unique(cell.com)
    cluster.counts <- lapply(all.clusters, function(x){
        
        curr.cluster.cells <- names(cell.com)[cell.com==x]
        curr.counts <- rowSums(as.matrix(gexp[,curr.cluster.cells]))
        curr.counts
        
    })
    
    
    cluster.counts <- Reduce(cbind, cluster.counts)
    rownames(cluster.counts) <- rownames(gexp)
    colnames(cluster.counts) <- all.clusters
    
    return(cluster.counts)
}


norm.DESeq <- function(gexp){
    ## return counts normalized by DESeq size factor
    #gexp: rows=genes, columns=cells
    
    ## size factor=median(gene counts / geometric mean(gene counts))
    lc <- log(gexp)
    lc[!is.finite(lc)] <- 0 #replace counts==0 
    loggeommeans <- rowMeans(lc) #geometric means
    
    allZero <- rowSums(gexp) == 0 ## genes with no counts (doesn't happen in merfish data)
    loggeommeans[allZero] <- -Inf
    
    sf <- apply(gexp, 2, function(x){
        exp(median((log(x) - loggeommeans)[is.finite(loggeommeans) & x>0]))
    })
    
#     print(dim(gexp))
    gexp.norm <- t(t(gexp)/sf)
    
    return(list(norm=gexp.norm, sf=sf))
}

run.all.norms <- function(gt_counts, captured_counts=NULL, log_fn=log10, vol=NULL){
    ## gt_counts: ground truth counts all cells 
    ## captured_counts: observed counts (e.g. after adjusting for volume captured and removing empty cells) (e.g. in sim)
    ## if capture_counts is NULL (e.g. using real data), run all norms on gt_counts 
    
    if (is.null(captured_counts)){
#         print('not sim')
        captured_counts <- gt_counts
        genes <- rownames(gt_counts)
        cells <- colnames(gt_counts)
    } else {
        cells <- intersect(colnames(gt_counts), colnames(captured_counts))
        genes <- intersect(rownames(gt_counts), rownames(captured_counts))
        
        gt_counts <- gt_counts[genes, cells]
        captured_counts <- captured_counts[genes, cells]
    }
    
    ## normalize
    # ground truth
    gt <- log_fn(gt_counts + 1)
    gt <- gt[genes,cells]
    # no norm
    nonorm <- log_fn(captured_counts + 1)
    # vol norm
    if (!is.null(vol)){
        volnorm <- t(t(captured_counts)/vol)
        volnorm <- log_fn(volnorm[,cells] + 1)
    } else {
        volnorm <- NA
    }
    
    # libnorm
    libnorm <- t(t(captured_counts)/colSums(captured_counts))*median(colSums(captured_counts), na.rm = T)
    libnorm <- log_fn(libnorm[,cells] + 1)
    # deseq
    deseq <- norm.DESeq(captured_counts)[['norm']]
    deseq <- deseq[,cells]
    
    all.norms <- list(
        gt=gt,
        nonorm=nonorm,
        volnorm=volnorm,
        libnorm=libnorm,
        deseq=deseq)
    
    return(all.norms)
    
}

run.all.pvals <- function(gt_gexp, other_gexp, cell_labs, plot_hist=TRUE, plot_qq=TRUE, plot_cv=TRUE, test.alt='greater'){
    ##gt_gexp: ground truth log normalized gene expression
    ##other_gexp: named list of other gene expression (should be same dim as gt_exp)
    ##cell_labs: named character vector of cell labels (currently supporting two groups)
    
    # get cells in each group
    print(unique(cell_labs))
    cells_A <- names(cell_labs)[cell_labs==unique(cell_labs)[1]]
    print(length(cells_A))
    cells_B <- names(cell_labs)[cell_labs==unique(cell_labs)[2]]
    print(length(cells_B))
    
    # get gt exp of each group
    gt_cells_A <- intersect(cells_A, colnames(gt_gexp))
    gt_cells_B <- intersect(cells_B, colnames(gt_gexp))
    gt_gexp_A <- gt_gexp[,gt_cells_A]
    gt_gexp_B <- gt_gexp[,gt_cells_B]
    
    # get gt pvals
    pvals_gt <- test_all_genes(gt_gexp_A, gt_gexp_B, test.alt = test.alt)
    
    # get other pvals
    pvals_all <- list(gt=pvals_gt)
    for (curr_other in names(other_gexp)){
        print(curr_other)
        curr_other_gexp <- other_gexp[[curr_other]]
        
        curr_cells_A <- intersect(cells_A, colnames(curr_other_gexp))
        curr_cells_B <- intersect(cells_B, colnames(curr_other_gexp))
        
        curr_other_gexp_A <- curr_other_gexp[,curr_cells_A]
        curr_other_gexp_B <- curr_other_gexp[,curr_cells_B]
        
        curr_pval <- test_all_genes(curr_other_gexp_A, curr_other_gexp_B, test.alt = test.alt)
        
        pvals_all[[curr_other]] <- curr_pval
    }
    
    #plot histograms of pvalues
    if (plot_hist){
        par(mfrow = c(1, length(pvals_all)))
        
        for (i in seq(length(pvals_all))){
            curr.norm <- names(pvals_all)[i]
            curr.norm.pvals <- pvals_all[[curr.norm]]
            
            hist(curr.norm.pvals, breaks=100,
                main = paste(curr.norm, 'pvals'))
            
        }
    }
    
    
    #plot QQ plots
    if (plot_qq){
        
        par(mfrow = c(1, length(other_gexp)))
        for (i in seq(length(other_gexp))){
            gt <- pvals_all[['gt']]
            other_name <- names(other_gexp)[i]
            other <- pvals_all[[other_name]]

            plot(gt[order(gt)], other[order(other)], 
                 main = paste(other_name, 'vs ground truth pvals'))
            abline(a=0,b=1, col='red', lwd=2)
        }
    }
    
    
    #plot CVs
    if (plot_cv){
        
        par(mfrow = c(1, length(other_gexp)))
        for (i in seq(length(other_gexp))){
            gt.sd <- apply(gt_gexp, 1, sd)
            gt.mean <- rowMeans(gt_gexp)
            gt.cv2 <- (gt.sd^2/gt.mean)
            gt.cv2 <- (gt.cv2-min(gt.cv2, na.rm=T))/(max(gt.cv2, na.rm=T)-min(gt.cv2, na.rm=T))
            
            other_name <- names(other_gexp)[i]
            other.sd <- apply(other_gexp[[other_name]], 1, sd)
            other.mean <- rowMeans(other_gexp[[other_name]])
            other.cv2 <- (other.sd^2/other.mean)
            other.cv2 <- (other.cv2-min(other.cv2, na.rm=T))/(max(other.cv2, na.rm=T)-min(other.cv2, na.rm=T))

            plot(gt.cv2, other.cv2,
                 main = paste(other_name, 'vs ground truth CV^2 (scaled)'))
            abline(a=0,b=1, col='red')
        }
    }
    
    
    return(pvals_all)  
             
}


scale.range <- function(vec){
    ## scales numbers in vector so that they're between 0 and 1
    vec.min <- min(vec, na.rm = T)
    vec.max <- max(vec, na.rm = T)
    
    vec.scaled <- (vec - vec.min)/(vec.max - vec.min)
    return(vec.scaled)
}
                    
                    
run.all.cvs <- function(all.norms){
    ## all.norms: named list of gene expression matrices w different normalizations
    
    all.cvs <- lapply(all.norms, function(gexp){
        gene.sd <- apply(gexp, 1, sd)
        gene.mean <- rowMeans(gexp)
        gene.cv2 <- (gene.sd^2/gene.mean)
#         gene.cv2 <- (gene.cv2-min(gene.cv2, na.rm=T))/(max(gene.cv2, na.rm=T)-min(gene.cv2, na.rm=T))
        names(gene.cv2) <- rownames(gexp)
        gene.cv2
    })
    
#     all.cvs <- Reduce(cbind, all.cvs)
    
    return(all.cvs)
}

                    
                    
rmse <- function(x,y) {
    #x and y are ordered vectors with corresponding pairs of values 
    
    inf.idx <- union(which(is.infinite(x)), which(is.infinite(y)))
    if (length(inf.idx)>0){
        x <- x[-inf.idx]
        y <- y[-inf.idx]
    }

    err <- x - y
    err.sq <- err^2
    mean.sq.err <- sum(err.sq)/length(err)
    rmse <- sqrt(mean.sq.err)
    return(rmse)
}                

                    
get.lfcs.fp.fn <- function(testing, ground.truth) {
    gt <- sign(ground.truth)
    x <- sign(testing)
    
    ## false positive: negative in gt and positive in x
    fp <- (gt == -1) & (x == 1)
    ## false negative: positive in gt and negative in x
    fn <- (gt == 1) & (x == -1)
    
    ## true positive: positive in x and positive in gt 
    tp <- (gt == 1) & (x == 1)
    ## true negative: negative in x and negative in gt 
    tn <- (gt == -1) & (x == -1)
    
    false <- list(false.positive = fp, false.negative = fn, true.positive = tp, true.negative = tn)
    return(false)
}
                    
get.pval.fp.fn <- function(testing, ground.truth, thresh) {
    gt <- ground.truth < thresh
    x <- testing < thresh
    
    ##false positive: positive in x (true) and negative in gt (false)
    fp <- x & !gt
    ##false negative: negative in x (false) and positive in gt (true)
    fn <- !x & gt
    
    ##true positive: positive in x (true) and positive in gt (true)
    tp <- x & gt
    ##true negative: negative in x (false) and negative in gt (false)
    tn <- !x & !gt
    
    false <- list(false.positive = fp, false.negative = fn, true.positive = tp, true.negative = tn)
    return(false)    
}                    
                    
