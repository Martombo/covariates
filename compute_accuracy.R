require(limma)
require(edgeR)
require(DESeq2)
require(sva)


# voom DEA function
voom_res = function(counts, group, design=NA, n_top=100000){
	if (is.na(design)){
		design = model.matrix(~group)
	}
	dge = DGEList(counts=counts, group=group)
	dge = calcNormFactors(dge)
	v = voom(dge, design)
	fit = lmFit(v, design)
	fit = eBayes(fit)
	results = topTable(fit, n=n_top, coef=2)
	return(results[order(row.names(results)),])
}

# true positive, false positive per row function
calc_pr = function(x){
	if (x[1] * x[7] <= 0){
		assign("fp", fp+1, 1)
	}else{
		assign("tp", tp+1, 1)
	}
	return(c(tp, fp))
}

# true positive, false positive rates computation function
aucc = function(res, n_genes, n_DEGs){
	assign("fp", 0, 1)
	assign("tp", 0, 1)
	pr = data.frame(t(apply(res, 1, calc_pr)))
	res$tpr = pr$X1 / n_DEGs
	res$fpr = pr$X2 / (n_genes - n_DEGs)
	return(res)
}

# sva correction function (based on counts)
sva_corr = function(counts, group, n_sv){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design=~group)
	dds = estimateSizeFactors(dds)
	dat = counts(dds, normalized=T)
	dat = dat[rowMeans(dat)>5,]
	mod = model.matrix(~group, colData(dds))
	mod0 = model.matrix(~1, colData(dds))
	svseq = svaseq(dat, mod, mod0, n.sv=n_sv)
	return(svseq)
}

# sva correction function (based on DESeq2 rld)
sva_corr_rld = function(counts, group){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design=~group)
	dds = estimateSizeFactors(dds)
	dat = rlogTransformation(dds)
	mod = model.matrix(~group, colData(dds))
	mod0 = model.matrix(~1, colData(dds))
	svseq = svaseq(dat, mod, mod0)
	return(svseq)
}


### run sva estimation ### #to be deleted?
estimate_sv = function(n_simulations, n_replicates=c(3,3), n_DEGs=3000, n_top=1000, n_sv=0){ 

	if (length(n_replicates) == 1) {
		n_replicates = rep(n_replicates, 2)
	}

	svs = list()
	group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
	for (n_simul in seq(n_simulations)){
		counts = read.table(paste0("counts", n_simul))
		svs[[n_simul]] = sva_corr(counts, group, n_sv=n_sv)
		n_genes = dim(counts)[1]
	}
	return(svs)
}

compute_accuracy = function(n_simulations, group, n_DEGs, design=NA){
	for (n_sim in seq(n_simulations)){
		counts = read.table(paste0("counts", n_sim))
		results = voom_res(counts, group, design)
		log2FCs = read.table(paste0("log2FCs", n_sim))
		names(log2FCs)[1] = "true_log2FC"
		results = merge(results, log2FCs, by=0)
		results$Row.names = NULL
		results = results[order(results$adj.P.Val), ]
		results = aucc(results, n_genes, n_DEGs)
		return(results)
	}
}
