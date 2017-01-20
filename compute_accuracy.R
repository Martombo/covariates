require(limma)
require(edgeR)
require(DESeq2)
require(sva)


# voom DEA function
voom_res = function(counts, group, design, n_top=100000){
	dge = DGEList(counts=counts, group=group)
	dge = calcNormFactors(dge)
	v = voom(dge, design)
	fit = lmFit(v, design)
	fit = eBayes(fit)
	results = topTable(fit, n=n_top, coef="groupB")
	return(results[order(row.names(results)),])
}

# DESeq2 DEA function
deseq2_res = function(counts, group, design){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design=~group)
	dds = DESeq(dds)
	results = results(dds)
}

# true positive, false positive rates function
aucc = function(res){
	fps = rep(1, dim(res)[1])
	fps[which(res$logFC * res$true_log2FC > 0)] = 0
	tps = -(fps - 1)
	res$fpr = cumsum(fps) / sum(fps)
	res$tpr = cumsum(tps) / sum(tps)
	return(res)
}

# sva correction function (based on counts)
sva_corr = function(counts, group, n_sv=0){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design=~group)
	dds = estimateSizeFactors(dds)
	dat = counts(dds, normalized=T)
	dat = dat[rowMeans(dat)>0,]
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

compute_accuracy = function(n_simulations, group, corr=F, n_sv=0){
	for (n_sim in seq(n_simulations)){
		counts = read.table(paste0("counts", n_sim))
		design = model.matrix(~group)
		if (corr){
			svs = sva_corr(counts, group, n_sv)
			if (length(svs$sv) > 1){
				design = model.matrix(~svs$sv+group)
			}
		}
		results = voom_res(counts, group, design)
		log2FCs = read.table(paste0("log2FCs", n_sim))
		names(log2FCs)[1] = "true_log2FC"
		results = merge(results, log2FCs, by=0)
		results$Row.names = NULL
		results = results[order(results$adj.P.Val), ]
		results = aucc(results)
		return(results)
	}
}
