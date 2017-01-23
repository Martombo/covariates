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

# true positive, false positive rates function
true_false_pos = function(res){
	fps = rep(1, dim(res)[1])
	fps[which(res$logFC * res$true_log2FC > 0)] = 0
	tps = -(fps - 1)
	res$fpr = cumsum(fps) / sum(fps)
	res$tpr = cumsum(tps) / sum(tps)
	return(res)
}

# AUC computation function 
auc = function(res, fpr_limit=0.1){
	n_genes = dim(res)[1]
	res = res[which(res$fpr < fpr_limit), ]
	auc = sum(res$tpr) / (fpr_limit * n_genes)
	return(auc)
}

# sva correction function (based on counts)
sva_corr = function(counts, group, n_sv=0){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design)
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
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design)
	dds = estimateSizeFactors(dds)
	dat = rlogTransformation(dds)
	mod = model.matrix(~group, colData(dds))
	mod0 = model.matrix(~1, colData(dds))
	svseq = svaseq(dat, mod, mod0)
	return(svseq)
}

### function to run the differential expression analysis on one dataset
	# includes the computation of true and false positive rates
perform_dea = function(n_sim, group, corr=0, n_sv=NA){
	counts = read.table(paste0("counts", n_sim))
	design = model.matrix(~group)
	if (!is.na(n_sv)){
		svs = sva_corr(counts, group, n_sv)
		corr = svs$sv
	}
	if (length(corr) > 1){
		design = model.matrix(~as.double(corr)+group)
	}
	results = voom_res(counts, group, design)
	log2FCs = read.table(paste0("log2FCs", n_sim))
	names(log2FCs)[1] = "true_log2FC"
	results = merge(results, log2FCs, by=0)
	results$Row.names = NULL
	results = results[order(results$adj.P.Val), ]
	results = true_false_pos(results)
	return(results)
}

### AUC computation on list of simulated datasets ###
compute_accuracy = function(n_simulations, n_replicates, corr=0, n_sv=NA){
	if (length(n_replicates) == 1) {
		n_replicates = rep(n_replicates, 2)
	}
	group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
	aucs = c()
	for (n_sim in seq(n_simulations)){
		results = perform_dea(n_sim, group, corr, n_sv)
		aucs = c(aucs, auc(results))
	}
	return(aucs)
}
