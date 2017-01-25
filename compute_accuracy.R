# voom DEA function
voom_res = function(counts, group, design, n_top=100000){
	require(limma)
	require(edgeR)
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

# RUVseq empirical control genes
ruv_corr_g = function(counts, group, n_sv=0){
	require(RUVSeq)
	require(DESeq2)
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), ~group)
	dds = DESeq(dds)
	res = as.data.frame(results(dds))
	res = res[order(res$padj),]
	gene_limit = round(dim(res)[1] * 0.2)
	control_genes = row.names(res)[gene_limit:dim(res)[1]]
	counts = as.matrix(counts)
	row.names(counts) = seq(dim(counts)[1])
	set = RUVg(counts, control_genes, k=n_sv)
	return(set[[1]])
}

# sva correction function (based on counts)
sva_corr = function(counts, group, n_sv=0){
	require(sva)
	require(DESeq2)
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), ~group)
	dds = estimateSizeFactors(dds)
	dat = counts(dds, normalized=T)
	dat = dat[rowMeans(dat)>5,]
	mod = model.matrix(~group, colData(dds))
	mod0 = model.matrix(~1, colData(dds))
	svseq = svaseq(as.matrix(dat), mod, mod0, n.sv=n_sv)
	return(svseq)
}


# function to determine design. estimates number of covariates, if not already given
determine_design = function(counts, group, sva_method, n_sv){
	require(sva)
	design = model.matrix(~group)
	if (is.null(n_sv)){
		n_sv = num.sv(as.matrix(counts), model.matrix(~group), sva_method)
	}
	if (n_sv > 0){
		if (!is.null(sva)){
			corr = sva_corr(counts, group, n_sv)$sv
		}else{
			if (ruv == "g"){
				corr = ruv_corr_g(counts, group, n_sv)
			}
		}
		corr = as.double(corr)
		if (is.null(dim(corr)[1])){
			design = model.matrix(~corr+group)
		}
		if (dim(corr)[1] == 2){
			design = model.matrix(~corr[,1]+corr[,2]+group)
		}else{
			design = model.matrix(~corr[,1]+corr[,2]+corr[,3]+group)
		}
	}
	return(design)
}

### function to run the differential expression analysis on one dataset
	# includes the computation of true and false positive rates
perform_dea = function(n_sim, group, corr=0, sva=NULL, ruv=NULL, n_sv=NULL, sva_method="be"){
	counts = read.table(paste0("counts", n_sim))
	design = determine_design(counts, group, sva_method, n_sv)
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
compute_accuracy = function(n_simulations, n_replicates, corr=0, sva=NULL, ruv=NULL){
	if (!is.null(sva) && !is.null(ruv)){
		stop("choose either sva or ruv correction")	
	}
	if (length(n_replicates) == 1) {
		n_replicates = rep(n_replicates, 2)
	}
	group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
	aucs = c()
	for (n_sim in seq(n_simulations)){
		results = perform_dea(n_sim, group, corr, sva, ruv)
		aucs = c(aucs, auc(results))
	}
	return(aucs)
}
