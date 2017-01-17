estimate_sv = function(n_simulations, n_replicates=3, n_DEGs=3000, n_top=1000, n_sv=0){ 

	require(limma)
	require(edgeR)
	require(DESeq2)
	require(sva)


# voom DEA function
	voom_res = function(counts_file, group, design, n_top){
		counts = as.matrix(read.table(counts_file, row.names=1))
		dge = DGEList(counts=counts, group=group)
		dge = calcNormFactors(dge)
		v = voom(dge, design)
		fit = lmFit(v, design)
		fit = eBayes(fit)
		return(topTable(fit, n=n_top, coef=2))
	}

# AUC computation function
	aucc = function(res, n_DEGs){
		res$gene_n = as.numeric(gsub("g", "", row.names(res)))
		n_top = length(res[,1])
		top_pos = min(n_DEGs, n_top)
		tpr=0; auc=0; fp=0
		for(n in seq(n_top)){
			if(res$gene_n[n] < n_DEGs){
				tpr = tpr + 1 / top_pos
				next
			}
			auc=auc+tpr; fp=fp+1
		}
		auc=auc+tpr; fp=fp+1
		return(auc / fp)
	}

# true positive percentage function
	tpr = function(res, n_DEGs){
		res$gene_n = as.numeric(gsub("g", "", row.names(res)))
		tp = sum(res$gene_n <= n_DEGs)
		return(tp / length(res$gene_n))
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


### run sva estimation ###
	svs = list()
	if (length(n_replicates) == 1) {
		n_replicates = rep(n_replicates, 2)
	}   
	group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
	for (n_simul in seq(n_simulations)){
		counts_file = paste0("counts", n_simul) 
		counts = read.table(counts_file, row.names=1)
		svs[[n_simul]] = sva_corr(counts, group, n_sv=n_sv)
		# top_res = voom_res(counts_file, group, design, n_top)
		# auc = aucc(top_res, n_DEGs, n_top)
	}
	return(svs)
}
