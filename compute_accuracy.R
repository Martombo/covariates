library(limma)
library(edgeR)
library(DESeq2)
library(sva)

n_simuls = 2
n_replicates = 3
group = rep(c("A","B"),each=n_replicates)
group_b = rep(seq(n_replicates),2)
design = model.matrix(~group)
n_DEGs = 1000
n_top = 1000


voom_res = function(counts_file, group, design, n_top){
		counts = as.matrix(read.table(counts_file, row.names=1))
		dge = DGEList(counts=counts, group=group)
		dge = calcNormFactors(dge)
		v = voom(dge, design)
		fit = lmFit(v, design)
		fit = eBayes(fit)
		return(topTable(fit, n=n_top, coef=2))
}

aucc = function(res, n_DEGs, n_top){
	res$gene_n = as.numeric(gsub("g", "", row.names(res)))
	tpr=0; auc=0; fp=0
	for(n in seq(length(res[,1]))){
		if(res$gene_n[n] < n_DEGs){
			tpr = tpr + 1 / n_DEGs
		}else{
			auc=auc+tpr; fp=fp+1
		}
	}
	auc=auc+tpr; fp=fp+1
	return(auc / fp)
}

sva_corr = function(counts, group){
	dds = DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group=group), design=~group)
	dds = estimateSizeFactors(dds)
	dat = counts(dds, normalized=TRUE)
	dat = dat[rowMeans(dat)>5,]
	mod = model.matrix(~group, colData(dds))
	mod0 = model.matrix(~1, colData(dds))
	svseq = svaseq(dat, mod, mod0)
	return(svseq)
}

for (simul in seq(n_simuls)){
	counts_file = paste0("counts", simul) 
	counts = read.table(counts_file, row.names=1)
	n_sv = sva_corr(counts, group)$n.sv
	print(n_sv)
}

#	top_res = voom_res(counts_file, group, design, n_top)
#	auc = aucc(top_res, n_DEGs, n_top)
#	top_res_b = voom_res(counts_file, group, design_b1, n_top)
#	auc_b = aucc(top_res_b, n_DEGs, n_top)
#	print(c(auc, auc_b))
