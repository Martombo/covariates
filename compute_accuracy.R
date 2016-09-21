library(limma)
library(edgeR)

counts_dir="~/Dropbox/projects/microprojects/covariates/counts/"
simuls = 2#seq(20)
batches = seq(0,1,by=0.2)
group = rep(c("A","B"),each=3)
group_b = rep(c(1,2,3),2)
design = model.matrix(~group)
design_b1 = model.matrix(~group+group_b)
n_DEGs = 1000
n_top = 1000


voom_res=function(counts_file,group,design,n_top){
		counts=as.matrix(read.table(counts_file,row.names=1))
		dge=DGEList(counts=counts,group=group)
		dge=calcNormFactors(dge)
		v=voom(dge,design)
		fit=lmFit(v,design)
		fit=eBayes(fit)
		return(topTable(fit,n=n_top,coef=2))
}
aucc=function(res,n_DEGs,n_top){
	res$gene_n=as.numeric(gsub("g","",row.names(res)))
	tpr=0; auc=0; fp=0
	for(n in seq(length(res[,1]))){
		if(res$gene_n[n]<n_DEGs){
			tpr=tpr+1/n_DEGs
		}else{
			auc=auc+tpr; fp=fp+1
		}
	}
	auc=auc+tpr; fp=fp+1
	return(auc/fp)
}


for(n_simuls in simuls){
	for(n_batch in batches){
		counts_file = paste0(counts_dir,"counts",n_simuls,"_",n_batch) 
		top_res = voom_res(counts_file,group,design,n_top)
		auc = aucc(top_res,n_DEGs,n_top)
		top_res_b = voom_res(counts_file,group,design_b1,n_top)
		auc_b = aucc(top_res_b,n_DEGs,n_top)
		print(c(auc,auc_b))
	}
}
