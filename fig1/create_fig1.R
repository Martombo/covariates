library(ggplot2)
source("../generate_counts.R")

n_genes = 10000
cov_width = 1
n_covar_genes = n_genes * cov_width
strong_effect = 1
weak_effect = 0.5
n_replicates = c(4, 4)
set.seed(3521)

# function to plot PCA of counts1 file
plotPCA1 = function(contrast, plot_file="PCAplot"){
	counts = read.table("counts1")
	colData = data.frame(contrast = contrast)
	dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~contrast)

	rld = rlogTransformation(dds)
	svg(paste0(plot_file, ".svg"), height=3, width=5)
		print(plotPCA(rld, intgroup=c("contrast"), ntop=10000) + ylim(-25,25))
	dev.off()
}

contrast = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))

simulate_dataset(n_replicates, data_file="counts0")
plotPCA1(contrast, "fig1_A")

cov_strengths = c(strong_effect, strong_effect, weak_effect, 0, strong_effect, 0, 0, 0)
simulate_dataset(n_replicates, cov_width=cov_width, cov_strengths=cov_strengths, data_file="counts0", n_covariates=1)
plotPCA1(contrast, "fig1_B")

cov_effect = data.frame(covar_log2FC=c(rnorm(n_covar_genes)*strong_effect, rep(0, (n_genes-n_covar_genes)*2), rnorm(n_covar_genes)*weak_effect))
cov_effect$factor = factor(rep(c(strong_effect, weak_effect), each=n_genes))
p = ggplot(cov_effect, aes(x=factor, y=covar_log2FC, fill=factor))
svg("fig1_C.svg", height=4, width=2)
	p + geom_boxplot() + guides(fill=F) + xlab("covariate factor") + ylab("covariate log2FC")
dev.off()
