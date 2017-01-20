library(ggplot2)
source("../modified_compcode.R")
source("../generate_counts.R")

n_genes = 10000
cov_width = 1
n_covar_genes = n_genes * cov_width
strong_effect = 1
weak_effect = 0.5
n_replicates = c(3, 3)
set.seed(3521)

# function to plot PCA of counts1 file
plotPCA1 = function(contrast, plot_file="PCAplot"){
	counts = read.table("counts1")
	colData = data.frame(contrast = contrast)
	dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~contrast)

	rld = rlogTransformation(dds)
	svg(paste0(plot_file, ".svg"))
		print(plotPCA(rld, intgroup=c("contrast"), ntop=3000))
	dev.off()
}

contrast = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))

simulate_dataset(n_replicates, n_genes=n_genes, cov_width=cov_width, data_file="counts0")
plotPCA1(contrast, "fig1_A")

simulate_dataset(n_replicates, n_genes=n_genes, cov_width=cov_width, cov_strengths=c(strong_effect, strong_effect, 0, weak_effect, 0, 0), data_file="counts0")
plotPCA1(contrast, "fig1_B")

cov_effect = data.frame(covar_log2FC=c(rnorm(n_covar_genes)*strong_effect, rep(0, (n_genes-n_covar_genes)*2), rnorm(n_covar_genes)*weak_effect))
cov_effect$legend = factor(rep(c(strong_effect, weak_effect), each=n_genes))
p = ggplot(cov_effect, aes(x=legend, y=covar_log2FC, fill=legend))
svg("fig1_C.svg")
	p + geom_boxplot()
dev.off()
