library(ggplot2)
source("../modified_compcode.R")
source("../generate_counts.R")

# function to plot PCA of counts1 file
plotPCA1 = function(contrast, plot_file="PCAplot"){
	counts = read.table("counts1", row.names = 1)
	colData = data.frame(contrast = contrast)
	dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~contrast)

	rld = rlogTransformation(dds)
	svg(paste0(plot_file, ".svg"))
		print(plotPCA(rld, intgroup=c("contrast"), ntop=3000))
	dev.off()
}

set.seed(1)
n_replicates = c(3, 3)
contrast = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))

simulate_dataset(n_replicates, cov_strengths=c(1,1,0,0,0,1))
plotPCA1(contrast, "PCA_covar")
