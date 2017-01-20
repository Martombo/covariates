library(ggplot2)
source("generate_counts.R")
source("compute_accuracy.R")

n_simulations = 1
n_replicates = c(2,4)

covariates = simulate_dataset(n_replicates, n_simulations, cov_strength_min=0, cov_strength_max=0.5, cov_width_min=1)
svs = estimate_sv(n_simulations, n_sv=1)

print(cbind(svs[[1]]$sv, covariates[[1]]$cov_strength1))

counts = read.table("counts1")
colData = data.frame(contrast = c(rep("A", n_replicates[1]), rep("B", n_replicates[2])))
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~contrast)

rld = rlogTransformation(dds)

svg("PCAplot.svg")
	plotPCA(rld, intgroup=c("contrast"), ntop=1000)
dev.off()

#dds = DESeq(dds)
#
#svg("dispPlot.svg")
#    plotDispEsts(dds)
#dev.off()
#
#save.image()
