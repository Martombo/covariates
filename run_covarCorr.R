library(ggplot2)
source("generate_counts.R")
source("compute_accuracy.R")

n_simulations = 1
n_replicates = 3

covariates = simulate_dataset(n_replicates, n_simulations, cov_strength_min=0, cov_strength_max=0.5, cov_width_min=1, data_file="counts0")
svs = estimate_sv(n_simulations, n_sv=1)

print(cbind(svs[[1]]$sv, covariates[[1]]$cov_strength1))

counts = read.table("counts1", row.names = 1)
colData = data.frame(contrast = rep(c("A","B"), each=3))
dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~contrast)

rld = rlogTransformation(dds)

svg("PCAplot.svg")
	plotPCA(rld, intgroup=c("contrast"), ntop=1000)
dev.off()

dds = DESeq(dds)

svg("dispPlot.svg")
    plotDispEsts(dds)
dev.off()

save.image()
