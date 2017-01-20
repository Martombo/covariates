library(ggplot2)
source("../modified_compcode.R")
source("../generate_counts.R")
source("../compute_accuracy.R")

n_genes = 10000
cov_width = 1
n_covar_genes = n_genes * cov_width
strong_effect = 1
weak_effect = 0.5
n_replicates = c(2, 2)
set.seed(3208)

group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))

data = simulate_dataset(n_replicates, n_genes=n_genes, cov_width=cov_width, data_file="counts0")
res1 = compute_accuracy(1, group)
res2 = compute_accuracy(1, group, corr=T)

auc_data = data.frame(fpr=res1$fpr, tpr=res1$tpr, dataset=factor("not_corr", levels=c("not_corr", "corr")))
auc_data=rbind(auc_data, data.frame(fpr=res2$fpr, tpr=res2$tpr, dataset="corr"))

p = ggplot(auc_data, aes(x=fpr, y=tpr, colour=dataset))

svg("fig3.svg")
	p + geom_point() + xlim(0,0.1)
dev.off()
