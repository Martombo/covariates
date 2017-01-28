library(ggplot2)
source("../generate_counts.R")
source("../compute_accuracy.R")

cov_width = 1
strong_effect = 1
weak_effect = 0.5
n_replicates = c(3, 3)
set.seed(6125)

group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))

data1 = simulate_dataset(n_replicates, cov_width=cov_width, data_file="counts0")
counts = read.table("counts1")
res1 = perform_dea(counts, 1, group, model.matrix(~group))
write.table(res1, "res1", quote=F)

cov_strengths = c(strong_effect, strong_effect, 0, weak_effect, 0, 0)
data2 = simulate_dataset(n_replicates, cov_width=cov_width, cov_strengths=cov_strengths, data_file="counts0", n_covariates=1)
counts = read.table("counts1")
res2 = perform_dea(counts, 1, group, model.matrix(~cov_strengths+group))
write.table(res2, "res2", quote=F)

auc_data = data.frame(FPR=res1$fpr, TPR=res1$tpr, dataset=factor("no_covariate", levels=c("covariate","no_covariate")))
auc_data = rbind(auc_data, data.frame(FPR=res2$fpr, TPR=res2$tpr, dataset="covariate"))

p = ggplot(auc_data, aes(x=FPR, y=TPR, colour=dataset))

svg("fig2.svg", height=5, width=5)
	p + geom_point() + xlim(0,0.1)
dev.off()
