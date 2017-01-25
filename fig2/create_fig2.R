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
res1 = perform_dea(1, group)

data2 = simulate_dataset(n_replicates, cov_width=cov_width, cov_strengths=c(strong_effect, strong_effect, 0, weak_effect, 0, 0), data_file="counts0")
res2 = perform_dea(1, group)

auc_data = data.frame(fpr=res1$fpr, tpr=res1$tpr, dataset=factor("no_covar", levels=c("covar","no_covar")))
auc_data = rbind(auc_data, data.frame(fpr=res2$fpr, tpr=res2$tpr, dataset="covar"))

p = ggplot(auc_data, aes(x=fpr, y=tpr, colour=dataset))

svg("fig2.svg")
	p + geom_point() + xlim(0,0.1)
dev.off()
