library(ggplot2)
source("../generate_counts.R")
source("../compute_accuracy.R")

n_simulations = 20
n_covariates = c(1,2,3)
set.seed(1234)

simulate_dataset(NA, n_simulations, n_covariates=0)
aucs = compute_all_accuracy(n_simulations, 0)
for (n_covar in n_covariates){
	simulate_dataset(NA, n_simulations, n_covariates=n_covar, depth=NA, n_degs=NA)
	aucs = rbind(aucs, compute_all_accuracy(n_simulations, n_covar))
}
save.image()

p = ggplot(aucs, aes(x=corr, y=auc, fill=corr))
svg("fig5.svg")
	p + geom_boxplot() + facet_wrap(~n_covar, ncol=1)
dev.off()
