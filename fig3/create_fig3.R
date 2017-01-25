library(ggplot2)
source("../generate_counts.R")
source("../compute_accuracy.R")

n_simulations = 20
n_replicates = c(3, 3)
cov_width = 0.2
covariate = c(0.4, 0.2, 0, 0.2, 0, 0)
set.seed(91208)

results_auc = data.frame(auc=rep(1,n_simulations*2), covar=rep(c("weak", "strong"), each=n_simulations), corr=factor("no_corr", levels=c("no_corr", "corr")))
data = simulate_dataset(n_replicates, n_simulations, cov_strengths=covariate, cov_width=cov_width)
res = compute_accuracy(n_simulations, n_replicates)
res_corr = compute_accuracy(n_simulations, n_replicates, corr=covariate)
results_auc = rbind(results_auc, data.frame(auc=res_corr/res, covar="weak", corr="corr"))

cov_width = 1
covariate = c(2, 1, 0, 2, 1, 0)
set.seed(91208)

data = simulate_dataset(n_replicates, n_simulations, cov_strengths=covariate, cov_width=cov_width)
res = compute_accuracy(n_simulations, n_replicates)
res_corr = compute_accuracy(n_simulations, n_replicates, corr=covariate)
results_auc = rbind(results_auc, data.frame(auc=res_corr/res, covar="strong", corr="corr"))

results_auc$covar = factor(results_auc$covar, levels=c("weak", "strong"))
p = ggplot(results_auc, aes(x=covar, y=auc, fill=corr))
svg("fig3.svg")
	p + geom_boxplot()
dev.off()
