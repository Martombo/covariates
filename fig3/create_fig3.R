library(ggplot2)
source("../generate_counts.R")
source("../compute_accuracy.R")

n_simulations = 20
n_replicates = c(3, 3)
cov_width = 0.2
cov_strengths = c(0.4, 0.2, 0, 0.2, 0, 0)
set.seed(91208)

results_auc = data.frame(AUC=rep(1,n_simulations*2), covariate=rep(c("weak", "strong"), each=n_simulations), correction=factor("not_corrected", levels=c("not_corrected", "corrected")))
data = simulate_dataset(n_replicates, n_simulations, cov_strengths=cov_strengths, cov_width=cov_width, n_covariates=1, data_file="counts0")
res = compute_accuracy(n_simulations, n_replicates)
res_corr = compute_accuracy(n_simulations, n_replicates, corr=cov_strengths)
results_auc = rbind(results_auc, data.frame(AUC=res_corr/res, covariate="weak", correction="corr"))

cov_width = 1
set.seed(91208)

cov_strengths = c(1, 1, 0, 0.5, 0, 0)
data = simulate_dataset(n_replicates, n_simulations, cov_strengths=cov_strengths, cov_width=cov_width, n_covariates=1, data_file="counts0")
res = compute_accuracy(n_simulations, n_replicates)
res_corr = compute_accuracy(n_simulations, n_replicates, corr=cov_strengths)
results_auc = rbind(results_auc, data.frame(AUC=res_corr/res, covariate="strong", correction="corr"))

results_auc$covariate = factor(results_auc$covariate, levels=c("weak", "strong"))
p = ggplot(results_auc, aes(x=covariate, y=AUC, fill=correction))
svg("fig3.svg", height=5, width=5)
	p + geom_boxplot()
dev.off()
