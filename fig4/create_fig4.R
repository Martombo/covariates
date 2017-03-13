library(ggplot2)
library(sva)
source("../generate_counts.R")

set.seed(5432)
n_samples = c(seq(3, 12, by=3))
n_covars = c(0,1,2)

est_covar = data.frame(n_samples=rep(n_samples,length(n_covars)))
est_covar$n_covars = rep(n_covars,each=length(n_samples))
est_covar$n_est_covars = as.double(rep(n_covars,each=length(n_samples)))
est_covar$legend = factor("true",levels=c("true","be","leek"))
for (n_sim_covar in n_covars){
	for (k in n_samples){
		n_replicates = c(k, k)
		data = simulate_dataset(n_replicates, 1, n_covariates=n_sim_covar, cov_width=1, cov_strength_max=2)
		group = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
		num_sv_be = num.sv(as.matrix(read.table("counts1")), model.matrix(~group))
		num_sv_leek = num.sv(as.matrix(read.table("counts1")), model.matrix(~group), method="leek")
		est_covar = rbind(est_covar, c(k, n_sim_covar, num_sv_be, "be"))
		est_covar = rbind(est_covar, c(k, n_sim_covar, num_sv_leek, "leek"))
	}
}

save.image()
#load(".RData")

est_covar$n_est_covars = as.double(est_covar$n_est_covars)
est_covar$n_samples = as.double(est_covar$n_samples)
hline_data = data.frame(n_covars=n_covars, n_est_covars=n_covars)
p = ggplot(est_covar, aes(x=n_samples, y=n_est_covars, colour=legend))

svg("fig4.svg",height=6,width=6)
	p + geom_point(size=2) + facet_wrap(~n_covars, ncol=1) + geom_hline(aes(yintercept=n_est_covars), hline_data) + xlab("n samples") + ylab("n estimated covariates") + scale_x_continuous(breaks=n_samples)
dev.off()
