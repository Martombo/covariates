simulate_dataset = function(n_replicates, n_simulations, n_genes=10000, n_covariates=1, cov_strength_min=0,
	cov_strength_max=2, cov_width_min=0, cov_width_max=1, min_counts=0, n_degs=3000, deg_min_fc=2, depth=3e+07){

	require(compcodeR)
	require(methods)

# apply covariate function
	apply_covariate = function(counts, covariates, n_sim, n_samples, n_genes, min_counts){

		for (n_cov in seq(length(covariates))){
			cov_strength = covariates[[n_sim]][, paste0("cov_strength", n_cov)]
			n_cov_genes = round(covariates[[n_sim]][, paste0("cov_width", n_cov)] * n_genes)

			for (n_rep in seq(n_samples)){
				cov_genes = sample(seq(n_genes))[seq(0, n_cov_genes[n_rep])]
				cov_effect = pmax(rnorm(n_cov_genes[n_rep], sd=cov_strength[n_rep]), -1)
				cov_var = round(counts[cov_genes, n_rep] * cov_effect)
				counts[cov_genes, n_rep] = counts[cov_genes, n_rep] + cov_var
			}
		}
		counts = counts[which(rowSums(counts)>=(min_counts)),]
		return(counts)
	}


### simulate datasets ###

	covariates = list()
	n_samples = n_replicates*2
	for (n_sim in seq(n_simulations)){

	# generate data with no bias
		dat = generateSyntheticData(dataset="simulated",
			n.vars=n_genes,
			samples.per.cond=n_replicates,
			seqdepth=depth,
			n.diffexp=n_degs,
			fraction.upregulated=0.5,
			between.group.diffdisp=F,
			effect.size=deg_min_fc,
			repl.id=1,
			filter.threshold.total=min_counts)

		# apply random covariates properties
		covariates[[n_sim]] = data.frame(row.names=seq(n_samples))
		for (n_cov in seq(n_covariates)){
			#covariates[[n_sim]][, paste0("cov_strength", n_cov)] = runif(n_samples, cov_strength_min, cov_strength_max)
			covariates[[n_sim]][, "cov_strength1"] = c(0,0,1,0,0,1)
			#covariates[[n_sim]][, paste0("cov_width", n_cov)] = runif(n_samples, cov_width_min, cov_width_max)
			covariates[[n_sim]][, "cov_width1"] = c(0,0,1,0,0,1)
		}
		cov_counts = apply_covariate(dat@count.matrix, covariates, n_sim, n_samples, n_genes, min_counts)
		write.table(cov_counts, file=paste0("counts", n_sim), col.names=F, quote=F)
	}
	return(covariates)
}
