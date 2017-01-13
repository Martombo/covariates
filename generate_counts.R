simulate_dataset = function(n_replicates, n_simulations, n_genes=10000, n_degs=3000, deg_min_fc=1.5,  n_covariates=1, 
			cov_strength_min=0, cov_strength_max=1, cov_width_min=0, cov_width_max=1, cov_decreasing_factor=0.8,
			min_counts=0, data_file=NA, depth=3e+07, relmeans="auto", dispersions="auto"){

	require(methods)
	require(DESeq2)
	source("compcodeR_mod.R")

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

# take depth relmeans and dispersions from a data_file, if specified
	if (!is.na(data_file)){
		counts_data = read.table(data_file, row.names = 1)
		counts_data = counts_data[which(rowSums(counts_data) > min_counts), ]
		n_genes = length(counts_data[,1])
		data_contrast = c(rep("A", n_replicates[1]), rep("B", n_replicates[2]))
		col_data = data.frame(contrast = data_contrast)
		dds = DESeqDataSetFromMatrix(countData = counts_data, colData = col_data, design = ~contrast)
		dds = estimateSizeFactors(dds)
		dds = estimateDispersionsGeneEst(dds)
		depth = sum(counts(dds)) / length(counts_data[1,])
		relmeans = mcols(dds)$baseMean
		dispersions = mcols(dds)$dispGeneEst
	}

	n_samples = sum(n_replicates)
	covariates = list()
	cov_correction = 1
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
			filter.threshold.total=min_counts,
			relmeans=relmeans,
			dispersions=dispersions)

	# apply random covariates properties
		covariates[[n_sim]] = data.frame(row.names=seq(n_samples))
		for (n_cov in seq(n_covariates)){
			covariates[[n_sim]][, paste0("cov_strength", n_cov)] = runif(n_samples, cov_strength_min, cov_strength_max) * cov_correction
			covariates[[n_sim]][, paste0("cov_width", n_cov)] = runif(n_samples, cov_width_min, cov_width_max) * cov_correction
			cov_correction = cov_correction * cov_decreasing_factor
		}
		cov_counts = apply_covariate(dat@count.matrix, covariates, n_sim, n_samples, n_genes, min_counts)
		write.table(cov_counts, file=paste0("counts", n_sim), col.names=F, quote=F)
	}
	return(covariates)
}
