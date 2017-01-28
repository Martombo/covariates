simulate_dataset = function(n_replicates=NA, n_simulations=1, n_genes=10000, n_degs=2000, deg_min_fc=1.3,  n_covariates=0, 
			cov_strength_min=0, cov_strength_max=2, cov_strengths=NA, cov_width=1, cov_decreasing_factor=1,
			min_counts=0, data_file=NA, data_reps=NA, depth=3e+07, relmeans="auto", dispersions="auto"){

	require(methods)
	require(DESeq2)
	source("../modified_compcode.R")

	if (is.na(n_replicates)){
		n_replicates = sample(3:12, 1)
	}
	if (is.na(depth)){
		depth = runif(1,1e+07,5e+07)
	}
	if (is.na(n_degs)){
		n_degs = runif(1, 1000, 3000)
	}

# apply covariate function
	apply_covariate = function(counts, cov_strength, cov_genes, n_samples, min_counts){
		cov_effect = rnorm(length(cov_genes), sd=1)
		for (n_sample in seq(n_samples)){
			cov_effect_sample = 2^(cov_effect * cov_strength[n_sample])
			cov_var = round(counts[cov_genes, n_sample] * cov_effect_sample)
			counts[cov_genes, n_sample] = counts[cov_genes, n_sample] + cov_var
		}
		counts = counts[which(rowSums(counts)>=(min_counts)),]
		return(counts)
	}


### simulate datasets ###

# take genes mean and dispersion from a data_file, if specified
	if (!is.na(data_file)){
		if (is.na(data_reps[1])){
			data_reps = n_replicates
		}
		counts_data = read.table(data_file)
		counts_data = counts_data[which(rowSums(counts_data) > min_counts), ]
		n_genes = length(counts_data[,1])
		data_contrast = c(rep("A", data_reps[1]), rep("B", data_reps[2]))
		col_data = data.frame(contrast = data_contrast)
		dds = DESeqDataSetFromMatrix(countData = counts_data, colData = col_data, design = ~contrast)
		dds = estimateSizeFactors(dds)
		dds = estimateDispersionsGeneEst(dds)
		depth = sum(counts(dds)) / length(counts_data[1,])
		relmeans = mcols(dds)$baseMean
		dispersions = mcols(dds)$dispGeneEst
	}

	n_samples = sum(n_replicates)
	strength_list = list()
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
		counts = dat@count.matrix
		true_log2FCs = dat@variable.annotations$truelog2foldchanges

	# apply covariates, random if not defined
		cov_correction = 1
		strength_list[[n_sim]] = data.frame(row.names=seq(n_samples))
		if (n_covariates > 0){
			for (n_cov in seq(n_covariates)){
				if (!is.na(cov_strengths[1])){
					strength_list[[n_sim]][, n_cov] = cov_strengths * cov_correction
				}else{
					strength_list[[n_sim]][, n_cov] = runif(n_samples, cov_strength_min, cov_strength_max) * cov_correction
				}
				cov_genes = sample(seq(n_genes))[1:(cov_width * n_genes)]
				cov_correction = cov_correction * cov_decreasing_factor
				counts = apply_covariate(counts, strength_list[[n_sim]][, n_cov], cov_genes, n_samples, min_counts)
			}
		}
		write.table(counts, paste0("counts", n_sim), row.names=F, col.names=F, quote=F)
		write.table(true_log2FCs, paste0("log2FCs", n_sim), row.names=F, col.names=F, quote=F)
	}
}
