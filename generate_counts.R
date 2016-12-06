library(compcodeR)
library(methods)

### set options ### ### ###

n_genes = 10000
batch_strength_min = 0
batch_strength_max = 2
batch_width_min = 0
batch_width_max = 1
n_replicates = 3
min_counts_gene = 0
n_degs = 1000
deg_min_fc = 1.5
depth = 1e+07

### batch simulation function ### ### ###
batch_simulation = function(dat, strength_min=0, strength_max=2, width_min=0, width_max=1, min_counts_gene=0){
	counts = dat@count.matrix
	n_replicates = length(counts[1,]) / 2
	n_genes = length(counts[,1])
	
	for (i in seq(n_replicates)){
		batch_strength = runif(1, strength_min, strength_max)
		n_batch_genes = round(runif(1, width_min, width_max) * n_genes)
		batch_genes = sample(seq(n_genes))[1:n_batch_genes]
		batch_effect = rep(1,n_genes)
		batch_effect[batch_genes] = abs(rnorm(n_batch_genes, sd=batch_strength)) + 1
		half_batch_genes = batch_genes[1:round(n_batch_genes / 2)]
		batch_effect[half_batch_genes] = 1 / batch_effect[half_batch_genes]
		dat@count.matrix[,i] = round(dat@count.matrix[,i]*batch_effect)
		dat@count.matrix[,n_replicates+i] = round(dat@count.matrix[,n_replicates+i]*batch_effect)
	}

	dat@count.matrix = dat@count.matrix[which(rowSums(dat@count.matrix)>=(min_counts_gene)),]
	return(dat)
}


### simulate! ### ### ###

for(k in seq(2)){

# generate data
	dat = generateSyntheticData(dataset="simulated",
		n.vars=n_genes,
		samples.per.cond=n_replicates,
		seqdepth=depth,
		n.diffexp=n_degs,
		fraction.upregulated=0.5,
		between.group.diffdisp=F,
		effect.size=deg_min_fc,
		repl.id=1,
		filter.threshold.total=min_counts_gene)
	write.table(dat@count.matrix, file=paste0("counts",k), col.names=FALSE, quote=FALSE)

# simulate batch
	dat_batch = batch_simulation(dat,
		strength_min=batch_strength_min,
		strength_max=batch_strength_max,
		width_min=batch_width_min,
		width_max=batch_width_max,
		min_counts_gene=min_counts_gene)
	write.table(dat_batch@count.matrix, file=paste0("counts", k), col.names=FALSE, quote=FALSE)
}
