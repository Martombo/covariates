generateSyntheticData = function (dataset, n.vars, samples.per.cond, n.diffexp, repl.id = 1, 
    seqdepth = 1e+07, minfact = 0.7, maxfact = 1.4, relmeans = "auto", 
    dispersions = "auto", fraction.upregulated = 1, between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, random.outlier.high.prob = 0, 
    random.outlier.low.prob = 0, single.outlier.high.prob = 0, 
    single.outlier.low.prob = 0, effect.size = 1.5, output.file = NULL) 
{
	require("compcodeR")
	require("edgeR")

    if (!is.null(output.file)) {
        if (!(substr(output.file, nchar(output.file) - 3, nchar(output.file)) == 
            ".rds")) {
            stop("output.file must be an .rds file.")
        }
    }
    uID <- paste(sample(c(0:9, letters, LETTERS), 10, replace = TRUE), 
        collapse = "")

########## modified part ##########

	if (length(samples.per.cond) == 1) {
		samples.per.cond = rep(samples.per.cond, 2)
	}
	condition <- c(rep(1, samples.per.cond[1]), rep(2, samples.per.cond[2]))
	n.samples <- sum(samples.per.cond)

##########   #########   ##########

    S1 <- which(condition == 1)
    S2 <- which(condition == 2)
    if (length(effect.size) == 1) {
        n.upregulated <- floor(fraction.upregulated * n.diffexp)
        if (fraction.upregulated != 0 & n.diffexp != 0) {
            genes.upreg <- 1:n.upregulated
        }
        else {
            genes.upreg <- NULL
        }
        if (fraction.upregulated != 1 & n.diffexp != 0) {
            genes.downreg <- (n.upregulated + 1):n.diffexp
        }
        else {
            genes.downreg <- NULL
        }
        genes.nonreg <- setdiff(1:n.vars, union(genes.upreg, 
            genes.downreg))
    }
    else {
        if (length(effect.size) != n.vars) {
            stop("The length of the effect.size vector must be the same as the number of simulated genes.")
        }
        else {
            genes.upreg <- which(effect.size > 1)
            genes.downreg <- which(effect.size < 1)
            genes.nonreg <- which(effect.size == 1)
            n.upregulated <- length(genes.upreg)
            n.diffexp <- length(genes.upreg) + length(genes.downreg)
            fraction.upregulated <- n.upregulated/n.diffexp
        }
    }
    differential.expression <- rep(0, n.vars)
    differential.expression[genes.upreg] <- 1
    differential.expression[genes.downreg] <- 1
    upregulation <- rep(0, n.vars)
    upregulation[genes.upreg] <- 1
    downregulation <- rep(0, n.vars)
    downregulation[genes.downreg] <- 1
    if (is.character(relmeans) | is.character(dispersions)) {
        mu.phi.estimates <- system.file("extdata", "Pickrell.Cheung.Mu.Phi.Estimates.rds", 
            package = "compcodeR")
        mu.phi.estimates <- readRDS(mu.phi.estimates)
        mu.estimates <- mu.phi.estimates$pickrell.cheung.mu
        phi.estimates <- mu.phi.estimates$pickrell.cheung.phi
        to.include <- sample(1:length(mu.estimates), n.vars, 
            replace = ifelse(n.vars > length(mu.estimates), TRUE, 
                FALSE))
        truedispersions.S1 <- phi.estimates[to.include]
        truemeans.S1 <- mu.estimates[to.include]
    }
    if (!is.character(relmeans)) {
        if (length(relmeans) != n.vars) 
            stop("The length of the relmeans vector must be the same as the number of simulated genes.")
        truemeans.S1 <- c(relmeans)
    }
    if (!is.character(dispersions)) {
        if (nrow(cbind(dispersions)) != n.vars) 
            stop("The number of provided dispersions must be the same as the number of simulated genes.")
        truedispersions.S1 <- cbind(dispersions)[, 1]
        if (ncol(cbind(dispersions)) > 1) {
            truedispersions.S2 <- cbind(dispersions)[, 2]
        }
        else {
            truedispersions.S2 <- truedispersions.S1
        }
    }
    nfacts <- runif(n.samples , min = minfact, max = maxfact)
    seq.depths <- nfacts * seqdepth
    overdispersed <- rep(1, n.vars)
    if (fraction.non.overdispersed > 0) {
        overdispersed[genes.upreg[1:round(fraction.non.overdispersed * 
            length(genes.upreg))]] <- 0
        overdispersed[genes.downreg[1:round(fraction.non.overdispersed * 
            length(genes.downreg))]] <- 0
        overdispersed[genes.nonreg[1:round(fraction.non.overdispersed * 
            length(genes.nonreg))]] <- 0
    }
    prob.S1 <- truemeans.S1
    prob.S2 <- rep(0, length(prob.S1))
    if (length(effect.size) == 1) {
        for (i in 1:n.vars) {
            if (i %in% genes.upreg) {
                prob.S2[i] <- (effect.size + rexp(1, rate = 1)) * 
                  prob.S1[i]
            }
            else {
                if (i %in% genes.downreg) {
                  prob.S2[i] <- 1/(effect.size + rexp(1, rate = 1)) * 
                    prob.S1[i]
                }
                else {
                  prob.S2[i] <- prob.S1[i]
                }
            }
        }
    }
    else {
        prob.S2 <- c(effect.size) * prob.S1
    }
    true.log2foldchange <- log2(prob.S2/prob.S1)
    sum.S1 <- sum(prob.S1)
    sum.S2 <- sum(prob.S2)
    if (is.character(dispersions)) {
        truedispersions.S2 <- truedispersions.S1
        if (between.group.diffdisp == TRUE) {
            for (i in 1:length(truedispersions.S2)) {
                sample.base <- phi.estimates[abs(log10(mu.estimates) - 
                  log10(prob.S2[i])) < 0.05]
                if (length(sample.base) < 50) {
                  sample.base <- phi.estimates[order(abs(log10(mu.estimates) - 
                    log10(prob.S2[i])))][1:500]
                }
                truedispersions.S2[i] <- sample(sample.base, 
                  1)
            }
        }
    }
    truedispersions.S1 <- truedispersions.S1 * overdispersed
    truedispersions.S2 <- truedispersions.S2 * overdispersed
    Z <- matrix(0, n.vars, length(S1) + length(S2))
    for (i in 1:n.vars) {
        for (j in 1:ncol(Z)) {
            if (j %in% S1) {
                if (overdispersed[i] == 1) {
                  Z[i, j] <- rnbinom(n = 1, mu = prob.S1[i]/sum.S1 * 
                    seq.depths[j], size = 1/truedispersions.S1[i])
                }
                else {
                  Z[i, j] <- rpois(n = 1, lambda = prob.S1[i]/sum.S1 * 
                    seq.depths[j])
                }
            }
            else {
                if (overdispersed[i] == 1) {
                  Z[i, j] <- rnbinom(n = 1, mu = prob.S2[i]/sum.S2 * 
                    seq.depths[j], size = 1/truedispersions.S2[i])
                }
                else {
                  Z[i, j] <- rpois(n = 1, lambda = prob.S2[i]/sum.S2 * 
                    seq.depths[j])
                }
            }
        }
    }
    random.outliers <- matrix(0, nrow(Z), ncol(Z))
    random.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
    if (random.outlier.high.prob != 0 | random.outlier.low.prob != 
        0) {
        for (i in 1:nrow(Z)) {
            for (j in 1:ncol(Z)) {
                tmp <- runif(1)
                if (tmp < random.outlier.high.prob) {
                  random.outliers[i, j] <- 1
                  random.outliers.factor[i, j] <- runif(1, min = 5, 
                    max = 10)
                }
                else if (tmp < random.outlier.low.prob + random.outlier.high.prob) {
                  random.outliers[i, j] <- (-1)
                  random.outliers.factor[i, j] <- 1/runif(1, 
                    min = 5, max = 10)
                }
            }
        }
        Z <- round(random.outliers.factor * Z)
    }
    has.single.outlier <- rep(0, n.vars)
    single.outliers <- matrix(0, nrow(Z), ncol(Z))
    single.outliers.factor <- matrix(1, nrow(Z), ncol(Z))
    if (single.outlier.high.prob != 0 | single.outlier.low.prob != 
        0) {
        has.single.outlier[genes.upreg[1:floor((single.outlier.high.prob + 
            single.outlier.low.prob) * length(genes.upreg))]] <- 1
        has.single.outlier[genes.downreg[1:floor((single.outlier.high.prob + 
            single.outlier.low.prob) * length(genes.downreg))]] <- 1
        has.single.outlier[genes.nonreg[1:floor((single.outlier.high.prob + 
            single.outlier.low.prob) * length(genes.nonreg))]] <- 1
        for (i in 1:nrow(Z)) {
            if (has.single.outlier[i] == 1) {
                the.sample <- sample(1:(ncol(Z)), 1)
                if (runif(1) < (single.outlier.high.prob/(single.outlier.high.prob + 
                  single.outlier.low.prob))) {
                  single.outliers[i, the.sample] <- 1
                  single.outliers.factor[i, the.sample] <- runif(1, 
                    min = 5, max = 10)
                }
                else {
                  single.outliers[i, the.sample] <- (-1)
                  single.outliers.factor[i, the.sample] <- 1/runif(1, 
                    min = 5, max = 10)
                }
            }
        }
        Z <- round(single.outliers.factor * Z)
    }
    rownames(Z) <- 1:n.vars
    n.random.outliers.up.S1 <- apply(random.outliers[, S1] > 
        0, 1, sum)
    n.random.outliers.up.S2 <- apply(random.outliers[, S2] > 
        0, 1, sum)
    n.random.outliers.down.S1 <- apply(random.outliers[, S1] < 
        0, 1, sum)
    n.random.outliers.down.S2 <- apply(random.outliers[, S2] < 
        0, 1, sum)
    n.single.outliers.up.S1 <- apply(single.outliers[, S1] > 
        0, 1, sum)
    n.single.outliers.up.S2 <- apply(single.outliers[, S2] > 
        0, 1, sum)
    n.single.outliers.down.S1 <- apply(single.outliers[, S1] < 
        0, 1, sum)
    n.single.outliers.down.S2 <- apply(single.outliers[, S2] < 
        0, 1, sum)
    nf <- calcNormFactors(Z)
    norm.factors <- nf * colSums(Z)
    common.libsize <- exp(mean(log(colSums(Z))))
    pseudocounts <- sweep(Z + 0.5, 2, norm.factors, "/") * common.libsize
    log2.pseudocounts <- log2(pseudocounts)
    M.value <- apply(log2.pseudocounts[, S2], 1, mean) - apply(log2.pseudocounts[, 
        S1], 1, mean)
    A.value <- 0.5 * (apply(log2.pseudocounts[, S2], 1, mean) + 
        apply(log2.pseudocounts[, S1], 1, mean))
    variable.annotations <- data.frame(truedispersions.S1 = truedispersions.S1, 
        truedispersions.S2 = truedispersions.S2, truemeans.S1 = prob.S1, 
        truemeans.S2 = prob.S2, n.random.outliers.up.S1 = n.random.outliers.up.S1, 
        n.random.outliers.up.S2 = n.random.outliers.up.S2, n.random.outliers.down.S1 = n.random.outliers.down.S1, 
        n.random.outliers.down.S2 = n.random.outliers.down.S2, 
        n.single.outliers.up.S1 = n.single.outliers.up.S1, n.single.outliers.up.S2 = n.single.outliers.up.S2, 
        n.single.outliers.down.S1 = n.single.outliers.down.S1, 
        n.single.outliers.down.S2 = n.single.outliers.down.S2, 
        M.value = M.value, A.value = A.value, truelog2foldchanges = true.log2foldchange, 
        upregulation = upregulation, downregulation = downregulation, 
        differential.expression = differential.expression)
    rownames(variable.annotations) <- rownames(Z)
    sample.annotations <- data.frame(condition = condition, depth.factor = nfacts)
    info.parameters <- list(n.diffexp = n.diffexp, fraction.upregulated = fraction.upregulated, 
        between.group.diffdisp = between.group.diffdisp, filter.threshold.total = filter.threshold.total, 
        filter.threshold.mediancpm = filter.threshold.mediancpm, 
        fraction.non.overdispersed = fraction.non.overdispersed, 
        random.outlier.high.prob = random.outlier.high.prob, 
        random.outlier.low.prob = random.outlier.low.prob, single.outlier.high.prob = single.outlier.high.prob, 
        single.outlier.low.prob = single.outlier.low.prob, effect.size = effect.size, 
        samples.per.cond = samples.per.cond, repl.id = repl.id, 
        dataset = dataset, uID = uID, seqdepth = seqdepth, minfact = minfact, 
        maxfact = maxfact)
    s <- apply(Z, 1, sum)
    keep.T <- which(s >= filter.threshold.total)
    Z.T <- Z[keep.T, ]
    variable.annotations.T <- variable.annotations[keep.T, ]
    filtering <- paste("total count >=", filter.threshold.total)
    cpm <- sweep(Z.T, 2, apply(Z.T, 2, sum), "/") * 1e+06
    m <- apply(cpm, 1, median)
    keep.C <- which(m >= filter.threshold.mediancpm)
    Z.TC <- Z.T[keep.C, ]
    variable.annotations.TC <- variable.annotations.T[keep.C, 
        ]
    filtering <- paste(filtering, "; ", paste("median cpm >=", 
        filter.threshold.mediancpm))
    rownames(Z.TC) <- paste("g", 1:nrow(Z.TC), sep = "")
    colnames(Z.TC) <- paste("sample", 1:ncol(Z.TC), sep = "")
    rownames(sample.annotations) <- colnames(Z.TC)
    rownames(variable.annotations.TC) <- rownames(Z.TC)
    data.object <- compData(count.matrix = Z.TC, variable.annotations = variable.annotations.TC, 
        sample.annotations = sample.annotations, filtering = filtering, 
        info.parameters = info.parameters)
    if (!is.null(output.file)) {
        saveRDS(data.object, file = output.file)
    }
    return(invisible(data.object))
}
