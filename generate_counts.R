### set options ### ### ###

nGenes=10000
reps=3
filterTS=0
eSize=1.5
bgdd=F
nDiff=1000
depth=1e+07

batchStInter=seq(0,0.5,by=0.1)

### simulate! ### ### ###

source("../old/functions.R")
library(compcodeR)

for(k in seq(20)){

# generate data
	dat=generateSyntheticData(dataset="simulated",
		n.vars=nGenes,
		samples.per.cond=reps,
		seqdepth=depth,
		n.diffexp=nDiff,
		fraction.upregulated=0.5,
		between.group.diffdisp=bgdd,
		effect.size=eSize,
		repl.id=1,
		filter.threshold.total=filterTS)
	write.table(dat@count.matrix,file=paste0("counts",k),col.names=FALSE,quote=FALSE)

# simulate batch
	for(batchSt in batchStInter){
		datB=batchEffect(dat,batchSt=batchSt,filterTS=filterTS,tripleEff=TRUE)
		write.table(datB@count.matrix,file=paste0("counts",k,"_",batchSt),col.names=FALSE,quote=FALSE)
	}
}
