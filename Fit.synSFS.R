# LT 23/06/2021

# analyze SFS of synonymous variants by codon

library(dplyr)
library(tidyr)

gnomad = read.delim(file='extract.synonymous.SFS.bycodon.tsv')
gnomad = subset(gnomad, is.finite(cat))

# remove ansum and anmax
SFS = gnomad[,-(4:5)] %>% spread(key=cat, value=obs, fill=0)

# get AN: group by codon, get sum of an and total number of variants
AN = gnomad %>% group_by(codon) %>% summarize(an=sum(ansum)/sum(obs), anmax=max(anmax))

# a quick look at prop of singletons:
summary((SFS[,-1]/rowSums(SFS[,-1]))[,1])
prop = cbind(SFS[,1],SFS[,-1]/rowSums(SFS[,-1]))
