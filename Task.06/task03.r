install.packages("learnPopGen")
library(learnPopGen)
coalescent.plot(n=50, ngen=100, col.order="alternating")
?coalescent.plot()
#1. Ten alleles, change the population of the simulation
#2. 4.6 generations
#3. 1, 
#4. The most fit allele makes it through the generations while the others go to fixation
#5. yes





library(coala) 
library(phytools)
model <- coal_model(sample_size = 5, loci_number = 10, loci_length = 500, ploidy = 2) + 
feat_mutation(10) +
feat_recombination(10) +
sumstat_trees() +
sumstat_nucleotide_div()
stats <- simulate(model, nsim =1)
Diversity <- stats$pi
Diversity
# No they are not,

Nloci <- length(stats$trees)
t1 <- read.tree(text=stats$trees[[1]][1])
plot(t1)
axisPhylo()
# Because of the 

Age1 <- max(nodeHeights(t1))
t2 <- read.tree(text=stats$trees[[2]][1])
plot(t2)
axisPhylo()
# It goes to 0.3, and it is not the same age

par(mfrow=c(1,2))
plot(t1)
axisPhylo()
plot(t2)
axisPhylo()
compare.chronograms(t1, t2)
# The do not match

tl_1 <- read.tree(text=stats$trees[[1]][1])
tl_2 <- read.tree(text=stats$trees[[1]][2])
compare.chronograms(tl_1, tl_2)
for(locus in 1:Nloci) {
ntree <- length(stats$trees[[locus]])
for(n in 1:ntree){
if(locus == 1&& n ==1){
outPhy <- read.tree(text=stats$trees[[locus]][n])
}
else{
outPhy <- read.tree(text=stats$trees[[locus]][n])
}
}
}
par(mfrow=c(1, 1))
densityTree(outPhy)
model3 <- coal_model(10,50) +
  feat_mutation(par_prior("theta", sample.int(100,1))) +
  sumstat_nucleotide_div()
stats <- simulate(model3, nsim = 40)
mean_pi <- sapply(stats, function(x) mean(x$pi))
theta <- sapply(stats, function(x) x$pars[["theta"]])

mean_pi
theta
plot(x=mean_pi, y=theta)
?lm
ModelR <- lm(theta ~ mean_pi)
summary(ModelR)
abline(ModelR, col='blue')