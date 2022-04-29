#NOTE: My comments / answers will begin with "##"

library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
#devtools::install_github("thej022214/hisse")
##Installed from CRAN instead
library(hisse)

#Let’s initially look just at diversification alone. Simulate a 30 taxon tree with only speciation, no extinction:

my.tree <- TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]]

#As always, plot it:

plot(my.tree)
##Pretty messy, but it's a tree for sure.

#One way to look at trees, and actually what many methods reduce to under the hood, is lineage through time plots.

ape::ltt.plot(my.tree)

#You should see it increasing exponentially. Let’s put it on a log scale:

ape::ltt.plot(my.tree, log="y")

#We can look at multiple trees:
  
yule.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, mu=0, complete=FALSE)
#stop("How to do a multiple ltt plot?")
  
ape::mltt.plot(yule.trees, log="y")
##works nicely

#We can also look at trees with birth and death.

bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.9, complete=FALSE)
ape::mltt.plot(bd.trees, log="y", legend=FALSE)

#And compare them:

depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)

#And zooming in on the final part of the plot:

depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range)) #ape rescales depths
plot(x=c(0, -5), y=c(200, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)

#So even though the the net diversification rate is the same, there are very different patterns: in theory, one can estimate both birth and death rates from these trees. In practice, of course, with rates that change over time due to mass extinctions or trait evolution, missing taxa, etc. it can practically be hard to tell these apart.

#TODO Try plotting some with other diversification parameters. What happens if speciation rate is much higher than extinction rate? How does the simulation change with different values, but keeping their difference constant? If their sum is constant? [remember to change to eval=TRUE to run]

my.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=stop(0.5), mu=stop(2), complete=TRUE)

?TreeSim::sim.bd.taxa

##I kept receiving an error of the format "Error in rexp(1, (length(leaves) * (lambda + mu))) : 0.5" apparently regardless of what values I set for lambda and mu, and whether Complete was set to TRUE or FALSE. Could not figure out how to proceed.

ape::mltt.plot(my.trees, log="y", legend=FALSE)

