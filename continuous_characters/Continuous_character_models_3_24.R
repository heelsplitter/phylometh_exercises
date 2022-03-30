library(ape) 
library(geiger)
library(OUwie)
library(phytools)

#Now get the tree and data. For these exercises, knowing uncertainty in your measurements can also be important. (remember for homework to change eval=FALSE to eval=TRUE).

tree <- read.tree("C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Continuous_character_models/RAxML_bipartitions.Porotheleaceae_multilocus_Trimal_2022_03_20_edit_3.tre")

continuous.data <- read.csv(file="C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Continuous_character_models/Tree_characters.csv", stringsAsFactors=FALSE) #death to factors.

#A function to clean data, make sure taxon names match between tree and data, etc.

rownames(continuous.data) <- continuous.data[,1]
continuous.data <- continuous.data[-1]

continuous_vector <- continuous.data[,1]
names(continuous_vector) <- rownames(continuous.data)

CleanData <- geiger::treedata(phy=tree, data=continuous_vector, sort=TRUE) 

#print(ape::Ntip(CleanData$phy))

viewdata <- phytools::contMap(chronos(CleanData$phy), CleanData$data[,1])

print(viewdata)

#First, start basic. What is the rate of evolution of your trait on the tree?

BM1 <- geiger::fitContinuous(CleanData$phy, CleanData$data, model="BM")
#print(BM1)
print(paste("The rate of evolution is", 1788.185755 , "in units of", "substitutions per site"))

#Important: What are the rates of evolution? In what units?

##The fitcontinuous output under a Brownian Motion model gives the rate of evolution (sigma^2) as 1788.185755. I am uncertain what units this would be in. My continuous character values are unitless (μm/μm) and I am not certain what the time unit would be for the rate. Before using chronos on my tree, I would have said substitutions per site, but I believe that was altered to make the phylogeny ultrametric.
  
OU1 <- geiger::fitContinuous(chronos(CleanData$phy), CleanData$data, model="OU")
#print(OU1)
par(mfcol = c(1,2))
plot(chronos(CleanData$phy), show.tip.label=FALSE)
ou.tree <- rescale(chronos(CleanData$phy), model="OU", 2.718282)
plot(ou.tree)

#How are the trees different?

##The tree differ in 2 major ways. Firstly, the rescaled OU tree appears to have changed one node to a polytomy that appears to be resolved in the non-rescaled tree ((Henningsomyces+Rectipilus)+(Hemimycena crispula+Phloeomana)+(core Porotheleaceae)). The other way in which these trees differ is in branch lengths. Specifically, in the rescaled tree, the branch lengths to most nodes from the root have been shortened, while some branch lengths from tips to these nodes have been increased.
  
#Compare trees

AIC.BM1 <- 288.664707
AIC.OU1 <- 241.785885

##The "best" model here is the one with the lowest AIC score, in this case the Ornstein–Uhlenbeck model OU1.

delta.AIC.BM1 <- AIC.BM1 - AIC.OU1
delta.AIC.OU1 <- AIC.OU1 - AIC.OU1

##OUwie runs##

#This takes longer than you may be used to.

#We’re a bit obsessive about doing multiple starts and in general performing a thorough numerical search. It took you 3+ years to get the data, may as well take an extra five minutes to get an accurate answer

#First, we need to assign regimes. The way we do this is with ancestral state estimation of a discrete trait. We can do this using ace() in ape, or similar functions in corHMM or diversitree. Use only one discrete char.

one.discrete.char <- read.csv(file="C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Continuous_character_models/Tree_characters_discrete.csv", stringsAsFactors=FALSE)

##MAKE VECTOR!!!

rownames(one.discrete.char) <- one.discrete.char[,1]
one.discrete.char <- one.discrete.char[-1]
discrete_vector <- one.discrete.char[,1]
names(discrete_vector) <- rownames(one.discrete.char)

CleanData2 <- geiger::treedata(phy=tree, data=discrete_vector, sort=TRUE) 

#print(ape::Ntip(CleanData2$phy))

reconstruction.info <- ace(CleanData2$data, chronos(CleanData2$phy), type="discrete", method="ML", CI=TRUE)
best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]

#Now add these labels to your tree.

tree_with_labels <- ape::chronos(CleanData$phy) # fast and dirty way to make ultrametric

tree_with_labels$node.label <- "1"

############


tree_with_labels$node.label <- as.character(best.states) # storing the best reconstruction at each internal node on the tree

tree_with_simmap <- phytools::make.simmap(ape::chronos(CleanData$phy), discrete_vector)

# make a data set for the species at the tips: name, regime (placeholder), continuous trait
ouwie_data <- data.frame(Species=rownames(CleanData$data), Reg=NA, Trait=CleanData$data)

for (ouwie_row in sequence(nrow(ouwie_data))) {
  ouwie_data$Reg[ouwie_row] <- discrete_vector[ouwie_data$Species[ouwie_row]] # for the ith species on the ouwie data, look up the same species name in the discrete tip data and assign that as the "regime"
}

oumv <- OUwie(tree_with_labels, ouwie_data, model="OUMV")

#Error in is.nloptr(ret) : at least one element in x0 > ub

##This may be an issue with polytomies in my tree, which OUwie apparently can't handle

print(oumv)
ou1 <- OUwie(tree_with_labels, ouwie_data, model="OU1")

bm1 <- OUwie(tree_with_labels, ouwie_data, model="BM1")

bms <- OUwie(tree_with_labels, ouwie_data, model="BMS")

oumva <- OUwie(tree_with_labels, ouwie_data, model="OUMVA")

oumva_simmap <- OUwie(tree_with_simmap, ouwie_data, model="OUMVA", simmap.tree=TRUE)

print(c(bm1$AICc, ou1$AICc, oumv$AICc, bms$AICc, oumva$AICc))

