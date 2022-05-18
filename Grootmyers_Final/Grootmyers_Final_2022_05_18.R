#INTRODUCTION:
##Meiospore ("spore" for simplicity) dimensions are often a very important feature in differentiating fungal species morphologically, and are also occasionally used to differentiate genera. 
##As such, I would like to used continuous trait analyses to investigate whether spore dimensions - average length [L], average width [W] and average length/width [Q]) - vary between supra-species level clades in my clade of interest, the mushroom family Porotheleaceae. 
##I would also like to investigate whether L, W and Q exhibit the same pattern of variation over the phylogeny or whether they differ, and whether Brownian Motion or Ornstein-Uhlenbeck models are better for analyzing spore dimensions as a whole or whether one model is better for a certain spore dimension and not others. 
##I will investigate these issues using a dataset I constructed, which includes a 3 DNA locus maximum likelihood tree generated using Raxml, and a spreadsheet containing spore dimensions for the tips (collections) of the tree.

#LOAD PACKAGES:
library(ape) 
library(geiger)
library(phytools)

#LOADING IN THE TREE:
##This is a RaxML ML tree for 42 individual collections I sequenced using a 3 DNA locus dataset. 
##The model used was GTR + I + G and it was run with 1000 bootstrap replicates.

tree <- read.tree("C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Final/RAxML_bipartitions.Porotheleaceae_multilocus_Trimal_2022_05_05_EDIT2.tre")

#LOADING IN THE MORPHOLOGICAL DATA: 
##This includes average spore length (L) for each tip (collection) and average L for each species assigned to each tip belonging to that species, average spore width (W) for each tip and average W for each species assigned to each tip belonging to that species, average spore quotient (Q) for each tip and average Q for each species assigned to each tip belonging to that species, and lastly the binary trait amyloidity which is the presence/absence of a blue reaction of spores in Melzer's reagent.
##All trait measurements are in units of µm except for Q and amyloidity, which are unit-less.

spore.data <- read.csv(file="C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Final/Tree_characters.csv", stringsAsFactors=FALSE)
rownames(spore.data) <- spore.data[,1]
spore.data <- spore.data[-1]

#CLEAN DATASET AND MATCH IT TO TREE:

clean <- treedata(phy=tree, data=spore.data, sort=TRUE) 

#MAKE OBJECT FOR CLEANED TREE:

tree2 <- clean$phy

#MAKE OBJECTS FOR SPORE TRAITS:

L <- clean$data[,4]
W <- clean$data[,6]
Q <- clean$data[,2]
Amyloid <- clean$data[,7]

#MAP SPORE TRAITS TO CLEANED TREE:

viewL <- contMap(tree2, L)
#Maps average L for species to the tree.
#Most clades have L around the average for the tree (~7), with the Hemimycena cucullata + Mycena acicula clade and Typhyla (outgroup) having longer spores, and Baeospora and Cystostereum having shorter spores.

viewW <- contMap(tree2, W)
#Maps average W for species to the tree.
#Most clades have W around the average for the tree (~4), with the Henningsomyces + Rectipilus clades and Pseudohydropus having wider spores, and Cystostereum and Baeospora having narrower spores. Phloeomana is heterogeneous, with one clade having slightly narrower spores and its sister clade (including P. speirea) having wider spores.

viewQ <- contMap(tree2, Q)
#Maps average Q for species to the tree.
#Q varies quite a bit across the phylogeny, with Typhula, the Mycena acicula clade and Phloeomana sp. MO482167 having higher length/width ratios (more cylindrical) than average (ellipsoid), and the Henningsomyces + Rectipilus clade, Phloeomana clavata, Pseudohydropus and Clitocybula having lower length/width ratios (more globose) than average.

viewA <- contMap(tree2, Amyloid)
#Maps a binary trait, amyloidity to the tree (no rxn = 0, positive rxn = 1). 
#Demonstrates that clades differ in it. 
##Not surprising given that this is an important characteristic in separating genera and families. 
#This is included as more of a test of methods and will not be analyzed further below.

#COMPARISON OF SPORE DIMENSIONS MAPPED TO TREE:
##Based on visual comparison of the contMap outputs above, it appears that Q is more variable across the tree than L or W. 
##I will attempt to quantify these differences by comparing rates of evolution for each character below.
##There is some overlap in patterns across the tree between Q and both L and W. 
##This is to be expected given that Q is L/W. 
##L and W mostly exhibit different different patterns, with the exceptions being Cystostereum and Baeospora, which have shorter and narrower spores than average. 
##Their spores are smaller in both dimensions. 
##In the case of Cystostereum, the Q is roughly average, so these merely an average ellipsoid spore "scaled down".

#COMPARE BM MODELS FOR SPORE DIMENSIONS:

BML <- fitContinuous(tree2, L, model="BM")
print(BML)
# sigsq = 5.511491
##Rate of evolution of L under a Brownian Motion model is 5.511491 µm per time step.

##NOTE: I will admit that I am uncertain what the "time step" in this model is. 
##Base pair seems unlikely given that the input tree is based on an analysis of an alignment 2874 base pairs long while the spore dimensions as a whole have a minimum of 1 and a maximum of 11.


BMW <- fitContinuous(tree2, W, model="BM")
print(BMW)
# sigsq = 1.218592
##Rate of evolution of W under a Brownian Motion model is 1.218592 µm per time step.
#AIC = -8.883342

BMQ <- fitContinuous(tree2, Q, model="BM")
print(BMQ)
# sigsq = 0.411925
##Rate of evolution of Q under a Brownian Motion model is 0.411925 (unit-less ratio) per time step.
#AIC = -54.436977

##Rate of L evolution is greater than that of W and Q in that order. However, rates for L and W are in units of µm per [TIME STEP], while the rate for Q is unitless, which complicates comparison. 

#COMPARE OU MODELS FOR SPORE DIMENSIONS:

OUL <- fitContinuous(tree2, L, model="OU")
#Ignored "non-ultrametric" error message. Making trees ultrametric caused all pairs of OU and BM models to have deltaAIC values of exactly 2
print(OUL)
# sigsq = 6.235734
#Rate of evolution of L under an Ornstein-Uhlenbeck model is 6.235734 µm per time step.
#AIC = 	AIC = 55.867312

OUW <- fitContinuous(tree2, W, model="OU")
print(OUW)
#sigsq = 1.218573
#Rate of evolution of L under an Ornstein-Uhlenbeck model is 1.218573 µm per time step.
#	AIC = -6.883798

OUQ <- fitContinuous(tree2, Q, model="OU")
print(OUQ)
#sigsq = 0.434688
#Rate of evolution of L under an Ornstein-Uhlenbeck model is 0.434688 (unit-less ratio) per time step.
#AIC = -52.553407

##Rate of L evolution is greater than that of W and Q in that order. 
##However, rates for L and W are in units of µm per time step, while the rate for Q is unit-less, which complicates comparison. 

##For both OU and BM models, rate of L evolution is greater than that of W and Q in that order. 
##Again, with the caveat that Q is unitless while L and W are not. 
##It is perhaps not surprising that rates for Q would be lower for Q given that it is L/W.
##It could be argued that perhaps L is less constrained evolutionary than W and Q, but how these dimensions are defined would complicate this. 
##L is defined here as the longest spore diameter, with W being that roughly perpendicular to it, and Q being L divided by W. W can not be higher than L, and may not be able to vary as much because of this. 
##W can be lower than L and could theoretically vary more in that direction than L, but it seems likely that it would be easier for W to reach its absolute minimum viable diameter than for L to reach its absolute maximum. 
##Similarly, Q can not be below 1 (bounded by 1 and positive infinity).

##Q is perhaps more useful for differentiating genera than L and W. 
##Based on visual comparision of the 3 dimensions mapped to the tree, Q is most heterogeneous across the tree. 
##Given that its rate of evolution is also lower than that of L and W, it would also make sense to use it over those to compare genera, just as we use relatively strongly conserved genes to compare taxa at higher taxonomic levels than individual or species.


#COMPARE AIC VALUES:

AIC.BML <- 54.500498
AIC.OUL <- 55.867312
AIC.BMW <- -8.883342
AIC.OUW <- -6.883798
AIC.BMQ <- -54.436977
AIC.OUQ <- -52.553407

##The "best" model here is the one with the lowest AIC score. 
##In this case a Brownian motion model is best for L, W, and Q.
###NOTE: I am not entirely sure why this would be the case either biologically or statistically.

#Subtract the "best" AIC model from the model AIC value which we are interested in comparing:
delta.AIC.BML <- AIC.BML - AIC.BML
#delta.AIC.BML = 0 (subtracted from itself)
delta.AIC.OUL <- AIC.OUL - AIC.BML
#delta.AIC.BML = 1.366814
##The Brownian Motion and Ornstein-Uhlenbeck models for L are in the same model "family" (deltaAIC < 2) and not significantly different.

delta.AIC.BMW <- AIC.BMW - AIC.BMW
#delta.AIC.BMW = 0 (subtracted from itself)
delta.AIC.OUW <- AIC.OUW - AIC.BMW
#delta.AIC.OUW = 1.999544
##The Brownian Motion and Ornstein-Uhlenbeck models for W are in the same model "family" (deltaAIC < 2) and not significantly different.

delta.AIC.BMQ <- AIC.BMQ - AIC.BMQ
#delta.AIC.BMQ = 0 (subtracted from itself)
delta.AIC.OUQ <- AIC.OUQ - AIC.BMQ
#delta.AIC.OUQ = 1.88357
##The Brownian Motion and Ornstein-Uhlenbeck models for Q are in the same model "family" (deltaAIC < 2) and not significantly different.

##Overall, the Brownian Motion and Ornstein-Uhlenbeck models for all 3 spore dimensions of interest (L, W and Q) are in the same model family. 
##Either type of model could be used roughly interchangeably for future analyses of spore dimensions.