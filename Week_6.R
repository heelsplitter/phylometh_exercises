library(geiger)
library(ape)
library(phangorn)

#You’ll need to get data into R in some way. Look at other phylometh assignments for how to get trees and data.

tree <- read.tree("C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Week_6_Discrete_Characters/RAxML_bipartitions.Porotheleaceae_multilocus_Trimal_2022_03_05_5_edit2.tre")
#plot.phylo(tree, type="fan", cex=0.5)

discrete.data <- read.csv(file="C:/Users/Owner/Desktop/EEB_587_Phylogenetic_Methods/phylometh_exercises/Week_6_Discrete_Characters/Tree_characters.csv", stringsAsFactors=FALSE) 
#death to factors.

rownames(discrete.data) <- discrete.data[,1]
discrete.data <- discrete.data[-1]
#print(discrete.data)


#Data are often not right in some way. They might not match the taxa in your tree, there may be missing data, etc. geiger::treedata is a great function for getting a tree and data that match, but your data may need other cleaning. Do it as a function so it’s repeatable.

CleanData <- treedata(phy=tree, data=discrete.data, sort=TRUE) 
#print(CleanData)

# Now write the code to use CleanData() to actually clean your data
#It’s critically important to LOOK at what you have. Are there weird values? Has the match between taxa and state gone correctly? Do you think you have binary data, but there’s actually only state 1? Especially as data sets grow (yay), and are assembled using scripts rather than by error-prone, non-reproducable hands (double yay), scientists are increasingly less likely to deeply look at our data. That’s bad – don’t be that person.
  
(VisualizeData <- CleanData) {

#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
  
# Now write the code to use VisualizeData() to actually look at your data
  
}
#First, let’s use parsimony to look at ancestral states:
  
cleaned.discrete.phyDat <- phangorn::phyDat(discrete.data, type="USER", levels=c(0,1)) #phyDat is a data format used by phangorn

anc.p <- phangorn::ancestral.pars(tree, cleaned.discrete.phyDat[,1])
plotAnc(tree, anc.p, 1)
#Do you see uncertainty? What does it mean?
  
#Now, plot the likelihood estimates.

anc.ml <- ancestral.pml(pml(tree, cleaned.discrete.phyDat), type="ml")
plotAnc(tree, anc.ml, 1)

#How does this differ from parsimony?

##Long branch attraction is less of a problem for maximum likelihood methods than for parsimony.
  
#Why does it differ from parsimony?

##Maximum likelihood differs from parsimony in being a model-based optimality criterion and generally having only structural parameters rather than incidental parameters.
  
#What does uncertainty mean?

##Uncertainty in the context of phylogenies is often equivalent to lack of resolution. Polytomies or nodes with low bootstrap support would be examples of this.
  
#Now, to the biological questions. For many of these, corHMM will be a useful package. Do the following analyses:
 

#How can you estimate transition rates between states? Do it.

##Substitution rate could be correlated to the apparent state-transition occuring at nodes. I am unable to test this with my data currently.

#How could you examine if transition rates are equal?

##You could construct a model where all transition rates are equal, and one or several where rates are allowed to vary and then run analyses, potentially maximum likelihood to compare the models. However, I am unsure what statistical parameters would be most useful in comparing the resulting outputs. I would expect that transition rates would not be equal between lineages with disparity in branch lengths, and these could possibly be compared as proxies.

#Think about the Lewis (2001) MKV model. Are your traits all variable? Will using this make sense for your data? Try using it. Do results change?

##For this analysis. I attempted to analyze a single variable trait, so I would assume that this would make sense for my data. However, given that I was unable to perform the prior analyses required for further testing, I am unable to answer the rest of this question.

#How could you test order of state evolution?

##I can imagine 2 different options, the viability of which I am unsure of. One would involve constructing models for all possible transition series for the traits being analyes. For example, if you have characters numbered 1 through 5, you could have models for 1->2->3->4->5, 2->4->1->3->5, etc. and possibly use ML to determine which of these fits the phylogeny best, but this seems computationally expensive and impractical. You could also infer a series of ancestral states following the nodes sequentially back to the root, but these would merely be inferred and it would be difficult to test this. I would assume there are statistical tests that are relevant here, but I suspect that they would be more useful in determining which models of state changes are more likely rather than actually testing the history of changes in the taxon being analyzed. We are still dealing with inferred historical events. If we examined phenotypes in bacterial or viral populations and froze various ancestor populations, we could then actually test transitions by referring back to frozen ancestors. Any methods used here would probably be good for eukaryotes.