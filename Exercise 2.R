```{r, eval=TRUE}

GetTreeFromOpenTree <- function() {
  library(rotl)
  library(ape)
  
  formica.id <- rotl::tnrs_match_names("Agaricomycetes")$ott_id
  
  # Now get Open Tree's current best estimate of the phylogeny for the group
  # They call this the tree of life; we can get the subtree for just this group.
  formica.tree <- rotl::tol_subtree(ott_id=formica.id)
  
  # Let's plot the tree:
  ape::plot.phylo("Agaricomycetes", type="fan", cex=0.2)
}

GetTreeFromOpenTree()