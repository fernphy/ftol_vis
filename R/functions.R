#' Get tips of a phylogenetic tree in their plotted order
#'
#' After re-rooting a tree, the order of tips when the tree
#' is plotted no longer match the order of $tip.label. Use
#' this function to get tips in the order they are plotted.
#' @param tree List of class "phylo"
#' @return Character vector
get_tips_in_ape_plot_order <- function(tree) {
  assertthat::assert_that(inherits(tree, "phylo"))
  # First filter out internal nodes
  # from the the second column of the edge matrix
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  # Use this vector to extract the tips in the right order
  tree$tip.label[ordered_tips]
}

#' Get a pair of tips that define a clade in a phylogenetic tree
#'
#' @param tree Phylogenetic tree, must be rooted.
#' @param node Number of a node in the tree.
#'
#' @return Character vector; a pair of tips whose MRCA is `node`
#'
get_spanning_tips <- function(tree, node) {
  # Tree must be rooted
  assertthat::assert_that(ape::is.rooted(tree))
  # Ladderize tree
  tree <- ape::ladderize(tree)
  # Get vector of tips in ladderized order
  all_tips_ladder_ord <- get_tips_in_ape_plot_order(tree)
  # Get spanning tips, not in ladderized order
  spanning_tips_unord <- phangorn::Descendants(tree, node, "tips") %>%
    magrittr::extract2(1) %>%
    magrittr::extract(tree$tip.label, .)
  # Put 'final' spanning tips in ladderized order
  final_tips <-
    all_tips_ladder_ord[all_tips_ladder_ord %in% spanning_tips_unord]
  # Return first and last spanning tips in ladderized order
  if (length(final_tips) > 2) return(
    c(dplyr::first(final_tips), dplyr::last(final_tips)))
  final_tips
}

#' Assess monophyly
#'
#' Wrapper around MonoPhy::AssessMonophyly()
#'
#' @param taxon_sampling Dataframe of taxa to assess for monophyly. Must
#' include column "species"
#' @param tree Phylogenetic tree.
#' @param og_taxa Character vector; pair of taxa to define the outgroup to
#' root the tree.
#' @param tax_levels Character vector; names of columns in `taxon_sampling`
#' to check for monophyly.
#'
#' @return List; results of MonoPhy::AssessMonophyly()
#'
assess_monophy <- function(
  taxon_sampling, tree,
  og_taxa = NULL,
  tax_levels) {
  tax_levels <- c("species", tax_levels) %>% unique()
  # Root tree
  if (!is.null(og_taxa)) {
    tree <- phytools::reroot(
      tree,
      getMRCA(tree, og_taxa)
    )
  }
  # Check monophyly
  taxon_sampling %>%
    assertr::verify("species" %in% colnames(.)) %>%
    select(species, all_of(tax_levels)) %>%
    as.data.frame() %>%
    MonoPhy::AssessMonophyly(tree, .)
}