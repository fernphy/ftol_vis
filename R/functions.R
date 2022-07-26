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

# Load fern taxonomy trimmed to tree, fix Dryopolystichum_phaeostigma
load_fern_taxonomy <- function(fern_tree) {
  fern_taxa <- ftolr::ftol_taxonomy %>%
    filter(species %in% fern_tree$tip.label)

  # Work-around until https://github.com/fernphy/ftol/issues/16 is fixed:
  # change Dryopolystichum_phaeostigma to Lomariopsidaceae
  # or Dryopteridaceae will show up as non-monophyletic
  new_dryopoly <-
    fern_taxa %>%
    filter(str_detect(species, "Lomario")) %>%
    slice(1) %>%
    mutate(species = "Dryopolystichum_phaeostigma")

  fern_taxa %>%
    filter(species != "Dryopolystichum_phaeostigma") %>%
    bind_rows(new_dryopoly) %>%
    arrange(outgroup, order, suborder, family, subfamily, genus)
}

#' Get monophy results for a specific taxonomic level
#'
#' Taxonomic level must match one of those used in input
#' taxonomic table for assess_monophy()
#'
#' @param monophy_list Output of assess_monophy().
#' @param tax_level Taxonomic level (e.g., "family").
#'
#' @return Tibble; monophyly results for that taxonomic level
get_monophy <- function(monophy_list, tax_level) {
  monophy_list[[tax_level]][["result"]] %>%
    rownames_to_column(tax_level) %>%
    as_tibble() %>%
    clean_names() %>%
    mutate(mrca = parse_number(mrca)) %>%
    # All families should be monophyletic or monotypic
    verify(all(monophyly %in% c("Monotypic", "Yes")))
}


#' Get a pair of tips spanning each monophyletic clade
#'
#' @param taxon_monophy Output of get_monophy().
#' @param tax_level Taxonomic level (e.g., "family").
#' @param tree Phylogenetic tree.
#'
#' @return Tibble with a column "rep_tips" including a pair of tips
#' that span each monophyletic clade. Species names separated by pipe symbol (|)
get_monophy_spanning_tips <- function(taxon_monophy, tax_level, tree) {
  taxon_monophy %>%
    select({{tax_level}}, monophyly, mrca) %>%
    mutate(
      rep_tips = map(mrca, ~get_spanning_tips(tree, .)),
      rep_tips = map_chr(rep_tips, ~paste(., collapse = "|"))
    )
}

#' Set clade colors
#'
#' @param spanning_tips Output of get_monophy_spanning_tips().
#' @param tax_level Taxonomic level (e.g., "family").
#'
#' @return Tibble
set_clade_colors <- function(spanning_tips, tax_level) {
  spanning_tips %>%
    filter(monophyly == "Yes") %>% # monotypic labels not showing up anyways
    mutate(
      color = Polychrome::createPalette(
        nrow(.), c("#00ffff", "#ff00ff", "#ffff00"),
        M = 5000)
    ) %>%
    transmute(rep_tips = rep_tips, color, label = {{tax_level}})
}

#' Write an iTOL annotation file to add colored strips to clades
#' 
#' https://itol.embl.de/help.cgi#strip
#'
#' @param clade_colors Output of set_clade_colors().
#' @param file Path to write annotation file.
#' @param dataset_lab Label to use for dataset.
#' @param dataset_col Color to use for dataset.
#' @return Path to output file
#'
write_itol_colored_strips <- function(clade_colors, file,
  dataset_lab, dataset_col) {
  # Write initial data as tsv file
  readr::write_tsv(clade_colors, file, col_names = FALSE)

  # Read back in, append header
  tsv_text <- readr::read_lines(file)

  readr::write_lines(
    c(
      "DATASET_COLORSTRIP",
      "SEPARATOR TAB",
      glue::glue("DATASET_LABEL\t{dataset_lab}"),
      glue::glue("COLOR\t{dataset_col}"),
      "COLOR_BRANCHES\t0", # don't color branches
      "SHOW_INTERNAL\t0", # show labels for internal nodes
      "STRIP_WIDTH\t40",
      "SHOW_STRIP_LABELS\t0", # hide text on labels
      "DATA",
      tsv_text
    ),
    file
  )

  # Return path to output file
  file
}

#' Write an iTOL annotation file to add colored branches to clades
#'
#' https://itol.embl.de/help.cgi#styleds
#'
#' @param clade_colors Output of set_clade_colors().
#' @param file Path to write annotation file.
#' @param dataset_lab Label to use for dataset.
#' @param dataset_col Color to use for dataset.
#' @return Path to output file
#'
write_itol_colored_branches <- function(clade_colors, file,
  dataset_lab, dataset_col) {
  # Write initial data as tsv file
  clade_colors %>%
    transmute(text = glue::glue("{rep_tips},branch,clade,{color},1,normal")) %>%
    readr::write_tsv(file, col_names = FALSE)

  # Read back in, append header
  tsv_text <- readr::read_lines(file)

  readr::write_lines(
    c(
      "DATASET_STYLE",
      "SEPARATOR COMMA",
      glue::glue("DATASET_LABEL,{dataset_lab}"),
      glue::glue("COLOR,{dataset_col}"), # light blue
      "DATA",
      tsv_text
    ),
    file
  )

  # Return path to output file
  file
}

#' Write an iTOL annotation file to add node labels
#'
#' See "Changing the leaf labels and internal node names" under iTOL help
#' https://itol.embl.de/help.cgi
#'
#' @param spanning_tips_list Output of get_spanning_tips().
#' @param file Path to write annotation file.
#'
write_itol_node_labels <- function(spanning_tips_list, file) {
  do.call(dplyr::bind_rows, spanning_tips_list) %>%
    # Restrict to only monophyletic nodes
    filter(monophyly == "Yes") %>%
    transmute(text = glue::glue("{rep_tips},{taxon}")) %>%
    write_tsv(file, col_names = FALSE)

  tsv_text <- readr::read_lines(file)

  readr::write_lines(
    c(
      "LABELS",
      "SEPARATOR COMMA",
      "DATA",
      tsv_text
    ),
    file
  )
  # Return path to output file
  file
}

# Write out a phylogenetic tree and return the path to the tree
write_tree <- function(tree, file) {
  ape::write.tree(tree, file)
  file
}

# Zip files, junking directories and deleting any existing output
zip_files <- function(zipfile, files) {
  if (fs::file_exists(zipfile)) fs::file_delete(zipfile)
  zip(zipfile, files, flags = "-j")
  zipfile
}