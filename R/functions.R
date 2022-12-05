# Load fern taxonomy trimmed to tree
load_fern_taxonomy <- function(fern_tree) {
  ftolr::ftol_taxonomy |>
    dplyr::filter(species %in% fern_tree$tip.label)
}

# Write out a phylogenetic tree and return the path to the tree
write_tree_file <- function(tree, file, ...) {
  ape::write.tree(tree, file, ...)
  file
}

write_csv_file <- function(data, file, ...) {
  readr::write_csv(x = data, file = file, ...)
  file
}

write_json_file <- function(json, file, ...) {
  jsontools::write_json(x = json, path = file, ...)
  file
}