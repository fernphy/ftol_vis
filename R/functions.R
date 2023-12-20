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

# Format JSON for taxonium tree viewer
format_ftol_json <- function(current_release, current_version) {
  tibble::tibble(
    source = glue::glue("GenBank release {current_release}"),
    title = glue::glue("FTOL v{current_version}")
  )
}

write_csv_file <- function(data, file, ...) {
  readr::write_csv(x = data, file = file, ...)
  file
}

write_json_file <- function(json, file, ...) {
  jsonlite::write_json(x = json, path = file, ...)
  file
}