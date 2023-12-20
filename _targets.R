# Format annotation files for plotting fern family and order on tree
# in taxonium

library(targets)
library(tarchetypes)
library(ftolr)

source("R/functions.R")

tar_plan(
  # Load tree
  fern_tree = ftolr::ft_tree(drop_og = TRUE),
  # Load taxa
  fern_taxa = load_fern_taxonomy(fern_tree),
  # Format tip metadata
  fern_tip_data = dplyr::select(
    fern_taxa, species, genus, subfamily, family, suborder, order),
  # Format json config
  tar_target(
    current_release,
    ft_data_ver("gb"),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    current_version,
    ft_data_ver("ftol"),
    cue = tar_cue(mode = "always")
  ),
  fern_config = format_ftol_json(
    current_release = current_release,
    current_version = current_version),
  # Write files
  # - tree
  tar_file(
    fern_tree_file,
    write_tree_file(fern_tree, file = "_targets/user/taxonium/ftol_tree.tree"),
  ),
  # - metadata
  tar_file(
    fern_tip_data_file,
    write_csv_file(
      fern_tip_data, file = "_targets/user/taxonium/ftol_data.csv"),
  ),
  # - json config
  tar_file(
    fern_config_file,
    write_json_file(
      fern_config, file = "_targets/user/taxonium/ftol_config.json")
  )
)
