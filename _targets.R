# Format iTOL annotation files for fern family and order

library(targets)
library(tarchetypes)

source("R/functions.R")
source("R/packages.R")

tar_option_set(
  imports = "ftolr"
)

tar_plan(
  # Load tree
  fern_tree = ft_tree(drop_og = TRUE),
  # Load taxa
  fern_taxa = load_fern_taxonomy(fern_tree),
  # Assess monophyly
  # - overall
  fern_monophy = assess_monophy(
    taxon_sampling = fern_taxa,
    tree = fern_tree,
    tax_levels = c("family", "suborder", "order")
  ),
  # - family
  family_monophy = get_monophy(
    fern_monophy, "family"
  ),
  # - suborder
  suborder_monophy = get_monophy(
    fern_monophy, "suborder"
  ),
  # - order
  order_monophy = get_monophy(
    fern_monophy, "order"
  ),
  # Get spanning tips for each monophyletic group
  family_tips = get_monophy_spanning_tips(
    family_monophy, "family", fern_tree
  ),
  suborder_tips = get_monophy_spanning_tips(
    suborder_monophy, "suborder", fern_tree
  ),
  order_tips = get_monophy_spanning_tips(
    order_monophy, "order", fern_tree
  ),
  # Set colors for monophyletic groups
  family_colors = set_clade_colors(
    family_tips, family
  ),
  suborder_colors = set_clade_colors(
    suborder_tips, suborder
  ),
  order_colors = set_clade_colors(
    order_tips, order
  ),
  # Write out iTOL annotation files
  # - family colored strips
  tar_file(
    family_colored_strips,
    write_itol_colored_strips(
      family_colors,
      dataset_lab = "Family strips",
      dataset_col = "#2986cc", # dark blue
      "_targets/user/itol/01_family_strip_cols.txt"
    )
  ),
  # - suborder colored strips
  tar_file(
    suborder_colored_strips,
    write_itol_colored_strips(
      suborder_colors,
      dataset_lab = "Suborder strips",
      dataset_col = "#8fce00", # dark green
      "_targets/user/itol/02_suborder_strip_cols.txt"
    )
  ),
  # - order colored strips
  tar_file(
    order_colored_strips,
    write_itol_colored_strips(
      order_colors,
      dataset_lab = "Order strips",
      dataset_col = "#e06666", # dark red
      "_targets/user/itol/03_order_strip_cols.txt"
    )
  ),
  # - family colored branches
  tar_file(
    family_colored_branches,
    write_itol_colored_branches(
      family_colors,
      dataset_lab = "Family branches",
      dataset_col = "#d0e0e3", # light blue
      "_targets/user/itol/04_family_branch_cols.txt"
    )
  ),
  # - suborder colored branches
  tar_file(
    suborder_colored_branches,
    write_itol_colored_branches(
      suborder_colors,
      dataset_lab = "Suborder branches",
      dataset_col = "#d9ead3", # light green
      "_targets/user/itol/05_suborder_branch_cols.txt"
    )
  ),
  # - order colored branches
  tar_file(
    order_colored_branches,
    write_itol_colored_branches(
      order_colors,
      dataset_lab = "Order branches",
      dataset_col = "#f4cccc", # light red
      "_targets/user/itol/06_order_branch_cols.txt"
    )
  ),
  # - tree
  tar_file(
    fern_tree_file,
    write_tree(
      fern_tree,
      file = glue::glue(
        "_targets/user/itol/ftol_con_dated_ogdrop_v{ft_data_ver()}.tree")
    )
  ),
  # Zip files
  # (only needed for batch uploader, which isn't working yet)
  tar_file(
    itol_zip_file,
    zip_files(
      "_targets/user/itol.zip",
      c(family_colored_strips, order_colored_strips,
      family_colored_branches, order_colored_branches, fern_tree_file)
    )
  )
)