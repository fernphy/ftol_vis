# Format iTOL annotation files for fern family and order
# use "colored strips" dataset type https://itol.embl.de/help.cgi#strip

library(ftolr)
library(tidyverse)
library(janitor)
library(assertr)

source("R/functions.R")

# Tree URL on iTOL:
# https://itol.embl.de/tree/1191064828238691658722151

# Load tree
fern_tree <- ft_tree(drop_og = TRUE)

# Load fern taxonomy
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

fern_taxa <-
  fern_taxa %>%
  filter(species != "Dryopolystichum_phaeostigma") %>%
  bind_rows(new_dryopoly) %>%
  arrange(outgroup, order, suborder, family, subfamily, genus)

# Assess monophyly at family and order level
fern_monophy <- assess_monophy(
  taxon_sampling = fern_taxa,
  tree = fern_tree,
  tax_levels = c("family", "order")
)

family_monophy <-
fern_monophy$family$result %>%
  rownames_to_column("family") %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(mrca = parse_number(mrca)) %>%
  # All families should be monophyletic or monotypic
  verify(all(monophyly %in% c("Monotypic", "Yes")))

order_monophy <-
fern_monophy$order$result %>%
  rownames_to_column("order") %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(mrca = parse_number(mrca)) %>%
  verify(all(monophyly %in% c("Monotypic", "Yes")))

# Obtain representative tips for plotting: pairs for monophyletic, single
# species for monotypic
family_tips <-
family_monophy %>%
  select(family, monophyly, mrca) %>%
  mutate(
    rep_tips = map(mrca, ~get_spanning_tips(fern_tree, .)),
    rep_tips = map_chr(rep_tips, ~paste(., collapse = "|"))
  )

order_tips <-
order_monophy %>%
  select(order, monophyly, mrca) %>%
  mutate(
    rep_tips = map(mrca, ~get_spanning_tips(fern_tree, .)),
    rep_tips = map_chr(rep_tips, ~paste(., collapse = "|"))
  )

# Specify colors for labels, write out datafile formatted for iTOL
family_tips_data <-
  family_tips %>%
  filter(monophyly == "Yes") %>% # monotypic labels not showing up anyways
  mutate(
    color = Polychrome::createPalette(
      nrow(.), c("#00ffff", "#ff00ff", "#ffff00"),
      M = 5000)
  ) %>%
  transmute(rep_tips = rep_tips, color, label = family)

order_tips_data <-
  order_tips %>%
  filter(monophyly == "Yes") %>%
  mutate(
    color = Polychrome::createPalette(
      nrow(.), c("#00ffff", "#ff00ff", "#ffff00"),
      M = 5000)
  ) %>%
  transmute(rep_tips = rep_tips, color, label = order)

write_tsv(family_tips_data, "family_cols.txt", col_names = FALSE)

tsv_text <- read_lines("family_cols.txt")

write_lines(
  c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tFamily",
    "COLOR\t#2986cc",
    "COLOR_BRANCHES\t1",
    "SHOW_INTERNAL\t1", # show labels for internal nodes
    "STRIP_WIDTH\t40",
    "SHOW_STRIP_LABELS\t0", # hide text on labels
    "DATA",
    tsv_text
  ),
  "family_cols.txt"
)

write_tsv(order_tips_data, "order_cols.txt", col_names = FALSE)

tsv_text <- read_lines("order_cols.txt")

write_lines(
  c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    "DATASET_LABEL\tOrder",
    "COLOR\t#e06666",
    "COLOR_BRANCHES\t0", # Don't color branches
    "SHOW_INTERNAL\t1", # show labels for internal nodes
    "STRIP_WIDTH\t40",
    "SHOW_STRIP_LABELS\t0", # show text on labels
    "STRIP_LABEL_SIZE_FACTOR\t2",
    "DATA",
    tsv_text
  ),
  "order_cols.txt"
)

# Write out internal node labels
# (though I can't get these to actually show up on the tree)
order_tips %>%
  rename(taxon = order) %>%
  bind_rows(
    rename(family_tips, taxon = family)
  ) %>%
  filter(monophyly == "Yes") %>%
  transmute(text = glue::glue("{rep_tips} = node {taxon}")) %>%
  write_tsv("node_labels.txt", col_names = FALSE)

tsv_text <- read_lines("node_labels.txt")

write_lines(
  c(
    "LABELS",
    "SEPARATOR SPACE",
    "DATA",
    tsv_text
  ),
  "node_labels.txt"
)

# Write out tree
ape::write.tree(fern_tree, "ftol_con_dated_ogdrop_v1.1.0.tre")
