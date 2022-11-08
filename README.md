# ftol_vis

Code to generate annotation files for displaying metadata on [fern tree of life (FTOL)](https://fernphy.github.io/) using the [Interactive Tree of Life viewer (iTOL)](https://itol.embl.de/shared/iwasaki_utokyo).

To generate the [annotation files](_targets/user/itol) (.txt), run `_targets.R` in R.

## Uploading data

Manually upload the tree in the iTOL user interface.

Then, open the tree and click on the "Datasets" tab. Select all dataset files (`*.txt` files in `_targets/user/itol`) and upload them.

A URL ending a long number will automatically be generated (e.g., [https://itol.embl.de/tree/1191064828398171667872880](https://itol.embl.de/tree/1191064828398171667872880)).

## License

Code: [MIT](LICENSE)