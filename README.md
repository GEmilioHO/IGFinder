[![DOI](https://zenodo.org/badge/320059139.svg)](https://zenodo.org/badge/latestdoi/320059139)
# IGFinder
**IGFinder**: Intronless genes finder for an input species. 

IGFinder is a script for identifying intronless genes for a particular organism. IGFinder works by using an ENSEMBL API. 

## Usage
python3 fetch_intronless_genes.py [optional species] (default *gallus gallus*) [/path/to/intron_UTR_db.txt]

## Files
<ul>
  <li> <b>fetch_intronless_genes.py</b>: &nbsp; Main script.</li>
  <li> <b>fetch_genes.py</b>: &nbsp; Script containing helper functions for <i>fetch_intronless_genes.py</i></li>
  <li> <b>intron_UTR_db.txt</b>: &nbsp; Required intron database of genes with introns in the UTR region for several species (NOT uploaded to GitHub). The corresponding Intron file can be downloaded from https://drive.google.com/file/d/1c1t4eMQZ4SSHD0-iUiHph8eivTDYUq0l/view?usp=share_link. </li>
</ul>

## Output
<ul>
  <li> <b>{species}-multi_exon_genes.tsv</b>: &nbsp; Contains data of multi-exon genes.</li>
  <li> <b>{species}-ui_single_exon_genes.tsv</b>: &nbsp; Contains data of UTR intron-containing single exon genes.</li>
  <li> <b>{species}-intronless_genes.tsv</b>: &nbsp; Contains data of intronless genes.</li>
</ul>

## Example
python3 fetch_intronless_genes.py "mus musculus" "/path/to/intron_UTR_db.txt"

## Citation
Aviña-Padilla, K., Ramírez-Rafael, J. A., Herrera-Oropeza, G. E., Muley, V. Y., Valdivia, D. I., Díaz-Valenzuela, E., García-García, A., Varela-Echavarría, A., & Hernández-Rosales, M. (2021). Evolutionary Perspective and Expression Analysis of Intronless Genes Highlight the Conservation of Their Regulatory Role. Frontiers in genetics, 12, 654256. https://doi.org/10.3389/fgene.2021.654256
