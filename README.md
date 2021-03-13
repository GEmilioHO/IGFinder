# IGFinder
**IGFinder**: Intronless genes finder for input species. 

IGFinder is a programme for identifying intronless genes for a particular organism. IGFinder works by using an ENSEMBL API. 

## Usage
python3 fetch_intronless_genes.py [optional species] (default *gallus gallus*) [/path/to/intron_UTR_db.txt]

## Files
<ul>
  <li> <b>fetch_intronless_genes.py</b>: &nbsp; Main script</li>
  <li> <b>fetch_genes.py</b>: &nbsp; Script containing functions required for <i>fetch_intronless_genes.py</i></li>
  <li> <b>intron_UTR_db.txt</b>: &nbsp; Required intron database of genes with introns in the UTR region for several species (NOT uploaded to GitHub). File should be downloaded from http://www.nextgenbioinformatics.org/IntronDB </li> 
</ul>

## Output
<ul>
  <li> <b>{species}-multi_exon_genes.tsv</b>: Contains data of multi-exon genes.</li>
  <li> <b>{species}-ui_single_exon_genes.tsv</b>: Contains data of UTR intron-containing single exon genes.</li>
  <li> <b>{species}-intronless_genes.tsv</b>: Contains data of intronless genes.</li>
</ul>

## Example
python3 fetch_intronless_genes.py "mus musculus" "/path/to/intron_UTR_db.txt"

## Citation
Aviña-Padilla, K., Rafael, J. A. R., Herrera-Oropeza, G. E., Muley, V., Valdivia, D., Valenzuela, E. D., ... & Hernández-Rosales, M. (2021). Evolutionary Perspective And Expression Analysis Of Intronless Genes Highlight The Conservation On Their Regulatory Role. bioRxiv.
