# intronless_genes
This script gets all the intronless genes of a particular organism.

### Files
fetch_intronless_genes.py &nbsp; -----> &nbsp; Main Python script <br>
fetch_genes.py &nbsp; -----> &nbsp; Script containing functions required for *fetch_intronless_genes.py*. <br>
intron_UTR_db.txt &nbsp; -----> &nbsp; Required intron database of genes with introns in the UTR region for the 7 species (NOT uploaded to GitHub). <br>
megs_segs_igs_bySpecie.zip &nbsp; -----> &nbsp; Output files. Three files (multi-exon genes, single-exon genes, and intronless genes) per specie (*danio rerio*, *mus musculus*, *gallus gallus*, *homo sapiens*, *monodelphis domestica*, *pan troglodytes*, and *rattus norvegicus*)

### Usage
python3 fetch_intronless_genes.py [optional species] (default *gallus gallus*)

### Example
python3 fetch_intronless_genes.py "mus musculus"
