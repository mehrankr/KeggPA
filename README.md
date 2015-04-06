![KeggPA](http://www.genome.jp/Fig/kegg128.gif)
# KeggPA
[Kegg pathway annotation database](http://www.genome.jp/kegg/pathway.html) overrepresentation analysis by Fisher's exact test.

This R package requires you to provide a pathway annotation dataset with the following columns:

Pathway.ID | Entrez | HGNC | Name.Pathway

path:hsa00010 | hsa:10327 | AKR1A1 | Glycolysis / Gluconeogenesis

...


Without providing such a dataframe for kegg.list argument, an old version will be used automatically.




## Examples

```R
require(KeggPA)
kpa = keggPA(hgnc=c("TP53", "BRCA1", "BRCA2", "VHL"), proj_name="example_genes")
```


### Note

This R packages is neither supported nor endorsed by Kegg pathway database.



