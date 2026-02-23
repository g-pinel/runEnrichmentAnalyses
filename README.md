# runEnrichmentAnalyses

### Easily run enrichment analyses against the Gene Ontology, KEGG, Wiki Pathways and Reactome databases in R

After performing differential expression analyses to obtain a list of differentially-expressed genes (DEGs), it is common to run enrichment analyses to identify groups of genes that participate in the same process or pathway and are enriched in one of the groups being compared.

The functions in this repository allow you to easily run enrichment analyses against Gene Ontology, KEGG, Wiki Pathways and Reactome databases using DESeq2's DeseqDataSet (DDS) objects as inputs.

Multiple DDS objects can be provided in batch to perform multiple sets of enrichments automatically.

### Dependencies

These function require the tidyverse, clusterProfiler, ReactomePA and organism database (e.g. org.Hs.eg.db for human data, or org.Mm.eg.db for mouse data) R packages, which can be installed in R with:

```
install.packages("tidyverse")
```
and
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "ReactomePA", "org.Hs.eg.db", "org.Mm.eg.db"))
```
<br>

## batchRunEnrichmentAnalyses function (run multiple sets of enrichment analyses in batch)

This is a wrapper function that passes all the arguments to the *runEnrichmentAnalyses* function looping over all the desired comparisons. It uses the following parameters:
<br><br>

- **dds_list**: a list of length 1 or more containing the input DDS objects. <ins>Must be a named list</ins>, as the object names will be used to name the result directories.
- **gene_id_format**: the type of gene IDs found in the DDS objects (string). If running Gene Ontology, must match a keytype listed when running (example for human data):
```
keytypes(org.Hs.eg.db)
```
- **counts**: a matrix with gene IDs as row names and samples as column names containing raw or gene-length-adjusted counts. Will be used to obtain the background universe.
- **entrez_df**: a data frame containing the gene IDs found in the DDS object (defined by  *gene_id_format*), and Entrez (NCBI) identifiers. <ins>Must contain Entrez IDs in a column called 'ENTREZID'</ins>.
- **org_db**: Annotation package for the target organism (e.g. org.Hs.eg.db for human data, or org.Mm.eg.db for mouse data).
- **outdir**: the path to the directory where all result files will be output.

### batchRunEnrichmentAnalyses call example
```
enrichments <- batchRunEnrichmentAnalyses(dds_list,
                                          gene_id_format = "ENSEMBL", 
                                          counts = counts_processed, 
                                          entrez_df = symbol2ensembl2entrez, 
                                          org_db = org.Mm.eg.db,
                                          outdir = "results/ORA/")
```

> [!NOTE]
> This function by default runs all types of supported enrichment analyses (Gene Ontology, KEGG, Wiki Pathways and Reactome). To disable some of them, set one/multiple of the following optional parameters to FALSE:<br><br>
> GOBP, KEGG, WIKIPATHWAYS, REACTOME
<br>

## runEnrichmentAnalyses function
This is the function that performs the enrichment analyses. Can be run separately from *batchRunEnrichmentAnalyses* if you only wish to run one set of enrichments. It performs the following:
- Creates a results directory/detects an existing one matching the input DDS object name inside the output directory.
- Creates a log file with useful information about the run.
- Retrieves significant DEGs for each group.
- Filters the counts matrix for the samples involved in the comparison and obtains the background universe of genes.
- Runs hypergeometric tests (classic enrichment analyses).
- Writes CSV files with the results in the results directory.
- Returns a list of data frames with results.

> [!IMPORTANT]
> Files corresponding to DEGs with a log2FoldChange > 0 will contain 'pos' in their filename, while those with log2FC < 0 will contain 'neg'.
