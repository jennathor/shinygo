# ShinyGO: Gene Set Enrichment Analysis and Visualization

**A graphical tool for functional enrichment analysis supporting over 600 plant and animal species plus 14,000 genomes from STRING-DB**

ShinyGO is a web-based application that provides comprehensive gene set enrichment analysis with interactive visualizations. It's designed to be accessible to all biologists, regardless of bioinformatics expertise.

## Quick Access

**Main Server:** [http://bioinformatics.sdstate.edu/go](http://bioinformatics.sdstate.edu/go)
**Mirror Server (if busy):** [http://ge-lab.org/go/](http://149.165.154.220/go/)
**Customized Version:** [http://bioinformatics.sdstate.edu/goc/](http://bioinformatics.sdstate.edu/goc/)
**Beta Version (v0.80):** [http://bioinformatics.sdstate.edu/go80/](http://bioinformatics.sdstate.edu/go80/)

## Overview

ShinyGO simplifies functional enrichment analysis by:
- Accepting most gene ID formats and automatically converting them to Ensembl IDs
- Auto-detecting species from your gene list
- Performing enrichment analysis across multiple pathway databases
- Generating publication-ready visualizations
- Providing interactive networks and hierarchical trees
- Integrating with KEGG pathway diagrams and STRING-DB

## Key Features

### Gene Analysis
- **Universal Gene ID Support**: Accepts Ensembl IDs, gene symbols, Entrez IDs, RefSeq IDs, and more
- **Automatic Species Detection**: Best-match algorithm identifies your species automatically
- **Custom Background Genes**: Upload your own background set for more accurate enrichment (highly recommended for RNA-seq data)
- **600+ Species Support**:
  - 215 main species from Ensembl
  - 177 metazoa species (Ensembl Metazoa)
  - 124 plant species (Ensembl Plants)
  - 33 protists, 1 bacteria
  - 14,094 additional species from STRING-DB v11.5 (bacteria, fungi, etc.)

### Pathway Databases

**Gene Ontology (GO)**:
- Biological Process (~15,796 gene sets)
- Cellular Component (~1,916 gene sets)
- Molecular Function (~4,605 gene sets)

**KEGG**: ~327 pathways (varies by species)

**Curated Databases** (for major species like human/mouse):
- BioCarta
- Reactome
- WikiPathways
- GeneSetDB
- And many more

### Enrichment Analysis Features
- **Statistical Testing**: Hypergeometric test with FDR (False Discovery Rate) correction
- **Flexible Sorting**: Sort by FDR, Fold Enrichment, gene count, pathway name, or combined metrics
- **Redundancy Removal**: Automatically filters similar pathways sharing 95% of genes
- **Pathway Size Filtering**: Customize minimum and maximum pathway sizes
- **Customizable Thresholds**: Adjust FDR cutoff, number of pathways displayed, and more

## Visualizations and Outputs

### 1. Enrichment Table
Interactive table displaying:
- Pathway names with links to source databases
- FDR values and fold enrichment
- Number of genes per pathway
- Complete gene lists for each enriched pathway

### 2. Enrichment Chart
Customizable publication-ready plots:
- Three chart types: lollipop, dotplot, barplot
- Adjustable axes (FDR, Fold Enrichment, Gene Count)
- Customizable colors, fonts, sizes, and themes
- Download in PDF, PNG, or SVG format

### 3. Hierarchical Clustering Tree
- Visualizes relationships between enriched pathways
- Groups pathways sharing many genes
- Circle size indicates statistical significance
- Helps identify major functional themes

### 4. Interactive Network Diagram
- Nodes represent enriched pathways
- Edges connect pathways sharing genes (adjustable threshold)
- Draggable, zoomable interface
- Download as HTML or static image
- Export nodes and edges data

### 5. KEGG Pathway Diagrams
- Your genes highlighted in red on official KEGG pathway maps
- Visualize genes in biological context
- Available when KEGG database is selected

### 6. Gene Characteristics Plots
Compare your genes with the genome background:

**Density Plots**:
- Coding sequence length
- Transcript length
- Genome span
- 5' UTR and 3' UTR lengths
- GC content

**Bar Plots**:
- Chromosome distribution
- Gene type distribution
- Exon count
- Transcript count

Statistical tests (t-test, chi-squared) show if your genes have special characteristics.

### 7. Genome Browser Plot
- Interactive chromosome-level visualization
- Identifies genomic regions enriched with your genes
- Sliding window analysis with statistical testing
- Zoom and explore regions of interest
- Available for Ensembl-annotated species

### 8. Promoter Motif Analysis
- Identifies enriched transcription factor binding motifs
- Analyzes 300bp or 600bp upstream sequences
- Indicates if TF gene is in your query list
- Statistical testing against genome background

### 9. STRING-DB Integration
- Protein-protein interaction (PPI) network retrieval
- Independent enrichment analysis via STRING API
- Interactive network visualization
- Links to STRING website for detailed exploration
- Multiple functional categories (GO, KEGG, Pfam, InterPro)

### 10. Gene Grouping
- Groups genes by high-level GO categories
- Available for GO Biological Process analysis
- Helps organize genes into functional themes

## How to Use ShinyGO

### Basic Workflow

1. **Paste Your Gene List**
   - Enter genes in the text area (one per line, or separated by spaces, tabs, or commas)
   - Most gene ID formats are accepted
   - Click "Demo genes" to try an example

2. **Verify Species**
   - ShinyGO auto-detects your species
   - Verify the suggested species is correct
   - Manually select if needed using the search box
      - Search by common name, scientific name, or NCBI taxonomy ID to find your species

3. **Add Background Genes (Recommended)**
   - Click "Background (recommended)"
   - Paste all genes from your experiment (e.g., all genes with detectable expression in RNA-seq)
   - This significantly improves accuracy

4. **Submit Analysis**
   - Click the "Submit" button
   - Wait for analysis to complete

5. **Explore Results**
   - Browse enrichment table
   - Switch between visualization tabs
   - Adjust parameters and resort as needed
   - Download results and figures

### Advanced Options

**Pathway Database Selection**:
- Choose specific databases (GOBP, GOCC, GOMF, KEGG, etc.)
- Or analyze all available databases

**FDR Cutoff**:
- Default: 0.05
- Adjust based on stringency needs
- Note: Really significant FDRs are typically 1E-5 to 1E-20

**Number of Pathways**:
- Display top 10-200 pathways
- Download all significant pathways regardless of display limit

**Pathway Size Limits**:
- Minimum: 2-30 genes (default: 2)
- Maximum: 1000-5000 genes (default: 2000)
- Smaller pathways can introduce noise
- Larger pathways may be less informative

**Sorting Options**:
- Sort by FDR (most significant first)
- Sort by Fold Enrichment (most enriched first)
- Select by FDR, sort by Fold Enrichment (recommended)
- Sort by average ranks
- Sort by gene count or pathway name

**Redundancy Removal**:
- Enabled by default
- Removes pathways sharing 95% of genes
- Keeps most significant pathway from redundant group

**Pathway Name Abbreviation**:
- Enabled by default
- Shortens long pathway names (e.g., "Positive regulation" → "Pos. reg.")

## Input Format

### Accepted Gene IDs
- Ensembl gene IDs (ENSG00000..., ENSMUSG00000..., etc.)
- Gene symbols (TP53, BRCA1, etc.)
- Entrez gene IDs
- RefSeq IDs
- UniProt IDs
- Most other common gene identifier formats

### Input Formatting
- Genes can be separated by: newlines, spaces, tabs, or commas
- Case-insensitive
- Duplicates are automatically removed
- Special characters are automatically cleaned

### Background Genes
- Same format as query genes
- Should represent all genes measured in your experiment
- Maximum: 30,000 genes
- Optional but highly recommended for RNA-seq data

## Output Files

All results are downloadable in multiple formats:

### Tables (CSV format)
- Enrichment results (top pathways and all significant pathways)
- Gene conversion table with genomic locations
- Detailed gene information
- Promoter motif analysis results
- Gene groupings by functional category
- Network nodes and edges

### Figures (PDF, PNG, SVG)
- Enrichment charts
- Hierarchical clustering trees
- Gene characteristic plots
- Genome browser plots

### Interactive Files
- Network diagrams (HTML format)
- Links to STRING-DB interactive networks

## Related Tools

**iDEP** ([http://bioinformatics.sdstate.edu/idep/](http://bioinformatics.sdstate.edu/idep/))
Integrated Differential Expression and Pathway analysis - for comprehensive RNA-Seq data analysis including preprocessing, differential expression, and pathway analysis.

**Customized ShinyGO** ([http://bioinformatics.sdstate.edu/goc/](http://bioinformatics.sdstate.edu/goc/))
Includes custom genomes requested by users. To request a new species/genome, fill out [this form](https://forms.gle/zLtLnqxkW187AgT76).

## Testing

ShinyGO has a comprehensive test suite covering all modules and core functions.

```bash
Rscript run_tests.R                            # run all tests
Rscript run_tests.R test-mod_02_enrichment.R   # run a single file
```

Current status: **323 pass | 16 skip | 2 pre-existing fails** (no full database required for most tests).

See [`TESTING.md`](TESTING.md) for the full test inventory, helper/fixture documentation, and known behavioral notes. See [`tests/README.md`](tests/README.md) for the file structure and patterns for adding new tests.

## Citation

**Required Citation:**
Ge SX, Jung D & Yao R. ShinyGO: a graphical gene-set enrichment tool for animals and plants. *Bioinformatics* 36:2628–2629, 2020. [https://doi.org/10.1093/bioinformatics/btz931](https://doi.org/10.1093/bioinformatics/btz931)

**Note**: Simply including the URL is not sufficient - please cite the paper.

**Additional Citations** (if applicable):
- If using KEGG diagrams, cite [Pathview](https://doi.org/10.1093/bioinformatics/btt285) and [KEGG](https://doi.org/10.1093/nar/gkaa970)
- If using STRING-DB features, cite [STRING](https://string-db.org/)

## Team

ShinyGO is developed and maintained by a small team at [South Dakota State University (SDSU)](https://www.sdstate.edu/):

- **Dr. Xijin Ge** (PI) - [Faculty page](https://www.sdstate.edu/directory/xijin-ge)
- **Jianli Qi** (Research Associate)
- **Graduate Students**: Emma Spors, Ben Derenge

We share a passion for developing user-friendly tools accessible to all biologists, especially those without access to bioinformatics support.

## Contact and Support

**Email**: [gelabinfo@gmail.com](mailto:gelabinfo@gmail.com?Subject=ShinyGO)
**Twitter**: Follow [@StevenXGe](https://twitter.com/StevenXGe) for updates
**GitHub**: [https://github.com/gexijin/shinygo](https://github.com/gexijin/shinygo)
**Issues/Bugs**: File bug reports or feature requests on [GitHub](https://github.com/gexijin/shinygo)

For questions, suggestions, or to contribute data, please email the team.

## Documentation

**Detailed Demo**: [Supplemental PDF](https://www.biorxiv.org/content/biorxiv/suppl/2018/05/04/315150.DC1/315150-1.pdf)
**Original Paper**: [Bioinformatics, 2020](https://doi.org/10.1093/bioinformatics/btz931)
**Version History**: See [CHANGELOG.md](CHANGELOG.md)

## Previous Versions

Previous versions remain accessible for reproducibility:

- [ShinyGO V0.76](http://bioinformatics.sdstate.edu/go76/) - Ensembl Release 104 (archived Sept 2, 2022)
- [ShinyGO V0.75](http://bioinformatics.sdstate.edu/go75/) - Ensembl Release 104 (archived April 4, 2022)
- [ShinyGO V0.74](http://bioinformatics.sdstate.edu/go74/) - Ensembl Release 104 (archived Feb 8, 2022)
- [ShinyGO V0.65](http://bioinformatics.sdstate.edu/go65/) - Ensembl Release 103 (archived Oct 15, 2021)
- [ShinyGO V0.61](http://bioinformatics.sdstate.edu/go61/) - Ensembl Release 96 (archived May 23, 2020)

See [CHANGELOG.md](CHANGELOG.md) for complete version history.

## Acknowledgments

Thank you to all users who have supported ShinyGO through feedback, citations, and support letters. Your engagement helps ensure the continued development and hosting of this tool.

## License

This project is intended for academic and research use. For commercial use inquiries, please contact the team.
