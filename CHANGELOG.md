# ShinyGO Version History and Changelog

This document contains the complete version history of ShinyGO, including all updates, bug fixes, and database changes.

## Version 0.80 (In Testing - May 1, 2023)

**Status**: Beta release at [http://bioinformatics.sdstate.edu/go80/](http://bioinformatics.sdstate.edu/go80/)

**Major Updates**:
- Updated to Ensembl Release 107
- Expanded species coverage to **620 species total**:
  - 215 main species
  - 177 metazoa species
  - 124 plant species
  - 33 protists
  - 1 bacteria
- Updated STRING-DB to version **11.5** with **14,094 species**

**Credits**: Thanks to Jenny's hard work on this major update.

---

## Version 0.77 (Current Production)

**Status**: Production release at [http://bioinformatics.sdstate.edu/go](http://bioinformatics.sdstate.edu/go)

**Updates**: Minor stability and performance improvements.

---

## Version 0.76.3 (October 26, 2022)

**New Features**:
- Added hover text/tooltips to improve user experience
- Changed plot styles for better readability

**Behavior Changes**:
- When users select "Sort by Fold Enrichment", the minimum pathway size is automatically raised to **10 genes** (from 2) to filter out noise from tiny gene sets

**Database**: Ensembl Release 104 (revised)

---

## Version 0.76.2 (September 28, 2022)

**Important Changes**:
- **KEGG is now the default pathway database**
- **Reverted to v0.76 gene counting method**: All protein-coding genes are used as the background by default (instead of pathway database genes)
- The v0.76.1 feature (using pathway database for gene counts) is now available as an **optional checkbox**: "Use pathway database for gene counts"

**Rationale**: Based on user feedback that when using smaller pathway databases like KEGG, the v0.76.1 method changed P-values substantially.

**Database**: Ensembl Release 104 (revised)

---

## Version 0.76.1 (September 3, 2022)

**New Feature**:
- Improved gene counting method for P-value calculation
- A gene must now match at least one pathway in the selected pathway database
- Genes not matching any pathway are ignored in hypergeometric distribution calculations
- Applies to both query and background genes

**Database**: Ensembl Release 104 (revised)

**Note**: This change was later made optional in v0.76.2 due to user feedback.

---

## Version 0.76 (April 19, 2022)

**Major Improvements**:
- Enhanced pathway filtering options
- Improved pathway sorting capabilities
- Better figure downloading functionality with multiple format options

**Download Options**:
- Added support for PDF, SVG, and high-resolution PNG downloads (added April 17, 2022)

**New Features** (added April 8, 2022):
- Option to remove redundant pathways (pathways sharing 95% of genes)
- Filters to remove extremely large or small pathways
- Changed interface to always show KEGG tab (previously hidden when not selected)

**Bug Fixes**:
- March 7, 2022: Fixed R library issue that affected KEGG diagrams for some organisms

**Database**: Ensembl Release 104 (revised)

**Previous Version**: [ShinyGO V0.75](http://bioinformatics.sdstate.edu/go75/) archived

---

## Version 0.75 (February 8, 2022 - Official Release)

**Critical Updates**:
- **R upgraded from 4.05 to 4.1.2** - Resolved STRING API issues
- Updated several Bioconductor packages
- Database updated to Ensembl Release 104
- STRING-DB updated to version 11

**Bug Fixes**:
- **February 26, 2022**: Fixed critical bug in Plot tab when background genes are used
  - Background genes were not correctly used to calculate distributions of various gene characteristics
  - **Action Required**: If these plots were important in your study, please re-analyze your genes

**Performance**:
- Much faster processing even with large sets of background genes
- Now recommends the use of background genes in enrichment analysis

**Testing Period**: Database updates available in testing mode since November 15, 2021

**Previous Version**: [ShinyGO V0.74](http://bioinformatics.sdstate.edu/go74/) archived

---

## Version 0.74 (February 8, 2022)

**Database Updates**:
- Ensembl Release 104
- STRING-DB v11

**New Features**:
- Added species from Ensembl Fungi
- Added species from Ensembl Protists

**Performance**: Improved speed with background genes

**Archive Date**: February 8, 2022

---

## Version 0.741 (October 23, 2021)

**New Features**:
- **Fully customizable enrichment chart**
  - Switch between bar, dot, or lollipop plots
  - Adjustable colors, fonts, sizes
  - Multiple theme options

**Gene Information**:
- Detailed gene information with external links on the Genes tab
- Enhanced gene annotation display

---

## Version 0.66 (June 6, 2021)

**Interface**:
- Adjusted and improved user interface

---

## Version 0.66 (June 2, 2021)

**New Feature**:
- Added support for **customized background genes**
- Users can now upload their own background gene sets for more accurate enrichment analysis

---

## Version 0.65 (October 25, 2021)

**New Features**:
- **Interactive genome plot** using Plotly
- **Identification of genomic regions** significantly enriched with user genes
- Sliding window analysis across chromosomes

**Archive Date**: October 15, 2021

---

## Version 0.65 (May 23, 2021)

**Database Updates**:
- Updated to Ensembl Release 103
- Updated to STRING-DB v11

**Archive Date**: October 15, 2021

---

## Version 0.61 (November 3, 2019)

**New Features**:
- Improved visualizations based on peer reviewer suggestions
- **Interactive networks** using visNetwork
- Enhanced user experience

**Archive Date**: May 23, 2020

---

## Version 0.60 (May 20, 2019)

**Database Updates**:
- Upgraded to Ensembl BioMart Release 96
- Added annotation from STRING-DB v10

**Archive Date**: November 6, 2019

---

## Version 0.51 (March 29, 2019)

**Database Updates**:
- Updated annotation to Ensembl Release 95

**Interface Changes**:
- Improved user interface
- Added demo gene lists
- Better error messages

**Archive Date**: May 20, 2019

---

## Version 0.50 (September 10, 2018)

**Database Updates**:
- Upgraded to Ensembl BioMart Release 92

**Archive Date**: March 29, 2019

---

## Version 0.42 (April 30, 2018)

**Updates**:
- Changed figure configurations for hierarchical clustering tree
- Improved tree visualization parameters

---

## Version 0.41 (April 27, 2018)

**Technical Updates**:
- Switched to ggplot2 for all visualizations
- Added grid and gridExtra packages for better plot layouts

**Archive Date**: July 11, 2018

**Database**: Ensembl Release 91

---

## Version 0.4 (April 24, 2018)

**Major New Features**:
- **STRING API integration** for protein-protein interaction networks
- **KEGG pathway diagrams** with user genes highlighted
- **Hierarchical clustering tree** visualization
- **Network visualization** for enriched pathways

This was a major milestone release that added many of the visualization features that define ShinyGO today.

---

## Critical Bug Fixes (Post-Release)

### January 19, 2023 - Serious Bug Fix (v0.76)

**Issue**: Some genes are represented by multiple gene IDs in Ensembl, causing them to be counted more than once in enrichment calculations.

**Impact**: Users who pasted Ensembl gene IDs to ShinyGO 0.76 between **April 4, 2022 and January 19, 2023** should **rerun their analysis**.

**Status**: Believed to be fixed.

**Recommendation**: Always double-check results with other tools such as:
- G:profiler
- Enrichr
- STRING-DB
- DAVID

**Note**: ShinyGO has not been thoroughly tested. Please always verify results with complementary tools.

---

## Database Version Reference

| ShinyGO Version | Ensembl Release | STRING-DB Version | Archive Date |
|-----------------|-----------------|-------------------|--------------|
| 0.80 (beta) | 107 | 11.5 | In testing |
| 0.77 (current) | 104 (revised) | 11.5 | Current |
| 0.76.3 | 104 (revised) | 11.0 | Sept 2, 2022 |
| 0.76.2 | 104 (revised) | 11.0 | Sept 28, 2022 |
| 0.76.1 | 104 (revised) | 11.0 | Sept 3, 2022 |
| 0.76 | 104 (revised) | 11.0 | April 4, 2022 |
| 0.75 | 104 | 11.0 | Feb 8, 2022 |
| 0.74 | 104 | 11.0 | Feb 8, 2022 |
| 0.65 | 103 | 11.0 | Oct 15, 2021 |
| 0.61 | 96 | 10.0 | May 23, 2020 |
| 0.60 | 96 | 10.0 | Nov 6, 2019 |
| 0.51 | 95 | - | May 20, 2019 |
| 0.50 | 92 | - | March 29, 2019 |
| 0.41 | 91 | - | July 11, 2018 |

---

## Accessing Previous Versions

All previous versions remain functional for reproducibility purposes:

- **V0.76**: [http://bioinformatics.sdstate.edu/go76/](http://bioinformatics.sdstate.edu/go76/)
- **V0.75**: [http://bioinformatics.sdstate.edu/go75/](http://bioinformatics.sdstate.edu/go75/)
- **V0.74**: [http://bioinformatics.sdstate.edu/go74/](http://bioinformatics.sdstate.edu/go74/)
- **V0.65**: [http://bioinformatics.sdstate.edu/go65/](http://bioinformatics.sdstate.edu/go65/)
- **V0.61**: [http://bioinformatics.sdstate.edu/go61/](http://bioinformatics.sdstate.edu/go61/)
- **V0.60**: [http://bioinformatics.sdstate.edu/go60/](http://bioinformatics.sdstate.edu/go60/)
- **V0.51**: [http://bioinformatics.sdstate.edu/go51/](http://bioinformatics.sdstate.edu/go51/)
- **V0.50**: [http://bioinformatics.sdstate.edu/go50/](http://bioinformatics.sdstate.edu/go50/)
- **V0.41**: [http://bioinformatics.sdstate.edu/go41/](http://bioinformatics.sdstate.edu/go41/)

---

## Important Notes

### Testing Status
ShinyGO has not been thoroughly tested. Users should always:
- Double-check results with other enrichment tools
- Validate significant findings
- Report any bugs or unexpected behavior to the development team

### Reproducibility
We maintain access to all previous versions to ensure scientific reproducibility. If you used a specific version in your publication, that version will remain accessible.

### Database Updates
- Ensembl databases are typically updated 1-2 times per year
- STRING-DB updates follow their release schedule
- Curated pathway databases (KEGG, Reactome, etc.) are updated periodically
- Some lag time is expected between external database releases and ShinyGO updates

---

## Future Development

The ShinyGO team continues to:
- Expand species coverage
- Update databases regularly
- Add new visualization features
- Improve performance and usability
- Fix bugs based on user feedback

For the latest updates, follow [@StevenXGe](https://twitter.com/StevenXGe) on Twitter.

---

## Reporting Issues

If you encounter bugs or have feature requests:
- **Email**: [gelabinfo@gmail.com](mailto:gelabinfo@gmail.com?Subject=ShinyGO)
- **GitHub Issues**: [https://github.com/gexijin/shinygo/issues](https://github.com/gexijin/shinygo/issues)

---

*Last Updated: Based on code review February 2026*
