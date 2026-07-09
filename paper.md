---
title: 'VizMut-lolliplot: A command-line pipeline for integrated genomic and protein variant visualization with population frequency context'
tags:
  - R
  - bioinformatics
  - variant visualization
  - lolliplot
  - genomics
  - population genetics
  - ALS
  - neurodegenerative diseases
authors:
  - name: Carlos Gómez-Pérez
    orcid: 0009-0006-5860-3121
    affiliation: 1
affiliations:
  - name: Department of Genetics, Instituto Nacional de Neurología y Neurocirugía Manuel Velasco Suárez, Mexico City, Mexico
    index: 1
date: 2026-07-08
bibliography: paper.bib
archive_doi: 10.5281/zenodo.21285302
repository-code: https://github.com/carlosgomezperez21/VizMut-lolliplot
---

# Summary

VizMut-lolliplot is an open-source command-line pipeline implemented in R for visualizing genetic variants as lolliplots at both the protein and genomic levels. The tool accepts a minimal input of variant identifiers in HGVS nomenclature and automatically retrieves genomic coordinates, ACMG pathogenicity classifications from ClinVar [@clinvar], protein domain annotations from UniProt, population allele frequencies from gnomAD v4 [@gnomad], and canonical transcript structures from NCBI using the MANE Select standard [@mane]. VizMut-lolliplot produces publication-ready figures integrating variant position, clinical significance, and population frequency context in a single visualization, addressing a critical need for researchers working with underrepresented populations in genomic studies.

# Statement of Need

Variant visualization tools are essential for interpreting genetic findings in clinical and research contexts. Existing tools such as lollipops [@lollipops], MutationMapper [@mutationmapper], and ProteinPaint [@proteinpaint] provide protein-level visualization but lack integration of genomic structure, UTR regions, population allele frequencies, and automatic data enrichment from public databases (Table 1). This fragmentation forces researchers to use multiple tools and perform manual data curation, increasing the risk of errors and limiting reproducibility.

This gap is particularly pronounced for researchers working with Latin American and admixed populations, where population-specific variant frequencies are critical for accurate pathogenicity interpretation. VizMut-lolliplot addresses these limitations by integrating in a single pipeline: (1) protein and genomic lolliplot visualization, (2) automatic enrichment from ClinVar, NCBI, and gnomAD v4, (3) population allele frequency pie charts across five continental populations (AFR, AMR, EAS, EUR, SAS), (4) UTR visualization, (5) MANE Select transcript support, and (6) variant filtering by ACMG classification, variant type, and phenotype.

**Table 1.** Feature comparison of VizMut-lolliplot with existing variant visualization tools. ✓ = fully supported; partial = limited support; — = not available; web = only available through web interface.

| Feature | lollipops | MutationMapper | ProteinPaint | VizMut-lolliplot |
|---|---|---|---|---|
| Protein lolliplot | ✓ | ✓ | ✓ | ✓ |
| Genomic structure plot | — | — | partial | ✓ |
| Multi-gene plot | — | — | — | ✓ |
| Chromosomal ideogram | — | — | — | ✓ |
| UTR visualization | — | — | — | ✓ |
| ACMG classification | — | partial | — | ✓ |
| gnomAD v4 frequencies | — | — | — | ✓ |
| Population pie charts | — | — | — | ✓ |
| UniProt features (auto) | — | partial | partial | ✓ |
| MANE Select transcript | — | — | — | ✓ |
| Enrichment from minimal input | — | — | — | ✓ |
| Liftover hg19 → hg38 | — | — | — | ✓ |
| Variant filtering | — | — | — | ✓ |
| Run log | — | — | — | ✓ |
| Command-line interface | ✓ | web | web | ✓ |
| Free / open source | ✓ | ✓ | ✓ | ✓ |

# Functionality

VizMut-lolliplot is invoked from the command line via `Rscript main.R` and supports three primary plot types (Figure 1):

**Protein lolliplot** (`--plot_type protein`): Variants are mapped to the protein sequence with functional features displayed in layers — domains, motifs, post-translational modifications (PTMs), and zinc fingers — retrieved automatically from UniProt using `--fetch_features TRUE` (Figure 2A).

**Single-gene lolliplot** (`--plot_type single_gene`): Variants are plotted over the genomic structure of a RefSeq transcript, including exons, introns, 5'UTR and 3'UTR regions. A chromosomal ideogram is displayed at the top of the figure. Population allele frequencies from gnomAD v4 are shown as proportional pie charts on each lollipop head when `--gnomad TRUE` is specified (Figure 2B).

**Multi-gene lolliplot** (`--plot_type multi_gene`): Multiple genes are displayed in stacked panels, each with its own ideogram and genomic structure.

The enrichment mode (`--enrich TRUE`) accepts a minimal two-column CSV (variant identifier and HGVS cDNA notation) and automatically retrieves coordinates, ACMG classification, phenotype, rsID, and population allele frequencies. Liftover from hg19 to hg38 is performed transparently using rtracklayer [@rtracklayer]. The `--show_only` flag allows filtering of variants by ACMG classification, variant type, or phenotype association before plotting.

![VizMut-lolliplot pipeline. Overview of the data flow from minimal CSV input through automatic enrichment from public databases to publication-ready plots.](pipeline_diagram.png)

**Figure 1.** VizMut-lolliplot pipeline overview.

![Example outputs. (A) Protein lolliplot of KMT2B with automatically retrieved UniProt domains, motifs, PTMs and zinc fingers. (B) Single-gene lolliplot of GJB1 showing exon structure, UTR regions and gnomAD v4 population frequency pie charts.](figures_panel.png)

**Figure 2.** Example outputs generated by VizMut-lolliplot.

# Case Study

VizMut-lolliplot was applied to visualize SOD1 variants identified in a cohort of Mexican ALS patients. Previously reported variants were plotted with gnomAD v4 population frequencies, demonstrating the predominance of European allele frequencies and the near-absence of pathogenic variants in non-European populations. Five novel variants identified for the first time in Mexican ALS patients (Gómez-Pérez et al., in preparation) were visualized separately; none are reported in gnomAD v4 or ClinVar, consistent with population-specific variants in admixed Latin American individuals.

# Performance

Benchmarking of the enrichment mode using GJB1 variants (n = 5–68) shows a mean processing time of 5 seconds for 5 variants and 22 seconds for 68 variants, with a plateau pattern at larger input sizes suggesting that network latency from external API calls is the primary bottleneck rather than computational complexity.

# Acknowledgements

The author thanks Dra. Petra Yescas Gómez and the Department of Genetics at INNNMVS for institutional support, and the treating physicians who provided clinical data for the ALS cohort used in the case study.

# References
