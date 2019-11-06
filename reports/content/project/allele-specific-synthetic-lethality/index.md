---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "KRAS Allele-Specific Synthetic Lethality"
summary: ""
authors: [~]
tags: ["synthetic lethality", "KRAS"]
categories: []
date: 2019-10-09T06:58:30-04:00

# Optional external URL for project (replaces project detail page).
external_link: ""

# Featured image
# To use, add an image named `featured.jpg/png` to your page's folder.
# Focal points: Smart, Center, TopLeft, Top, TopRight, Left, Right, BottomLeft, Bottom, BottomRight.
image:
  caption: ""
  focal_point: ""
  preview_only: false

# Custom links (optional).
#   Uncomment and edit lines below to show custom links.
# links:
# - name: Follow
#   url: https://twitter.com
#   icon_pack: fab
#   icon: twitter

url_code: ""
url_pdf: ""
url_slides: ""
url_video: ""

# Slides (optional).
#   Associate this project with Markdown slides.
#   Simply enter your slide deck's filename without extension.
#   E.g. `slides = "example-slides"` references `content/slides/example-slides.md`.
#   Otherwise, set `slides = ""`.
slides: ""
---

## Synthetic Lethality

Synthetic lethality refers to a genetic interaction whereby one mutation sensitizes an organism to another normally-benign mutation.
This overlaps with "dependency" whereby the first mutation pushes the organism into a state of reliance on specific functions.
Blocking (by mutation or inhibition) this function would cause the mutant organism to die, but not effect the wild-type organisms.
Thus, much of current cancer therapy research is focused on finding synthetic lethal interactions with common oncogenes.

## Allele-specific synthetic lethality

The goal of this portion of the Comutation project is to identify synthetic lethal interactions with specific *KRAS* alleles.
From a theory point-of-view, these interactions would provide information on the differences between the alleles.
From a therapeutics point-of-view, the identification of allele-specific interactions could be used to further personalize treatment strategies towards the pursuit of precision medicine.

## Data

For this analysis, we used the Cancer Dependency Mapping (DepMap) data.
In this ongoing study, the cell lines from the Cancer Cell Line Encyclopedia (CCLE) were analyzed in a genome-wide, CRISPR-Cas9, loss-of-function screen whereby each gene was systemically knocked-out and the loss of growth of the cell line was measured.
Thankfully, many omics data are available for the cell lines, including RNA expression, mutations, and CNA.

(Of note, there was an RNAi screen, though it was not as extensive as the CRISPR screen, and thus, was not used for this analysis.)

## Results

### Identification synthetic lethal interactions

The first test employed searched for individual genes with allele-specific synthetic lethality.
Any depletion effects of knocking-out a gene that could be explained by the gene's own RNA expression were removed first, as this was likely the cause.
Then, an ANOVA was used to find genes that had differential effects across cell lines when grouped by their *KRAS* allele.
The genes with a statistically significant (p < 0.01) difference are shown below.
Below the heatmaps are specific genes with strong allele-specific effects.

**COAD**

![](/img/graphs/10_10_linear-modeling-syn-let_pheatmaps/COAD_CRISPR_pheatmap.svg)

![](/img/graphs/10_10_linear-modeling-syn-let_boxplots/COAD-IDH1.svg)

![](/img/graphs/10_10_linear-modeling-syn-let_boxplots/COAD-PIP5K1A.svg)

![](/img/graphs/10_10_linear-modeling-syn-let_boxplots/COAD-WDR26.svg)

**LUAD**

![](/img/graphs/10_10_linear-modeling-syn-let_pheatmaps/LUAD_CRISPR_pheatmap.svg)

Because there are too many genes to show their names for LUAD, I identified enriched functions in the clusters.
The functions that are a single gene are either transcription factors or kinases with their targets enriched in the cluster (the number at the end of each bar indicates the cluster).

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/functional-enrichment_LUAD.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-ARMC6.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-EPN2.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-FOXJ3.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-HSBP1.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-MAPK6.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-P4HTM.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-PTBP3.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-TP53BP1.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/LUAD-STAMBP.svg" width=300 />

**PAAD**

![](/img/graphs/10_10_linear-modeling-syn-let_pheatmaps/PAAD_CRISPR_pheatmap.svg)

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/PAAD-JUN.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/PAAD-EGLN2.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/PAAD-GPT2.svg" width=300 />

<img src="/img/graphs/10_13_linear-modeling-syn-let_fxnal-enrichment/PAAD-MRPS25.svg" width=300 />


### Functional enrichment of synthetic lethal interactions

I use GSEA to identify functions that were generally enriched by *KRAS* allele.
Select results are shown below, but the full results are available in the notebook page from [2019-10-18]({{< ref "/post/2019-10-18-2019-10-18.md" >}}).
Note that a **positive** NES corresponds to **stronger** depletion effect.



![](/img/graphs/10_37_gsea-depmap-analysis/gsea-results-COAD-select.svg)

G12V de-enriched for nonsense-mediated decay

![](/img/graphs/10_37_gsea-depmap-output/COAD_G12V/enplot_REACTOME_NONSENSE_MEDIATED_DECAY_NMD_66.svg)

G12V enriched in TP53-regulated metabolism

![](/img/graphs/10_37_gsea-depmap-output/COAD_G12V/enplot_REACTOME_TP53_REGULATES_METABOLIC_GENES_57.svg)

G13D de-enriched for complement

![](/img/graphs/10_37_gsea-depmap-output/COAD_G13D/enplot_REACTOME_COMPLEMENT_CASCADE_117.svg)

G13D de-enriched for the ETC

![](/img/graphs/10_37_gsea-depmap-output/COAD_G13D/enplot_REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_69.svg)


**LUAD**

![](/img/graphs/10_37_gsea-depmap-analysis/gsea-results-LUAD-select.svg)

G12C enriched in the Bard1 pathway

![](/img/graphs/10_37_gsea-depmap-output/LUAD_G12C/enplot_PID_BARD1_PATHWAY_33.svg)

G12C enriched in the Faconi anemia pathway

![](/img/graphs/10_37_gsea-depmap-output/LUAD_G12C/enplot_REACTOME_FANCONI_ANEMIA_PATHWAY_3.svg)

G12C de-enriched for $\beta$ alanine metabolism

![](/img/graphs/10_37_gsea-depmap-output/LUAD_G12C/enplot_KEGG_BETA_ALANINE_METABOLISM_66.svg)


**PAAD**

![](/img/graphs/10_37_gsea-depmap-analysis/gsea-results-PAAD-select.svg)

G12R enriched in JNK (c-Jun kinases) phosphorylation and activation mediated by activated human TAK1

![](/img/graphs/10_37_gsea-depmap-output/PAAD_G12R/enplot_REACTOME_JNK_C_JUN_KINASES_PHOSPHORYLATION_AND_ACTIVATION_MEDIATED_BY_ACTIVATED_HUMAN_TAK1_57.svg)

G12R enriched in NF-kB Signaling Pathway

![](/img/graphs/10_37_gsea-depmap-output/PAAD_G12R/enplot_BIOCARTA_NFKB_PATHWAY_45.svg)

G12R enriched in G2/M DNA damage checkpoint

![](/img/graphs/10_37_gsea-depmap-output/PAAD_G12R/enplot_REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT_24.svg)

G12V de-enriched in the hedgehog signaling pathway

![](/img/graphs/10_37_gsea-depmap-output/PAAD_G12V/enplot_HALLMARK_HEDGEHOG_SIGNALING_72.svg)

G12V de-enriched in the regulation of cholesterol biosynthesis by SREBP (SREBF)

![](/img/graphs/10_37_gsea-depmap-output/PAAD_G12V/enplot_REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF_69.svg)


### No obvious differences in the general distribution of depletion effects per allele

Of note, there are no obvious differences in the overall distribution of depletion effects per *KRAS* allele.

![](/img/graphs/10_45_overall-dependencies-by-allele/dependency-overview_COAD.svg)

![](/img/graphs/10_45_overall-dependencies-by-allele/dependency-overview_LUAD.svg)

![](/img/graphs/10_45_overall-dependencies-by-allele/dependency-overview_PAAD.svg)

