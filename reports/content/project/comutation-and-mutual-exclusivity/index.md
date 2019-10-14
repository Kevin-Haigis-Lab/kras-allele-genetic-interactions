---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Comutation and Mutual Exclusivity"
summary: "Testing for genetic interactions in human tumor samples based on their mutational patterns."
authors: [~]
tags: ["comutation", "KRAS"]
categories: []
date: 2019-10-09T06:58:09-04:00

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

## Overview

**Our hypothesis is that genes with comutate or be mutually exclusive with specific *KRAS* alleles depending on the biological properties of the allele.**
Comutation would occur when the additional mutation enhances the growth phenotype of the tumor.
And example of this would be *KRAS*, *TP53*, and *CDKN2A* in PAAD.
On the other hand, another gene and *KRAS* would mutate in a mutually exclusive pattern because doing so leads to a reduction in growth.
A very strong example of this is the mutually exclusive activation mutations in *BRAF* and *KRAS*.


## Analysis

### Comutation

To study comutation, we conducted a one-sided Fisher's exact test for association between each *KRAS* allele and all other genes.

### Mutual Exclusivity

The Row-Column test for mutual exclusivity (the "RC-test") was developed by the [Raphael Lab](http://compbio.cs.brown.edu/projects/wext/) in 2016 and described in [*A weighted exact test for mutually exclusive mutations in cancer* (Leiserson *et al.*, 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27587696).
The software supplied with the paper searched for any groups of genes demonstrating mutually exclusive relationships, but took far too long to compute over all of our samples.
Therefore, I implemented the RC-test in my own R package (in dev, [WeXT](https://github.com/jhrcook/wext)) and used [Snakemake](https://snakemake.readthedocs.io/en/stable/) to orchestrate the runs on O2 in parallel.
Again, the test was run between each *KRAS* allele and all other genes.

### High-level visualization

To provide a higher-level visualization of the results, the genetic interactions were displayed as a network.
These are meant to be interpreted as "clouds" where the density of the cloud provides an idea of the number of interactions.

I also created [UpSet plots]() to show the size and overlap of the sets.
These can be interpreted as high-dimensional Venn diagrams.

### Analysis of biological properties of the genetic interactors

I used *a priori* gene lists to identify any genes of interest that have detectable genetic interactions with the *KRAS* alleles.
The lists were notable KEGG pathways, the physical interactors of *KRAS* from a BioID experiment conducted in [Kovalski *et al.*, 2019](https://www.ncbi.nlm.nih.gov/pubmed/?term=30639242), and the [COSMIC Cancer Gene Consensus](https://cancer.sanger.ac.uk/census).

I also used Enrichr ([Chen *et al.*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23586463)) to identify gene sets enriched in the interactors for each *KRAS* allele.



## Results

The following figure shows the number of genes with a genetic interaction that was detected for each *KRAS* allele.
For comutation, there had to be at least 3 comutation events to be included.
For mutual exclusivity, the other gene had to be mutated at least 10 times in the cancer.
For both tests, a p-value of 0.01 was used for the declaration of statistical significance.

The plot shows that the number of detectable interactions increases with the number of samples, as expected.

![](/img/graphs/20_35_rc-fisher-comparison/rc_fisher_comparison_specific.svg)

The rest of the summary is subdivided by cancer.


### COAD

Below are the high-level visualizations of the genetic interactions.
The alleles tend to be connected by mutual exclusive interactions (instead of comutation interactions).
Each allele also has their own set of comutation and mutually exclusive interactors.

![](/img/graphs/20_40_highlivel-genetic-interactions/genetic_interaction_network_COAD.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/COAD_comutation_upset.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/COAD_exclusivity_upset.svg)


The network below shows all of the interactions with genes in the *a priori* lists.
As expected, *BRAF* is strongly mutually exclusive with many of the alleles and *PIK3CA* and *APC* comutate with several alleles, too.
Surprisingly, G12D and G13D are mutually exclusive with *TP53*.
A146T has two interesting comutation interactions with *PORCN* and *MAP3K1*.
G12V has a very strong comutation interaction with *SMAD4*.
G12C has a strong comutation interaction with *RAPGEF2*.
G12S (bottom left) has a comutation interaction with *SMAD3*.

![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_COAD_allLists.svg)

The KEGG pathway gene set was very informative.
Notably, G13D was enriched for the mTOR pathway, including *RICTOR* as shown in the above network.
The genes with interactions with G12D are enriched in adhesion pathways ("Regulation of actin cytoskeleton" and "Focal adhesion", "Leukocyte TEM"), possibly linked to the enrichment for the Rap1 signaling pathway.
Many of the alleles have interactions with genes involved in ERBB signaling.

![](/img/graphs/20_45_fxnal-enrich-genetic-interactions/enrichr_COAD_KEGG_2019_Human.svg)

Also of note, but not shown graphically here, the physical interactors of SMAD2, a transcription factor, were found to be enriched in the G12D genetic interactors by two different gene sets.

---

### LUAD