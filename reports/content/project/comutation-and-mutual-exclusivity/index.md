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

I also created [UpSet plots](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4720993/) to show the size and overlap of the sets.
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

The plot below shows the distribution of number of mutually exclusive events by p-value of the association.

![](/img/graphs/20_35_rc-fisher-comparison/rc_mutations_distribition.svg)

The rest of the summary is subdivided by cancer.

[COAD](#coad)  
[LUAD](#luad)  
[MM](#mm)  
[PAAD](#paad)   


### COAD

Below are the high-level visualizations of the genetic interactions.
The alleles tend to be connected by mutual exclusive interactions (instead of comutation interactions).
Each allele also has their own set of comutation and mutually exclusive interactors.

![](/img/graphs/20_40_highlivel-genetic-interactions/genetic_interaction_network_COAD_thick.svg)

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

A more direct visualization of the mutually exclusive and comutation interactions is via an oncoplot.
On a heatmap, it shows a sample per column and a gene per row.
The cell is filled with a color if the sample is mutated in the gene.
The arrangement of the samples highlights the genetic interactions by showing either the overlap or exclusivity.

I have shown them here for a select few of the alleles with interesting interactions.
For the paper, a more selective portion of the interactions will be shown.

![](/img/graphs/20_50_rainfall-plots/COAD_A146T_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/COAD_G12D_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/COAD_G13D_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/COAD_G13D_exclusivity_oncostrip_allInteractors.svg)

---

### LUAD

Below are the high-level visualizations of the genetic interactions.
The alleles tend to be connected by mutual exclusive interactions (instead of comutation interactions).
This is especially true for G12C, G12V, and G12D.
Each allele also has their own set of comutation and mutually exclusive interactors.
There are many genes that are mutually exclusive with one allele and comutate with another.

![](/img/graphs/20_40_highlivel-genetic-interactions/genetic_interaction_network_LUAD_thick.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/LUAD_comutation_upset.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/LUAD_exclusivity_upset.svg)


The genetic interaction networks with genes in the *a prior* lists are divided in to the KEGG pathways (top) and the CGC (bottom).
There are very strong mutually exclusive relationships between *EGFR*, *BRAF*, and *TP53*, as expected.
Also, there is a strong comutation interaction between G12C and *STK11*.
In the CGC network, it appears that G12C is mutually exclusive with *CTNNB1* ($\beta$-catenin).
G12C also is mutually exclusive with three different genes involved in voltage-gated calcium channels.
G12C is mutually exclusive with *JAG2* and comutates with *NOTCH1* suggesting a special relationship with the Notch signaling pathway.
The three alleles linked to the smoking mutational signature, G12A/C/V, and Q61L have strong comutation interactions with ATM.
There are several comutation interactions that are unexpected, namely between G12D and G12R with *NKX2-1*, Homeobox protein Nkx-2.1.
This gene is a "transcription factor that binds and activates the promoter of thyroid specific genes," and "may play a role in lung development and surfactant homeostasis," ([UniProt P43699](https://www.uniprot.org/uniprot/P43699)).
Another interesting relationship is that *KRAP1* comutates with G13D and G13C, but is mutually exclusive with G12V and G12C.
This could indicate an codon-specific property of hyperactivated *KRAS* in LUAD.

![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_LUAD_kegg.svg)

![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_LUAD_cgc.svg)

The results of an enrichment analysis using a gene set of PPI hub proteins are shown first.
The physical interactors with $\beta$*-catenin* and SMARC4 are enriched in the interactors of G12C. 
The interactors of PKC-$\alpha$ (gene *PKRCA*) are enriched in G12V genetic interactions.
PKC-$\alpha$ is a "calcium-activated, phospholipid- and diacylglycerol (DAG)-dependent serine/threonine-protein kinase that is involved in positive and negative regulation of cell proliferation, apoptosis, differentiation, migration and adhesion, tumorigenesis, cardiac hypertrophy, angiogenesis, platelet function and inflammation, by directly phosphorylating targets such as RAF1, BCL2, CSPG4, TNNT2/CTNT, or activating signaling cascade involving MAPK1/3 (ERK1/2) and RAP1GAP," ([UniProt P17252](https://www.uniprot.org/uniprot/P17252)).

![](/img/graphs/20_45_fxnal-enrich-genetic-interactions/enrichr_LUAD_PPI_Hub_Proteins.svg)

The results of an enrichment analysis using a gene set of GO Biological Processes are shown below.
G12A seems to have an enrichment for genes involved in P53-mediated cell cycle arrest.
There are several Ras-related processes enriched in the G12C interactors.

![](/img/graphs/20_45_fxnal-enrich-genetic-interactions/enrichr_LUAD_GO_Biological_Process_2018.svg)

Below are the oncoplots for the more common alleles in LUAD.

![](/img/graphs/20_50_rainfall-plots/LUAD_G12D_exclusivity_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/LUAD_G12D_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/LUAD_G12C_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/LUAD_G12C_exclusivity_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/LUAD_G12V_exclusivity_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/LUAD_G12V_comutation_oncostrip_allInteractors.svg)


---

### MM

Below are the high-level visualizations of the genetic interactions.
There were not very many hits for multiple myeloma and very little overlap of the interacting genes between *KRAS* alleles.

![](/img/graphs/20_40_highlivel-genetic-interactions/genetic_interaction_network_MM_thick.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/MM_comutation_upset.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/MM_exclusivity_upset.svg)

Below is the network of genetic interactions for genes in the *a priori* lists.
There is a detectable mutually exclusive interaction between Q61H and *BRAF*, though it is documented that MM can be highly multi-clonal.
Therefore, there are often subpopulations of *KRAS* mutants and *BRAF* mutants within a single disease.
This is likely true for KRAS and other genes, making the interpretation of these results difficult. 

![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_MM_allLists.svg)

There was little evidence for functional enrichment in these genes.

---

### PAAD

Below are the high-level visualizations of the genetic interactions.
There were many allele-specific comutation and mutually exclusive interactions detected in PAAD.
Further or more strict thresholds may be applied in future analyses to reduce false-positives.

![](/img/graphs/20_40_highlivel-genetic-interactions/genetic_interaction_network_PAAD_thick.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/PAAD_comutation_upset.svg)

![](/img/graphs/20_40_highlivel-genetic-interactions/PAAD_exclusivity_upset.svg)


Below are the interaction networks for genes in the *a priori* gene lists.

The first is of genes found to directly interact with *KRAS* in the CGC.
Along with *BRAF*, there were some unexpected genes found to be mutually exclusive with multiple alleles.
Interestingly, *RNF43* was foound to be mutually exclusive with G12R, but comutates with G12D. *RNF43* is "E3 ubiquitin-protein ligase that acts as a **negative regulator of the Wnt signaling pathway** by mediating the ... degradation of ... Frizzled," ([UniProt Q68DV7](https://www.uniprot.org/uniprot/Q68DV7)).
In addition, *TP53* is mutually exclusive with G12V, but comutates with Q61H and G12D.
This is interesting because of the general perception that both KRAS and TP53 mutations are quintessential steps of PAAD progression.
*SMARCA4* is strongly mutually exclusive with G12R.


![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_PAAD_cgc.svg)

The following is the network of genetic interactions with genes found to directly interact with *KRAS*.
One of the strongest interactions is between G12V and *ATP2B1*, the protein product of which helps regulate intracellular calcium levels at the plasma membrane.
Another strong interaction is the mutually exclusive relationship between G12D and *NF1*.
The obvious mechanism to propose for this is the same as that for *BRAF*.


![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_PAAD_BioID.svg)

The next genetic network is of the genes from select KEGG pathways.
An interaction that really pops out is the strong comutation between Q61K and *MAP3K13* and compared to the mutually exclusive relationships between G12V and G12R and the same gene.
A similar case of conflicting interactions is found with *MAP2K4* where it comutates with Q61R and is mutually exclusive with G12D.
G12D strongly comutates with *TGFBR2* - this transmembrane S/T-kinase activates SMAD2/4 when it binds with the cytokines TGFB1/2/3 ([UniProt P37173](https://www.uniprot.org/uniprot/P37173)).
*JAG2* is found to comutate with G12R, yet is mutually exclusive with G12V.

![](/img/graphs/20_43_apriori-lists-genetic-interactions/goi_overlap_genetic_interactions_network_PAAD_kegg.svg)

Two gene sets with interesting results from the enrichment analysis are shown below.
The top is of kinases whose substrates are overrepresented in the associations, and the bottom is of Panther pathways.

There are a lot of substrates for the kinase PKA C-$\alpha$ (gene *PRKACA*) that have genetic interactions with G12D and G12V.
This kinase is involved in many cellular functions including Hedgehog and NF$\kappa$B signaling.
The substrates of Akt1 are also overrepresented in the genes with genetic interactions with G12D.

![](/img/graphs/20_45_fxnal-enrich-genetic-interactions/enrichr_PAAD_KEA_2015.svg)

Many genes involved in Wnt signaling have genetic interactions with several alleles.
Genes of the EGFR pathway have many interactions with G12D.

![](/img/graphs/20_45_fxnal-enrich-genetic-interactions/enrichr_PAAD_Panther_2016.svg)

Below are the oncoplots for the more interesting alleles in PAAD.

![](/img/graphs/20_50_rainfall-plots/PAAD_G12D_bothInteractions_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/PAAD_G12D_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/PAAD_G12R_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/PAAD_G12R_exclusivity_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/PAAD_Q61H_comutation_oncostrip_allInteractors.svg)

![](/img/graphs/20_50_rainfall-plots/PAAD_Q61H_exclusivity_oncostrip_allInteractors.svg)
