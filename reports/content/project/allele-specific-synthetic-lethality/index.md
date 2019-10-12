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
Thus, a current goal of much of cancer research is focussed on finding synthetic lethal interactions with common oncogenes.

## Allele-specific synthetic lethality

The goal of this portion of the Comutation project is to identify synthetic lethal interactions with specific *KRAS* alleles.
From a theory point-of-view, these interactions would provide information on the differences between the alleles.
From a therapeutics point-of-view, the identification of allele-specific interactions could be used to further personalize treatment strategies.

## Data

For this analysis, we used the Cancer Dependency Mapping (DepMap) data.
In this ongoing study, the cell lines from the Cancer Cell Line Encyclopedia (CCLE) were analyzed in a genome-wide, CRISPR-Cas9, loss-of-function screen whereby each gene was systemically knocked-out and the loss of growth of the cell line was measured.
Thankfully, many omics data are available for the cell lines, including RNA expression, mutations, and CNA.
