# To-Do

## General project tasks

1. Update time-line
2. Mock-up figures
3. Write


## Urgent


## Doc

* add logs and test output to website


## Analyses

- Change colors of GSEA to fit standards for dependency analysis.
    + Also make ranking plot like in "Defining DepMap" paper figure...
- Analyze comutation + dependency subnetworks:
    + include genes that interact with a few hits
    + identify specific motifs of *a priori* interest
    + statistics on submodules and clusters
- Lollipop plots for specific genes across multiple alleles
    + My API is in place and I just need to feed it genes to run.
    + I need to switch it over to using `cancer_full_coding_muts_df`.
- Make my own "lollipop" plot for KRAS mutation frequency


- Rethink Figure  1
    + It may be a better choice to merge it into Figure 2
        * I am expecting to have to completely re-do Figure 2 anyways, and I know where the data for it is.
- overlap between hits from comutation and synthetic lethal
- functional enrichment from genetic interactions
- are any of the syn. let. hits classified as pan-nonessential or pan-essential?
    + do any of them *become* lethal in an allele?
- Expected vs. observed frequency -> save to color-scaled table


## Questions

---

## Thoughts

Looking at the PPI subnetworks pulled out from the overlap of comutation and mutual exclusivity, there seem to be some striking clusters when looking for patterns in shape and color. 
Perhaps these are worth looking into as a whole for finishing up the paper.
For instance.


## Final global changes

- Change all `cache()` to `ProjectTemplate::cache()`
- A single script or command line to run the entire project.