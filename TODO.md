# To-Do

## General project tasks

1. Update time-line
2. Mock-up figures
3. Write


## Urgent


## Doc

* add logs and test output to website


## Analyses

- Use multiple PPI networks for the analysis of the overlapping interactions.
    + BioPlex2, INTACT, HINT, STRING
    + OP suggests including other genes that interact with multiple of the hits
- Analyze comutation + dependency subnetworks:
    + include genes that interact with a few hits
    + identify specific motifs of *a priori* interest
    + statistics on submodules and clusters

- Rethink Figure  1
    + It may be a better choice to merge it into Figure 2
        * I am expecting to have to completely re-do Figure 2 anyways, and I know where the data for it is.
- overlap between hits from comutation and synthetic lethal
- look into specific genes with genetic interactions
    + what other plots are available from 'maftools'?
    + TP53 type of mutations for claims of reduced comutation (are the mutations that do occur considered the hotspots?)
- functional enrichment from genetic interactions
- are any of the syn. let. hits classified as pan-nonessential or pan-essential?
    + do any of them *become* lethal in an allele?
- Look into Synthetic Lethal Consortium for list of KRAS synthetic lethals
    + Tim's Cell Reports paper may be a good source to cross reference.
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