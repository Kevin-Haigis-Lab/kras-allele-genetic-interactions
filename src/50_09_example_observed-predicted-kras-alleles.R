# This is a small-scale example of predicting the KRAS allele frequency from
# the mutational signature. It takes a single sample, a few fake mutational
# signatures and artificial contributions to a few fake KRAS alleles.

#### ---- Example analysis with a single tumor sample ---- ####


#>  mutational signature x tumor sample
mutsig_by_tumor <- matrix(c(0.3, 0.2, 0.5), ncol = 1)
colnames(mutsig_by_tumor) <- "T"
rownames(mutsig_by_tumor) <- paste0("sig", 1:3)
mutsig_by_tumor
# Each column should sum to 1.
all(colSums(mutsig_by_tumor) == 1)


#> mutational signature x KRAS allele
mutsig_by_kras <- matrix(
    c(
        0.1, 0.1, 0.1, 0.7,
        0.5, 0.1, 0.2, 0.2,
        0.2, 0.3, 0.2, 0.3
    ),
    byrow = TRUE,
    nrow = 3)
colnames(mutsig_by_kras) <- c(letters[1:3], "R")
rownames(mutsig_by_kras) <- rownames(mutsig_by_tumor)
mutsig_by_kras
# All rows should sum to one.
all(rowSums(mutsig_by_kras) == 1)


#> mutational signature x KRAS allele for the sample
mutsig_by_kras_for_T <- matrix(
    rep(mutsig_by_tumor, 4), nrow = nrow(mutsig_by_tumor)
) * mutsig_by_kras
mutsig_by_kras_for_T

# likelihood per allele by summing over the mutational signatures.
likelihood_per_allele <- colSums(mutsig_by_kras_for_T)
likelihood_per_allele

# Liklihood of getting KRAS allele "a" given will get a KRAS mutation.
res1 <- likelihood_per_allele[[1]] / sum(likelihood_per_allele[1:3])



## Same analysis, except removing the remainder and normalizing the alleles.
mutsig_by_kras_noR <- mutsig_by_kras[, -4]
mutsig_by_kras_noR <- t(apply(mutsig_by_kras_noR, 1, function(r) { r / sum(r) }))
mutsig_by_kras_noR

#> mutational signature x KRAS allele for the sample
mutsig_by_kras_for_T2 <- matrix(
    rep(mutsig_by_tumor, 3), nrow = nrow(mutsig_by_tumor)
) * mutsig_by_kras_noR
mutsig_by_kras_for_T2

# likelihood per allele by summing over the mutational signatures.
likelihood_per_allele2 <- colSums(mutsig_by_kras_for_T2)
likelihood_per_allele2

# Liklihood of getting KRAS allele "a" given will get a KRAS mutation.
res2 <- likelihood_per_allele2[[1]] / sum(likelihood_per_allele2[1:3])


# The second analysis does have the same result as the first.
message(glue("probability of KRAS allele 'a' without removing remainder: {round(res1, 4)}"))
message(glue("   probability of KRAS allele 'a' with removing remainder: {round(res2, 4)}"))

message("Conclusion: keep the remainder until the end.")
