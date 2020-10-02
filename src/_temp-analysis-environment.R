library(mustashe)
library(jhcutils)
library(glue)
library(ggridges)
library(patchwork)
library(ggtext)
library(tidyverse)

for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
  if (str_detect(f, "global|enrich")) {
    next
  }
  source(f)
}

options(dplyr.summarise.inform = FALSE)


theme_set(
  theme_bw(base_size = 11, base_family = "Arial") %+replace%
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.ticks = element_blank()
    )
)

load("~/Downloads/kras_trinucleotide_contexts.RData")
load("~/Downloads/mutational_signature_spectra.RData")


#### ---- Diagram for figure 1 ---- ####

kras_trinucleotide_contexts
mutational_signature_spectra


example_ts <- tibble(
  signature = c("1", "5", "18"),
  composition = c(0.6, 0.1, 0.3)
)
example_ts


example_ts %>%
  ggplot(aes(x = "a", y = composition, fill = signature)) +
  geom_col(position = "stack") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text.x = element_blank()
  )


mutational_signature_spectra %>%
  mutate(signature = fct_inorder(signature)) %>%
  filter(signature %in% example_ts$signature) %>%
  left_join(
    kras_trinucleotide_contexts %>%
      filter(kras_allele %in% c("G12C", "G12D")) %>%
      select(kras_allele, tricontext),
    by = c("tricontext")
  ) %>%
  ggplot(aes(x = tricontext, y = composition)) +
  facet_wrap(~ signature, ncol = 1) +
  geom_col(aes(fill = kras_allele)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = short_allele_pal,
    drop = TRUE,
    na.value = "grey75"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.text = element_blank(),
    panel.grid.major.x = element_blank()
  )
