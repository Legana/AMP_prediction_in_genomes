Taxonomic composition of Swiss-Prot AMPs
================

UniProt AMPs were obtained from www.uniprot.org with the keyword
“antimicrobial” [KW-0929](https://www.uniprot.org/keywords/KW-0929).
On 14 January 2021 this returned 3,339 Reviewed proteins and 41,477
Unreviewed proteins.

All 3,339 reviewed entries were downloaded as a tab separated file which
included the sequence and taxonomic lineage (CLASS and ALL) as
additional columns. Note that the CLASS column was not used in the code
analysis and therefore is not strictly necessary. However, it helped me
as a visual tool to group the organisms together.

``` r
swissprot_amps <- read_tsv("data/uniprot-keyword Antimicrobial+[KW-0929] -filtered-reviewed yes.tab") %>%
    rename(Class = `Taxonomic lineage (CLASS)`) %>%
    rename(Taxonomic_lineage = `Taxonomic lineage (all)`)
```

Organisms were grouped according to their taxonomic class as much as
possible. However, to simplify, some cases were grouped together by
colloquial terms for cases where multiple Class arrangements were
present or where using a different taxonomic level was easier, e.g. the
order Testudines (which includes the turtles) and order Squamata (which
includes lizards and snakes) were combined into the “Reptiles” group.
Furthermore, Ascidiacea (which includes sea squirts), Amphioxus
(lancelets), Cnidaria (includes jellyfishes), Echinodermata (includes
sea stars), Mollusca (includes octopods), Merostomata (horseshoe crabs)
and Crustacea (includes crabs) were all grouped into “Marine
invertebrates” as their individually known AMP count was low.

``` r
swissprot_amps_grouped <- swissprot_amps %>%
mutate(Taxonomic_lineage = case_when(
   str_detect(Taxonomic_lineage, "Mammalia") ~ "Mammals",
   str_detect(Taxonomic_lineage, "Insecta") ~ "Insects",
   str_detect(Taxonomic_lineage, "Arachnida") ~ "Arachnids",
   str_detect(Taxonomic_lineage, "Bacteria") ~ "Bacteria",
   str_detect(Taxonomic_lineage, "Viruses") ~ "Bacteriophages",
   str_detect(Taxonomic_lineage, "Amphibia") ~ "Amphibians",
   str_detect(Taxonomic_lineage, "Viridiplantae") ~ "Plants",
   str_detect(Taxonomic_lineage, "Fungi") ~ "Fungi",
   str_detect(Taxonomic_lineage, "Annelida") ~ "Worms",
   str_detect(Taxonomic_lineage, "Trematoda") ~ "Worms",
   str_detect(Taxonomic_lineage, "Nematoda") ~ "Worms",
   str_detect(Taxonomic_lineage, "lamprey") ~ "Fishes",
   str_detect(Taxonomic_lineage, "Actinopteri") ~ "Fishes",
   str_detect(Taxonomic_lineage, "Testudines") ~ "Reptiles",
   str_detect(Taxonomic_lineage, "Squamata") ~ "Reptiles", 
   str_detect(Taxonomic_lineage, "Aves") ~ "Birds", 
   str_detect(Taxonomic_lineage, "Ascidiacea") ~ "Marine invertebrates", 
   str_detect(Taxonomic_lineage, "Amphioxus") ~ "Marine invertebrates", 
   str_detect(Taxonomic_lineage, "Cnidaria") ~ "Marine invertebrates", 
   str_detect(Taxonomic_lineage, "Echinodermata") ~ "Marine invertebrates", 
   str_detect(Taxonomic_lineage, "Mollusca") ~ "Marine invertebrates",
   str_detect(Taxonomic_lineage, "Merostomata") ~ "Marine invertebrates",
   str_detect(Taxonomic_lineage, "Crustacea") ~ "Marine invertebrates",
   str_detect(Taxonomic_lineage, "Amoebozoa") ~ "Amoebae",
   str_detect(Taxonomic_lineage, "Archaea") ~ "Archaea",
   str_detect(Taxonomic_lineage, "Chilopoda") ~ "Centipedes",
                                        TRUE ~ Taxonomic_lineage))
```

``` r
swissprot_amp_counts <- swissprot_amps_grouped %>%
  group_by(Taxonomic_lineage) %>% 
  summarise(amp_count = n()) %>% 
  mutate(percentage = round(amp_count / sum(amp_count) * 100, digits = 1))
```

Amphibians make up the majority of AMPs in the Swiss-Prot database
(\~24%), followed by mammals (\~19%), arthropods (\~20%) and plants
(16%). Marine invertebrates, which include several large and diverse
taxonomic groups, represented here by cnidarians, echinoderms and
crustaceans, comprise only 3.5% of the currently known AMPs. It is clear
that there is substantial taxonomic bias in the organism groups that
represent the AMPs in Swiss-Prot.

![](04_taxonomic_composition_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

**Figure 4.1:** Antimicrobial peptides (AMPs) present within groups of
organisms as obtained from SwissProt