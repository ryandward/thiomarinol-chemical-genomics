# Thiomarinol chemical genomics: figure code and data

Analysis code and data behind Figures 2 and 3 of:

> Johnson RA, Chan AN, **Ward RD**, McGlade CA, Hatfield BM, Peters JM, Li B.
> *Inhibition of Isoleucyl-tRNA Synthetase by the Hybrid Antibiotic Thiomarinol.*
> **J Am Chem Soc** 143(31):12003-12013 (2021).
> doi: [10.1021/jacs.1c02622](https://doi.org/10.1021/jacs.1c02622)
>
> Free full text: [PMC8479755](https://pmc.ncbi.nlm.nih.gov/articles/PMC8479755/) ·
> PMID [34342433](https://pubmed.ncbi.nlm.nih.gov/34342433/)

The experiment is a CRISPR interference (CRISPRi) knockdown screen of *Bacillus subtilis* essential
genes against a panel of chemical treatments. Each guide-condition pair gets an S-score: how far
that knockdown strain's fitness departs from what the rest of the library predicts. Negative means
the knockdown is sensitised to the compound, positive means it is protected.

- **Figure 2.** Median S-score for every guide, ordered by position in the genome, with the top and
  bottom 2.5% labelled. `ileS` is the most sensitised knockdown under both mupirocin (2A) and
  PAC-holo (2B), the paper's central genetic evidence that both compounds hit isoleucyl-tRNA
  synthetase.
- **Figure 3.** Pearson correlation between the 36 treatments' chemical-gene profiles, clustered by
  average linkage. PAC-holo correlates with mupirocin at r ≈ 0.37, while holomycin, thiolutin and
  gliotoxin form a tight separate cluster at r ≈ 0.75.

## Reproducing the figures

Requires R (tested on 4.6.1) with `data.table`, `ggplot2`, `ggrepel`, `pheatmap`, `viridis` and
`svglite`. Run from the repository root:

```sh
Rscript R/01_figure2_outlier_plots.R        # figures/{Mupirocin,PAC-holo,...}.svg
Rscript R/02_figure3_correlation_heatmap.R  # figures/MedianChemCorByCondition.{svg,png,tsv}
```

Rendering is not pixel-identical across R and package versions: `ggrepel` positions outlier labels
differently, so labels and their leader lines shift slightly.

## The data

`data/raw/Bsu_CRISPRi_S-scores_w_Bo_chemicals.tsv` is the measurement table everything derives from.
One row is a library-wide profile of one chemical at one concentration, scored across 300 guides.
There are 104 such profiles covering 38 conditions: 36 chemical treatments, plus UV and growth in
minimal S7-50 medium, which the published analysis excludes.

The same table is published as Data Set 1 of the paper's Supplemental Datasets, under slightly
different row labels. Prefer that version if you are starting fresh, since it records the numeric
concentration for each row where the labels here keep only the dose letter.

The comparison panel is the published reference compendium, as the paper states: the four compounds
of this study were screened alongside "31 other growth perturbing chemicals that have published
chemical-gene interactions". Concretely, 93 of the 104 profiles, covering 34 of the 38 conditions,
come from the *B. subtilis* chemical-genomics screen of:

> Peters JM, Colavin A, Shi H, Czarny TL, Larson MH, Wong S, Hawkins JS, Lu CHS, Koo B-M,
> Marta E, Shiver AL, Whitehead EH, Weissman JS, Brown ED, Qi LS, Huang KC, Gross CA.
> *A Comprehensive, CRISPR-based Functional Analysis of Essential Genes in Bacteria.*
> **Cell** 165(6):1493-1506 (2016). doi: [10.1016/j.cell.2016.05.003](https://doi.org/10.1016/j.cell.2016.05.003)

published there as Table S2. Those 93 rows match Table S2 to better than 1e-6 across the 280 guides
the two tables name identically. The remaining 11 profiles are the four compounds this study
contributed: holomycin, thiolutin and gliotoxin at three concentrations each, and PAC-holo at two.
That is what the "w_Bo_chemicals" in the filename refers to.

Both figures are new to the 2021 paper. Mupirocin's profile is drawn from the reference set, which
is what makes it a benchmark rather than a fresh measurement to be compared against.

Twenty guides are named differently between the two tables, for example `acpS50_1` and `acpS50_2`
in Table S2 against `acpS` and `acpS48` here, and `ftsZ226` and `ftsZ1134` against `ftsZ` and
`ftsZ1132`. Cross-reference by value rather than by name.

Raw S-scores run from -16.25 to 14.49, with four missing values. After taking medians over
concentration the range narrows to -9.05 to 9.52.

Most genes carry a single guide, but seven loci carry several at different knockdown strengths,
labelled with a numeric suffix: `ylaN`, `ylaN25`, `ylaN85`, `ylaN187` and `ylaN190` all target
BSU14840. These are distinct strains with distinct phenotypes, so a table with one row per guide is
longer than the number of genes it covers.

### Concentration series

The bracketed letter in a condition label indexes concentration. Data Set 1 writes the dose out in
full, `mupirocin [A] 0.01` through `mupirocin [G] 0.07` in ug/mL.

The median behind both figures is taken over a chemical's concentration series. S-scores arrive
already averaged over biological replicates, done upstream when they were computed (Peters et al.
2016, above); neither those replicates nor the colony-size measurements behind them are in this
repository.

Dose coverage is uneven and worth checking before relying on any single condition. Trimethoprim has
10 concentrations and mupirocin 7, but PAC-holo has 2, and six conditions were assayed at a single
concentration: actinomycin D, cefaclor, chlorhexidine, nalidixic acid, and the two non-chemical
conditions the analysis drops. For those, the median is one measurement. A condition's column can
therefore rest on very different amounts of data, and nothing in the correlation matrix reflects
that.

## Layout

```
data/raw/
  Bsu_CRISPRi_S-scores_w_Bo_chemicals.tsv   104 chemical-concentration profiles x 300 guides
  conditions.tsv                            38 condition names
  BSUtoGene.tsv                             332 locus tag to guide labels
  annotations/peptidoglycan.tsv             39 genes, the blue points in Fig 2B
  annotations/cell_shape.tsv                43 genes, the yellow points in Fig 2B
data/processed/
  median_s_scores_by_condition.tsv          297 guides x 36 conditions
  median_s_scores_by_guide.tsv              306 rows, ordered by locus tag, as plotted in Fig 2
  correlation_by_condition.tsv              36 x 36 Pearson correlation (Fig 3)
  heatmap_correlation_second_order.tsv      36 x 36 correlation of the correlation profiles
R/
  00_common.R                               paths, and the raw to median-matrix step
  01_figure2_outlier_plots.R
  02_figure3_correlation_heatmap.R
figures/                                    output
```

`data/processed/` is derived from `data/raw/` by the scripts here. It is committed so the matrices
can be used directly without running R.

The gene lists in `data/raw/annotations/` set the colouring in Figure 2B. They come from the GO
analysis published as Data Set 3, which was run in STRING rather than in this repository.

The Figure 2 panels share a y-axis of -7.1 to 4.9, wide enough for all five panels. Pointing the
script at another condition may exceed it, since the full matrix reaches -9.05; the script warns if
that happens.

## Licence

Code in `R/` is MIT licensed (`LICENSE`). The data accompanies the published papers; if you use it,
cite Johnson et al. 2021, and cite Peters et al. 2016 as well for any of the 34 conditions that
originate there.
