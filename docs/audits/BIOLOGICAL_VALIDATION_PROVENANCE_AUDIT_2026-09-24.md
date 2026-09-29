# Biological-validation provenance audit

- Audit date: 2026-09-24
- Git branch: `fsf-manuscript-revision`
- HEAD commit: `7facc863c6e3f0b3375ccd37e7e8bbccf647b7f0`
- Audit mode: read-only. The provenance audit did not edit or regenerate results, change algorithms or rankings, modify Figure 5 or literature evidence, commit, or push.

## Candidate-shortlist provenance

`candidate_shortlist_top5.tsv` is constructed by `scripts/current/build_current_biological_validation.py` from these locked authorities, joined by `feature_id`:

- `results/current_fsf_v1/current_fsf_feature_metrics.tsv`
- `/media/saji/5E06441D0643F5152/FSF_R_PACKAGE_ARTIFACTS/FSF_v1_annotation_master.tsv`

Candidates are grouped by `(condition, validation_class)`. For Low Stability, `validation_class` is the single nondirectional `Low Stability` class. Otherwise, `validation_class` is `signal_class`, retaining Up, Constant, or Down directionality. The 34 observed validation groups comprise DES=4, GAM=10, HT=3, LT=4, OSM=3, and UV=10.

A feature is eligible when it has a named identity (a present eggNOG preferred name, or a present non-domain-only eggNOG description); its name/description does not match the builder's vague-annotation exclusions; `annotation_status` is present; `annotation_source_count >= 1`; and at least one of eggNOG GO, eggNOG KEGG KO, eggNOG KEGG pathway, InterPro, or Pfam annotation is present.

Within each group, candidates are sorted by this exact hierarchy:

1. Strong candidate before merely eligible candidate. Strong requires named identity, usable annotation status, at least two annotation sources, and at least two supported functional-annotation categories.
2. `annotation_source_count`, descending.
3. `functional_support_count`, descending.
4. Preferred name present, YES first.
5. KEGG pathway annotation present, YES first.
6. KEGG KO annotation present, YES first.
7. GO annotation present, YES first.
8. InterPro annotation present, YES first.
9. Pfam annotation present, YES first.
10. SSI, descending.
11. `feature_id`, ascending lexical order as the deterministic final tie-breaker.

The first five candidates per group receive ranks 1 through 5. Literature evidence and manual author decisions do not participate in top-five generation. Significant GO/KEGG enrichment, term membership, enrichment ratios, p-values, and adjusted p-values do not participate directly in ranking. Generic presence/absence of GO and KEGG annotation does participate as described above.

## Relationship to Figure 5

`FUNCTIONAL_ENRICHMENT_LINK=PARTIAL`: the shortlist prioritizes generic GO/KEGG annotation coverage but does not use the significant enrichment evidence displayed in Figure 5.

`GENE_TERM_TRACEABILITY=YES`: the significant Figure 5 GO and KEGG source tables contain `contributing_feature_ids` for each significant term. Together with `condition`, `stability_region`, `gene_set_id`, term ID, and the retained class fields, these values permit deterministic gene-to-significant-term tracing. The Figure 5 source tables are byte-identical to their corresponding significant-enrichment tables.

## SHA-256 record

### Candidate-shortlist authorities and output

```text
2a0a064573fd02647284e0fbfbb3b9d862b42b6a56ab5ed98ec7a816fd227028  scripts/current/build_current_biological_validation.py
e7d21744681122f041ad38623bb4cb4807ab3d254132fffd6269ce0ef0b1ab89  results/current_fsf_v1/current_fsf_feature_metrics.tsv
25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584  /media/saji/5E06441D0643F5152/FSF_R_PACKAGE_ARTIFACTS/FSF_v1_annotation_master.tsv
763913ae0780934c5a1c75c1fef9813d8f180b424c5eabad513cb88160c9b37d  results/current_fsf_v1/manuscript/biological_validation/candidate_shortlist_top5.tsv
```

### Figure 5 enrichment generation and authorities

```text
6a783f6d85dae71473d9dc2c5db053723deb4808084fc1f2cdf81956be3a651d  scripts/lib/current_fsf_annotation_authority.R
48c548930a406159fe20f5608ced47342c2e2d6c20ed2b5eb7502ccb34ab0063  scripts/current/build_current_enrichment_gene_sets.R
c8afacee3de3cf057ff372748b1425d79a00867d63cbd1fa15cfc855a116e09f  scripts/current/run_current_go_enrichment.R
ac156d37f94cd2b063d5210e00428dcc0068b8b8947ff032b0129f2e3e6d2c9f  scripts/current/run_current_kegg_enrichment.R
9820de30384aa26949885f097a777396a510a0f5eb061088d78459d6c849cd39  scripts/current/build_current_manuscript_enrichment_outputs.R
b05a21fe6c145eb1c7a794dc60797139f6a215d1251577981aec4f728144aabb  scripts/current/build_current_manuscript_tables.R
0ded08c5f67abdcf6dbbc2d77bcd47530319ecdeeb7b6ad13cd402a65d5dcdc2  results/current_fsf_v1/manuscript/enrichment/go/go_enrichment_significant.tsv
780415dbe1a73d9246b16b58b58892b24f20d89eb1d3cf9f6458479d2c2f856b  results/current_fsf_v1/manuscript/enrichment/kegg/kegg_enrichment_significant.tsv
0ded08c5f67abdcf6dbbc2d77bcd47530319ecdeeb7b6ad13cd402a65d5dcdc2  results/current_fsf_v1/manuscript/source_data/figure5_go_source.tsv
780415dbe1a73d9246b16b58b58892b24f20d89eb1d3cf9f6458479d2c2f856b  results/current_fsf_v1/manuscript/source_data/figure5_kegg_source.tsv
7a2f35a602cf4703ee341726f9dc86e207cba5e063a964a58614c73b80bf98f7  scripts/current/figures/build_current_manuscript_figures.R
e106bba69f1a7069463c29b390a23876bc15a4e15cf57c1c99cc38e5e6e807cb  results/current_fsf_v1/manuscript/figures/biological_theme_prototypes/Figure5C_DES_biological_themes_barplot_AUDIT.png
436180835ad7cf0b15def1aad2ff2b8184b577d4950c05e0ebd6d9c96e9ad671  results/current_fsf_v1/manuscript/figures/main/Figure5.png
e7738da021371a342c72bd5d67955d2bb0b901abf20952ecb6afb6d65ddf613e  results/current_fsf_v1/manuscript/figures/main/Figure5.pdf
```

## Final audit conclusions

```text
CURRENT_TOP5_LOGIC=Top five annotation-rich eligible genes per condition × collapsed/directional FSF validation class, ranked by annotation strength/coverage, then descending SSI, then ascending feature_id
FUNCTIONAL_ENRICHMENT_LINK=PARTIAL
GENE_TERM_TRACEABILITY=YES
MANUAL_JUDGMENT_IN_TOP5=NO
LITERATURE_IN_TOP5=NO
SAFE_TO_DESIGN_FINAL_SHORTLIST=YES
```
