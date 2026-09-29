# Current FSF biological validation

This workflow has three explicit layers:

1. Automatic deterministic candidate generation from the locked FSF and annotation authorities.
2. Manual author decisions and literature evidence stored in TSV registries.
3. Generated validated shortlist, architecture table, summary, and manifest.

`candidate_shortlist_top5.tsv` is a literature-review queue, not biological validation.

Low Stability is collapsed to one non-directional class before ranking.

The builder never searches the web or invents literature evidence.

PMID/DOI records and evidence statements are explicit version-controlled inputs.

Modes:

- `--initialize`
- `--build`
- `--verify`
