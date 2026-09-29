# FSF biological validation

Status: **COMPLETE**

Generated from locked FSF/annotation authorities plus explicit manual author decisions and literature-evidence registries. No literature interpretation is hidden in the builder.

## Condition-level architecture evidence

| Condition | Low Stability | Transitional | Stable | Highly Stable | Evidence | Claim(s) |
|---|---:|---:|---:|---:|---:|---|
| DES | 38.43% | 0.00% | 0.00% | 61.57% | 1 | Ramazzottius varieornatus displays rapid-desiccation tolerance and comparatively limited transcriptional regulation during entry into anhydrobiosis, establishing a same-species transcriptomic context for desiccation tolerance. |
| GAM | 7.14% | 13.87% | 40.34% | 38.64% | 1 | Time-course RNA-seq of Ramazzottius varieornatus exposed to 500 Gy cobalt-60 gamma radiation demonstrated a measurable transcriptional response after irradiation and provides direct same-species gamma-response context. |
| HT | 0.00% | 0.00% | 0.00% | 100.00% | 1 | Heat exposure at 35 degrees C produced extensive transcriptomic remodeling in Ramazzottius varieornatus, involving more than one-third of the transcriptome and establishing direct same-species heat-stress context. |
| LT | 7.97% | 0.00% | 0.00% | 92.03% | 1 | Ramazzottius varieornatus survived extreme freezing while transcriptomic analysis detected no significant differential expression during cooling, providing direct same-species evidence for largely constitutive freeze-tolerance biology. |
| OSM | 0.00% | 0.00% | 0.00% | 100.00% | 1 | Ramazzottius varieornatus tolerated extreme increases in external osmolality, formed osmobiotic tuns, and showed a comparatively modest transcriptomic shift during osmotic stress, providing direct same-species osmobiosis context. |
| UV | 19.41% | 44.10% | 20.75% | 15.74% | 1 | Time-series transcriptomic analysis of ultraviolet-exposed Ramazzottius varieornatus identified stress-induced transcriptional responses and candidate factors associated with UV tolerance and anhydrobiosis cross-tolerance. |

## Validated gene shortlist

| Condition | FSF class | Feature | Gene | SSI | Decision | PMID/DOI | Manuscript use |
|---|---|---|---|---:|---|---|---|
| DES | Low Stability | RvY_00025-1 | CPEB1 | 0.5 | RETAIN_SUPPORTING | 15731006; 10.1242/jcs.01692 | Section 3.4 representative biological example |
| GAM | Stable Down | RvY_00648-1 | SOD1 | 0.875 | RETAIN_SUPPORTING | 31805397; 10.1016/j.freeradbiomed.2019.11.037 | Section 3.4 representative biological example |
| GAM | Highly Stable Constant | RvY_00174-1 | RPA2 | 1 | RETAIN_SUPPORTING | 11731442 | Section 3.4 representative biological example |
| HT | Highly Stable Constant | RvY_00096-1 | SMT3 | 1 | RETAIN_SUPPORTING | 21683690; 10.1016/j.bbrc.2011.06.025 | Section 3.4 representative biological example |
| LT | Highly Stable Up | RvY_02543-1 | TRPA1 | 1 | RETAIN_SUPPORTING | 19144922; 10.1073/pnas.0808487106 | Section 3.4 representative biological example |
| OSM | Highly Stable Up | RvY_00637-1 | SRSF3 | 1 | RETAIN_SUPPORTING | 31505169; 10.1016/j.bbamcr.2019.118557 | Section 3.4 representative biological example |
| UV | Transitional Down | RvY_01668-1 | PPM1D | 0.714285714285714 | RETAIN_SUPPORTING | 19015127; 10.1093/nar/gkn888 | Section 3.4 representative biological example |
| UV | Stable Up | RvY_00608-1 | SQSTM1 | 0.857142857142857 | RETAIN_SUPPORTING | 23340736; 10.1038/jid.2013.26 | Section 3.4 representative biological example |

## Guardrails

- Low Stability remains one non-directional FSF class.
- Annotation prioritizes literature review but is not independent validation.
- Retained genes require explicit author decisions and explicit PMID/DOI evidence.
- Condition relevance does not automatically prove the exact FSF direction.
