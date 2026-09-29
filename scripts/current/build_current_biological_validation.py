#!/usr/bin/env python3
from __future__ import annotations

import csv
import filecmp
import hashlib
import os
import re
import shutil
import sys
import tempfile

from collections import Counter, defaultdict
from pathlib import Path


FSF_SHA = {
    "e7d21744681122f041ad38623bb4cb4807ab3d254132fffd6269ce0ef0b1ab89",
    "d4384daa37a2532cbb601bab87556c2a5e8ad3d87a73690f6ddc4f04cce95a57",
}

ANN_SHA = (
    "25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584"
)

N_FSF = 91122
N_FEATURES = 15187

# After Low Stability is correctly collapsed to one non-directional class:
#
# DES = 4
# GAM = 10
# HT  = 3
# LT  = 4
# OSM = 3
# UV  = 10
#
# Total = 34 biologically interpretable validation groups.
N_GROUPS = 34

# Five candidates per validation group.
N_TOP5 = 170


CONDS = (
    "DES",
    "GAM",
    "HT",
    "LT",
    "OSM",
    "UV",
)


REGIONS = (
    "Low Stability",
    "Transitional",
    "Stable",
    "Highly Stable",
)


CLASS_ORDER = (
    "Low Stability",
    "Transitional Up",
    "Transitional Constant",
    "Transitional Down",
    "Stable Up",
    "Stable Constant",
    "Stable Down",
    "Highly Stable Up",
    "Highly Stable Constant",
    "Highly Stable Down",
)


VAGUE = re.compile(
    r"hypothetical|"
    r"uncharacteri[sz]ed|"
    r"unknown function|"
    r"function unknown|"
    r"domain of unknown function|"
    r"\bDUF[0-9]+\b|"
    r"anonymous|"
    r"unnamed|"
    r"weak fragment|"
    r"protein fragment",
    re.I,
)


DOMAIN = re.compile(
    r"^domain$|"
    r"^.* domain$|"
    r"^.* repeat$|"
    r"^.* repeats$|"
    r"^.*-like domain$|"
    r"^protein of unknown function",
    re.I,
)


MANUAL = (
    "gene_author_decisions.tsv",
    "gene_literature_evidence.tsv",
    "condition_literature_evidence.tsv",
)


GENERATED = (
    "condition_architecture.tsv",
    "candidate_shortlist_top5.tsv",
    "validated_gene_shortlist.tsv",
    "validated_condition_architecture.tsv",
    "biological_validation_summary.md",
    "README.md",
    "biological_validation_manifest.tsv",
)


def csv_limit():
    value = sys.maxsize

    while True:
        try:
            csv.field_size_limit(value)
            return
        except OverflowError:
            value //= 10


def txt(value):
    if value is None:
        return ""

    return str(value).strip()


def present(value):
    value = txt(value)

    return (
        bool(value)
        and value.lower()
        not in {
            "-",
            "na",
            "n/a",
            "none",
            "null",
            "missing",
        }
    )


def num(value, default=-1e99):
    try:
        return float(txt(value))
    except Exception:
        return default


def sha(path):
    path = Path(path)

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def read_tsv(path):
    path = Path(path)

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise AssertionError(
                f"Missing header: {path}"
            )

        return (
            list(reader.fieldnames),
            list(reader),
        )


def write_tsv(path, fields, rows):
    path = Path(path)

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: txt(
                        row.get(
                            field,
                            "",
                        )
                    )
                    for field in fields
                }
            )


def require(header, fields, label):
    missing = [
        field
        for field in fields
        if field not in header
    ]

    if missing:
        raise AssertionError(
            (
                label,
                "MISSING_FIELDS",
                missing,
            )
        )


def paths():
    repo = Path(
        os.environ.get(
            "FSF_REPO_ROOT",
            "",
        )
    ).expanduser().resolve()

    annotation = Path(
        os.environ.get(
            "FSF_ANNOTATION_MASTER",
            "",
        )
    ).expanduser().resolve()

    fsf = (
        repo
        / "results"
        / "current_fsf_v1"
        / "current_fsf_feature_metrics.tsv"
    )

    root = (
        repo
        / "results"
        / "current_fsf_v1"
        / "manuscript"
        / "biological_validation"
    )

    if (
        not repo.is_dir()
        or not fsf.is_file()
        or not annotation.is_file()
    ):
        raise SystemExit(
            "Set valid FSF_REPO_ROOT "
            "and FSF_ANNOTATION_MASTER."
        )

    return (
        repo,
        fsf,
        annotation,
        root,
    )


def authorities(fsf, annotation):

    if sha(fsf) not in FSF_SHA:
        raise AssertionError(
            "FSF_SHA256_MISMATCH"
        )

    if sha(annotation) != ANN_SHA:
        raise AssertionError(
            "ANNOTATION_SHA256_MISMATCH"
        )

    fsf_header, fsf_rows = read_tsv(
        fsf
    )

    annotation_header, annotation_rows = read_tsv(
        annotation
    )

    require(
        fsf_header,
        (
            "feature_id",
            "condition",
            "dominant_state",
            "ssi",
            "stability_region",
            "signal_class",
        ),
        "FSF",
    )

    require(
        annotation_header,
        (
            "feature_id",
            "locus_tag",
            "protein_id",
            "eggnog_preferred_name",
            "eggnog_description",
            "eggnog_go",
            "eggnog_kegg_ko",
            "eggnog_kegg_pathway",
            "interpro_descriptions",
            "pfam_descriptions",
            "annotation_source_count",
            "annotation_status",
        ),
        "ANNOTATION",
    )

    if (
        len(fsf_rows) != N_FSF
        or len(annotation_rows)
        != N_FEATURES
    ):
        raise AssertionError(
            (
                "ROW_COUNT",
                len(fsf_rows),
                len(annotation_rows),
            )
        )

    annotation_map = {}

    for row in annotation_rows:

        feature_id = txt(
            row["feature_id"]
        )

        if (
            not feature_id
            or feature_id
            in annotation_map
        ):
            raise AssertionError(
                (
                    "ANNOTATION_KEY",
                    feature_id,
                )
            )

        preferred_name = txt(
            row[
                "eggnog_preferred_name"
            ]
        )

        description = txt(
            row[
                "eggnog_description"
            ]
        )

        has_preferred_name = present(
            preferred_name
        )

        has_description = present(
            description
        )

        named_identity = (
            (
                has_preferred_name
                or (
                    has_description
                    and not DOMAIN.search(
                        description
                    )
                )
            )
            and not VAGUE.search(
                preferred_name
                + " | "
                + description
            )
        )

        annotation_flags = {
            key: present(
                row[key]
            )
            for key in (
                "eggnog_go",
                "eggnog_kegg_ko",
                "eggnog_kegg_pathway",
                "interpro_descriptions",
                "pfam_descriptions",
            )
        }

        support_count = sum(
            annotation_flags.values()
        )

        source_text = txt(
            row[
                "annotation_source_count"
            ]
        )

        source_count = int(
            float(
                source_text
                or 0
            )
        )

        usable_status = present(
            row[
                "annotation_status"
            ]
        )

        eligible = (
            named_identity
            and usable_status
            and source_count >= 1
            and support_count >= 1
        )

        strong = (
            named_identity
            and usable_status
            and source_count >= 2
            and support_count >= 2
        )

        annotation_map[
            feature_id
        ] = {
            "locus_tag":
                txt(
                    row[
                        "locus_tag"
                    ]
                ),
            "protein_id":
                txt(
                    row[
                        "protein_id"
                    ]
                ),
            "gene":
                preferred_name,
            "description":
                description,
            "kegg_ko":
                txt(
                    row[
                        "eggnog_kegg_ko"
                    ]
                ),
            "kegg_pathway":
                txt(
                    row[
                        "eggnog_kegg_pathway"
                    ]
                ),
            "interpro":
                txt(
                    row[
                        "interpro_descriptions"
                    ]
                ),
            "pfam":
                txt(
                    row[
                        "pfam_descriptions"
                    ]
                ),
            "source_count":
                source_count,
            "annotation_status":
                txt(
                    row[
                        "annotation_status"
                    ]
                ),
            "support_count":
                support_count,
            "has_pref":
                has_preferred_name,
            "has_go":
                annotation_flags[
                    "eggnog_go"
                ],
            "has_ko":
                annotation_flags[
                    "eggnog_kegg_ko"
                ],
            "has_path":
                annotation_flags[
                    "eggnog_kegg_pathway"
                ],
            "has_interpro":
                annotation_flags[
                    "interpro_descriptions"
                ],
            "has_pfam":
                annotation_flags[
                    "pfam_descriptions"
                ],
            "eligible":
                bool(
                    eligible
                ),
            "strong":
                bool(
                    strong
                ),
        }

    fsf_features = {
        txt(
            row[
                "feature_id"
            ]
        )
        for row in fsf_rows
    }

    if (
        len(fsf_features)
        != N_FEATURES
        or fsf_features
        != set(
            annotation_map
        )
    ):
        raise AssertionError(
            "FEATURE_UNIVERSE_MISMATCH"
        )

    if {
        txt(
            row[
                "condition"
            ]
        )
        for row in fsf_rows
    } != set(
        CONDS
    ):
        raise AssertionError(
            "CONDITION_SET_MISMATCH"
        )

    fsf_keys = [
        (
            txt(
                row[
                    "condition"
                ]
            ),
            txt(
                row[
                    "feature_id"
                ]
            ),
        )
        for row in fsf_rows
    ]

    if (
        len(fsf_keys)
        != len(
            set(
                fsf_keys
            )
        )
    ):
        raise AssertionError(
            "DUPLICATE_CONDITION_FEATURE"
        )

    return (
        fsf_rows,
        annotation_map,
    )


def architecture(fsf_rows):

    counts = Counter(
        (
            txt(
                row[
                    "condition"
                ]
            ),
            txt(
                row[
                    "stability_region"
                ]
            ),
        )
        for row in fsf_rows
    )

    totals = Counter(
        txt(
            row[
                "condition"
            ]
        )
        for row in fsf_rows
    )

    output = []

    for condition in CONDS:

        if (
            totals[
                condition
            ]
            != N_FEATURES
        ):
            raise AssertionError(
                (
                    "CONDITION_FEATURE_COUNT",
                    condition,
                    totals[
                        condition
                    ],
                )
            )

        region_counts = {
            region:
                counts[
                    (
                        condition,
                        region,
                    )
                ]
            for region in REGIONS
        }

        total = totals[
            condition
        ]

        output.append(
            {
                "condition":
                    condition,
                "total_features":
                    total,
                "low_stability_count":
                    region_counts[
                        "Low Stability"
                    ],
                "transitional_count":
                    region_counts[
                        "Transitional"
                    ],
                "stable_count":
                    region_counts[
                        "Stable"
                    ],
                "highly_stable_count":
                    region_counts[
                        "Highly Stable"
                    ],
                "low_stability_proportion":
                    f"{region_counts['Low Stability']/total:.8f}",
                "transitional_proportion":
                    f"{region_counts['Transitional']/total:.8f}",
                "stable_proportion":
                    f"{region_counts['Stable']/total:.8f}",
                "highly_stable_proportion":
                    f"{region_counts['Highly Stable']/total:.8f}",
            }
        )

    return output


ARCH_FIELDS = (
    "condition",
    "total_features",
    "low_stability_count",
    "transitional_count",
    "stable_count",
    "highly_stable_count",
    "low_stability_proportion",
    "transitional_proportion",
    "stable_proportion",
    "highly_stable_proportion",
)


CAND_FIELDS = (
    "condition",
    "validation_class",
    "stability_region",
    "dominant_state",
    "signal_class",
    "candidate_rank",
    "feature_id",
    "locus_tag",
    "protein_id",
    "ssi",
    "p_up",
    "p_down",
    "p_const",
    "stability_deviation",
    "eggnog_preferred_name",
    "eggnog_description",
    "eggnog_kegg_ko",
    "eggnog_kegg_pathway",
    "interpro_descriptions",
    "pfam_descriptions",
    "annotation_source_count",
    "annotation_status",
    "functional_support_count",
    "has_go",
    "has_kegg_ko",
    "has_kegg_pathway",
    "has_interpro",
    "has_pfam",
    "candidate_strength",
    "directionality_interpretation",
)


DEC_FIELDS = (
    "condition",
    "validation_class",
    "feature_id",
    "candidate_rank",
    "eggnog_preferred_name",
    "author_decision",
    "manuscript_use",
    "author_notes",
)


GENE_EVID_FIELDS = (
    "evidence_id",
    "condition",
    "validation_class",
    "feature_id",
    "pmid",
    "doi",
    "citation",
    "evidence_role",
    "supports_condition_relevance",
    "supports_directionality",
    "evidence_notes",
)


COND_EVID_FIELDS = (
    "evidence_id",
    "condition",
    "pmid",
    "doi",
    "citation",
    "evidence_claim",
    "evidence_scope",
    "author_decision",
    "author_notes",
)


def candidates(
    fsf_rows,
    annotation_map,
):

    groups = defaultdict(
        list
    )

    for row in fsf_rows:

        feature_id = txt(
            row[
                "feature_id"
            ]
        )

        annotation = annotation_map[
            feature_id
        ]

        if not annotation[
            "eligible"
        ]:
            continue

        region = txt(
            row[
                "stability_region"
            ]
        )

        # ----------------------------------------------------
        # SCIENTIFIC RULE:
        # Low Stability is one non-directional class.
        # Internal dominant_state is retained only for
        # provenance and never creates directional
        # Low-Stability biological classes.
        # ----------------------------------------------------

        if (
            region
            == "Low Stability"
        ):
            validation_class = (
                "Low Stability"
            )
        else:
            validation_class = txt(
                row[
                    "signal_class"
                ]
            )

        record = {
            "condition":
                txt(
                    row[
                        "condition"
                    ]
                ),
            "validation_class":
                validation_class,
            "stability_region":
                region,
            "dominant_state":
                txt(
                    row[
                        "dominant_state"
                    ]
                ),
            "signal_class":
                txt(
                    row[
                        "signal_class"
                    ]
                ),
            "feature_id":
                feature_id,
            "ssi":
                txt(
                    row[
                        "ssi"
                    ]
                ),
            "p_up":
                txt(
                    row.get(
                        "p_up",
                        "",
                    )
                ),
            "p_down":
                txt(
                    row.get(
                        "p_down",
                        "",
                    )
                ),
            "p_const":
                txt(
                    row.get(
                        "p_const",
                        "",
                    )
                ),
            "stability_deviation":
                txt(
                    row.get(
                        "stability_deviation",
                        "",
                    )
                ),
            "locus_tag":
                annotation[
                    "locus_tag"
                ],
            "protein_id":
                annotation[
                    "protein_id"
                ],
            "eggnog_preferred_name":
                annotation[
                    "gene"
                ],
            "eggnog_description":
                annotation[
                    "description"
                ],
            "eggnog_kegg_ko":
                annotation[
                    "kegg_ko"
                ],
            "eggnog_kegg_pathway":
                annotation[
                    "kegg_pathway"
                ],
            "interpro_descriptions":
                annotation[
                    "interpro"
                ],
            "pfam_descriptions":
                annotation[
                    "pfam"
                ],
            "annotation_source_count":
                annotation[
                    "source_count"
                ],
            "annotation_status":
                annotation[
                    "annotation_status"
                ],
            "functional_support_count":
                annotation[
                    "support_count"
                ],
            "has_go":
                (
                    "YES"
                    if annotation[
                        "has_go"
                    ]
                    else "NO"
                ),
            "has_kegg_ko":
                (
                    "YES"
                    if annotation[
                        "has_ko"
                    ]
                    else "NO"
                ),
            "has_kegg_pathway":
                (
                    "YES"
                    if annotation[
                        "has_path"
                    ]
                    else "NO"
                ),
            "has_interpro":
                (
                    "YES"
                    if annotation[
                        "has_interpro"
                    ]
                    else "NO"
                ),
            "has_pfam":
                (
                    "YES"
                    if annotation[
                        "has_pfam"
                    ]
                    else "NO"
                ),
            "candidate_strength":
                (
                    "STRONG"
                    if annotation[
                        "strong"
                    ]
                    else "ELIGIBLE"
                ),
            "directionality_interpretation":
                (
                    "NON_DIRECTIONAL"
                    if region
                    == "Low Stability"
                    else
                    "DIRECTIONAL_CLASS_ALLOWED"
                ),
        }

        record[
            "_rank"
        ] = (
            -int(
                annotation[
                    "strong"
                ]
            ),
            -annotation[
                "source_count"
            ],
            -annotation[
                "support_count"
            ],
            -int(
                annotation[
                    "has_pref"
                ]
            ),
            -int(
                annotation[
                    "has_path"
                ]
            ),
            -int(
                annotation[
                    "has_ko"
                ]
            ),
            -int(
                annotation[
                    "has_go"
                ]
            ),
            -int(
                annotation[
                    "has_interpro"
                ]
            ),
            -int(
                annotation[
                    "has_pfam"
                ]
            ),
            -num(
                row[
                    "ssi"
                ]
            ),
            feature_id,
        )

        groups[
            (
                record[
                    "condition"
                ],
                validation_class,
            )
        ].append(
            record
        )

    if (
        len(groups)
        != N_GROUPS
    ):
        raise AssertionError(
            (
                "VALIDATION_GROUP_COUNT",
                len(groups),
            )
        )

    condition_order = {
        value: index
        for index, value
        in enumerate(
            CONDS
        )
    }

    class_order = {
        value: index
        for index, value
        in enumerate(
            CLASS_ORDER
        )
    }

    output = []

    for key in sorted(
        groups,
        key=lambda item: (
            condition_order[
                item[0]
            ],
            class_order[
                item[1]
            ],
        ),
    ):

        rows = sorted(
            groups[
                key
            ],
            key=lambda row:
                row[
                    "_rank"
                ],
        )

        if len(rows) < 5:
            raise AssertionError(
                (
                    "TOO_FEW_CANDIDATES",
                    key,
                    len(rows),
                )
            )

        for rank, row in enumerate(
            rows[:5],
            1,
        ):
            row = dict(
                row
            )

            row[
                "candidate_rank"
            ] = rank

            output.append(
                row
            )

    if (
        len(output)
        != N_TOP5
    ):
        raise AssertionError(
            (
                "TOP5_ROW_COUNT",
                len(output),
            )
        )

    return output


def init_manual(
    root,
    candidate_rows,
):

    decisions = (
        root
        / "gene_author_decisions.tsv"
    )

    if not decisions.exists():

        write_tsv(
            decisions,
            DEC_FIELDS,
            (
                {
                    "condition":
                        row[
                            "condition"
                        ],
                    "validation_class":
                        row[
                            "validation_class"
                        ],
                    "feature_id":
                        row[
                            "feature_id"
                        ],
                    "candidate_rank":
                        row[
                            "candidate_rank"
                        ],
                    "eggnog_preferred_name":
                        row[
                            "eggnog_preferred_name"
                        ],
                    "author_decision":
                        "PENDING",
                    "manuscript_use":
                        "",
                    "author_notes":
                        "",
                }
                for row
                in candidate_rows
            ),
        )

    gene_evidence = (
        root
        / "gene_literature_evidence.tsv"
    )

    if not gene_evidence.exists():

        write_tsv(
            gene_evidence,
            GENE_EVID_FIELDS,
            [],
        )

    condition_evidence = (
        root
        / "condition_literature_evidence.tsv"
    )

    if not condition_evidence.exists():

        write_tsv(
            condition_evidence,
            COND_EVID_FIELDS,
            [],
        )


def manual(
    root,
    candidate_rows,
):

    candidate_map = {
        (
            row[
                "condition"
            ],
            row[
                "validation_class"
            ],
            row[
                "feature_id"
            ],
        ):
            row
        for row in candidate_rows
    }

    decision_header, decisions = read_tsv(
        root
        / "gene_author_decisions.tsv"
    )

    require(
        decision_header,
        DEC_FIELDS,
        "GENE_DECISIONS",
    )

    decision_keys = []

    for row in decisions:

        key = (
            txt(
                row[
                    "condition"
                ]
            ),
            txt(
                row[
                    "validation_class"
                ]
            ),
            txt(
                row[
                    "feature_id"
                ]
            ),
        )

        decision_keys.append(
            key
        )

        if (
            key
            not in candidate_map
        ):
            raise AssertionError(
                (
                    "STALE_DECISION_KEY",
                    key,
                )
            )

        current = candidate_map[
            key
        ]

        if (
            txt(
                row[
                    "candidate_rank"
                ]
            )
            != txt(
                current[
                    "candidate_rank"
                ]
            )
            or
            txt(
                row[
                    "eggnog_preferred_name"
                ]
            )
            != txt(
                current[
                    "eggnog_preferred_name"
                ]
            )
        ):
            raise AssertionError(
                (
                    "STALE_DECISION_ROW",
                    key,
                )
            )

        decision = (
            txt(
                row[
                    "author_decision"
                ]
            )
            or "PENDING"
        )

        if decision not in {
            "PENDING",
            "RETAIN_HEADLINE",
            "RETAIN_SUPPORTING",
            "REJECT",
        }:
            raise AssertionError(
                (
                    "BAD_GENE_DECISION",
                    key,
                )
            )

    if (
        len(decision_keys)
        != len(
            set(
                decision_keys
            )
        )
        or
        set(
            decision_keys
        )
        != set(
            candidate_map
        )
    ):
        raise AssertionError(
            "GENE_DECISION_COVERAGE_MISMATCH"
        )

    gene_header, gene_evidence = read_tsv(
        root
        / "gene_literature_evidence.tsv"
    )

    require(
        gene_header,
        GENE_EVID_FIELDS,
        "GENE_EVIDENCE",
    )

    evidence_ids = []

    for row in gene_evidence:

        evidence_id = txt(
            row[
                "evidence_id"
            ]
        )

        evidence_ids.append(
            evidence_id
        )

        key = (
            txt(
                row[
                    "condition"
                ]
            ),
            txt(
                row[
                    "validation_class"
                ]
            ),
            txt(
                row[
                    "feature_id"
                ]
            ),
        )

        if (
            not evidence_id
            or key
            not in candidate_map
            or not (
                present(
                    row[
                        "pmid"
                    ]
                )
                or
                present(
                    row[
                        "doi"
                    ]
                )
            )
            or not present(
                row[
                    "citation"
                ]
            )
            or not present(
                row[
                    "evidence_role"
                ]
            )
        ):
            raise AssertionError(
                (
                    "BAD_GENE_EVIDENCE",
                    evidence_id,
                    key,
                )
            )

        if (
            txt(
                row[
                    "supports_condition_relevance"
                ]
            )
            not in {
                "YES",
                "PARTIAL",
                "NO",
            }
            or
            txt(
                row[
                    "supports_directionality"
                ]
            )
            not in {
                "YES",
                "NO",
                "NOT_CLAIMED",
            }
        ):
            raise AssertionError(
                (
                    "BAD_GENE_EVIDENCE_FLAGS",
                    evidence_id,
                )
            )

    if (
        len(evidence_ids)
        != len(
            set(
                evidence_ids
            )
        )
    ):
        raise AssertionError(
            "DUPLICATE_GENE_EVIDENCE_ID"
        )

    condition_header, condition_evidence = read_tsv(
        root
        / "condition_literature_evidence.tsv"
    )

    require(
        condition_header,
        COND_EVID_FIELDS,
        "CONDITION_EVIDENCE",
    )

    evidence_ids = []

    for row in condition_evidence:

        evidence_id = txt(
            row[
                "evidence_id"
            ]
        )

        evidence_ids.append(
            evidence_id
        )

        if (
            not evidence_id
            or txt(
                row[
                    "condition"
                ]
            )
            not in CONDS
            or not (
                present(
                    row[
                        "pmid"
                    ]
                )
                or
                present(
                    row[
                        "doi"
                    ]
                )
            )
            or not present(
                row[
                    "citation"
                ]
            )
            or not present(
                row[
                    "evidence_claim"
                ]
            )
        ):
            raise AssertionError(
                (
                    "BAD_CONDITION_EVIDENCE",
                    evidence_id,
                )
            )

        if (
            txt(
                row[
                    "author_decision"
                ]
            )
            not in {
                "RETAIN_CONCORDANCE",
                "RETAIN_CONTEXT_ONLY",
                "REJECT",
            }
        ):
            raise AssertionError(
                (
                    "BAD_CONDITION_DECISION",
                    evidence_id,
                )
            )

    if (
        len(evidence_ids)
        != len(
            set(
                evidence_ids
            )
        )
    ):
        raise AssertionError(
            "DUPLICATE_CONDITION_EVIDENCE_ID"
        )

    return (
        decisions,
        gene_evidence,
        condition_evidence,
        candidate_map,
    )


def md(value):
    return (
        txt(value)
        .replace(
            "\n",
            " ",
        )
        .replace(
            "|",
            "\\|",
        )
    )


def final(
    root,
    candidate_rows,
    architecture_rows,
):

    (
        decisions,
        gene_evidence,
        condition_evidence,
        candidate_map,
    ) = manual(
        root,
        candidate_rows,
    )

    evidence_by_gene = defaultdict(
        list
    )

    for row in gene_evidence:

        key = (
            txt(
                row[
                    "condition"
                ]
            ),
            txt(
                row[
                    "validation_class"
                ]
            ),
            txt(
                row[
                    "feature_id"
                ]
            ),
        )

        evidence_by_gene[
            key
        ].append(
            row
        )

    retained_genes = []

    for decision_row in decisions:

        decision = (
            txt(
                decision_row[
                    "author_decision"
                ]
            )
            or "PENDING"
        )

        if decision not in {
            "RETAIN_HEADLINE",
            "RETAIN_SUPPORTING",
        }:
            continue

        key = (
            txt(
                decision_row[
                    "condition"
                ]
            ),
            txt(
                decision_row[
                    "validation_class"
                ]
            ),
            txt(
                decision_row[
                    "feature_id"
                ]
            ),
        )

        evidence = [
            row
            for row
            in evidence_by_gene[
                key
            ]
            if txt(
                row[
                    "supports_condition_relevance"
                ]
            )
            in {
                "YES",
                "PARTIAL",
            }
        ]

        if not evidence:
            raise AssertionError(
                (
                    "RETAINED_WITHOUT_EVIDENCE",
                    key,
                )
            )

        if (
            key[1]
            == "Low Stability"
            and any(
                txt(
                    row[
                        "supports_directionality"
                    ]
                )
                != "NOT_CLAIMED"
                for row
                in evidence
            )
        ):
            raise AssertionError(
                (
                    "LOW_STABILITY_DIRECTION_CLAIM",
                    key,
                )
            )

        base = candidate_map[
            key
        ]

        retained_genes.append(
            {
                "condition":
                    key[0],
                "validation_class":
                    key[1],
                "feature_id":
                    key[2],
                "candidate_rank":
                    base[
                        "candidate_rank"
                    ],
                "eggnog_preferred_name":
                    base[
                        "eggnog_preferred_name"
                    ],
                "eggnog_description":
                    base[
                        "eggnog_description"
                    ],
                "ssi":
                    base[
                        "ssi"
                    ],
                "author_decision":
                    decision,
                "manuscript_use":
                    txt(
                        decision_row[
                            "manuscript_use"
                        ]
                    ),
                "author_notes":
                    txt(
                        decision_row[
                            "author_notes"
                        ]
                    ),
                "literature_evidence_count":
                    len(
                        evidence
                    ),
                "evidence_ids":
                    ";".join(
                        txt(
                            row[
                                "evidence_id"
                            ]
                        )
                        for row
                        in evidence
                    ),
                "pmids":
                    ";".join(
                        dict.fromkeys(
                            txt(
                                row[
                                    "pmid"
                                ]
                            )
                            for row
                            in evidence
                            if present(
                                row[
                                    "pmid"
                                ]
                            )
                        )
                    ),
                "dois":
                    ";".join(
                        dict.fromkeys(
                            txt(
                                row[
                                    "doi"
                                ]
                            )
                            for row
                            in evidence
                            if present(
                                row[
                                    "doi"
                                ]
                            )
                        )
                    ),
                "evidence_roles":
                    " | ".join(
                        txt(
                            row[
                                "evidence_role"
                            ]
                        )
                        for row
                        in evidence
                    ),
                "directionality_support":
                    ";".join(
                        sorted(
                            set(
                                txt(
                                    row[
                                        "supports_directionality"
                                    ]
                                )
                                for row
                                in evidence
                            )
                        )
                    ),
            }
        )

    validated_gene_fields = (
        "condition",
        "validation_class",
        "feature_id",
        "candidate_rank",
        "eggnog_preferred_name",
        "eggnog_description",
        "ssi",
        "author_decision",
        "manuscript_use",
        "author_notes",
        "literature_evidence_count",
        "evidence_ids",
        "pmids",
        "dois",
        "evidence_roles",
        "directionality_support",
    )

    write_tsv(
        root
        / "validated_gene_shortlist.tsv",
        validated_gene_fields,
        retained_genes,
    )

    condition_evidence_by_condition = defaultdict(
        list
    )

    for row in condition_evidence:

        if txt(
            row[
                "author_decision"
            ]
        ) in {
            "RETAIN_CONCORDANCE",
            "RETAIN_CONTEXT_ONLY",
        }:

            condition_evidence_by_condition[
                txt(
                    row[
                        "condition"
                    ]
                )
            ].append(
                row
            )

    architecture_map = {
        row[
            "condition"
        ]:
            row
        for row in architecture_rows
    }

    validated_conditions = []

    for condition in CONDS:

        evidence = (
            condition_evidence_by_condition[
                condition
            ]
        )

        architecture_row = (
            architecture_map[
                condition
            ]
        )

        validated_conditions.append(
            {
                **architecture_row,
                "retained_evidence_count":
                    len(
                        evidence
                    ),
                "evidence_ids":
                    ";".join(
                        txt(
                            row[
                                "evidence_id"
                            ]
                        )
                        for row
                        in evidence
                    ),
                "pmids":
                    ";".join(
                        dict.fromkeys(
                            txt(
                                row[
                                    "pmid"
                                ]
                            )
                            for row
                            in evidence
                            if present(
                                row[
                                    "pmid"
                                ]
                            )
                        )
                    ),
                "dois":
                    ";".join(
                        dict.fromkeys(
                            txt(
                                row[
                                    "doi"
                                ]
                            )
                            for row
                            in evidence
                            if present(
                                row[
                                    "doi"
                                ]
                            )
                        )
                    ),
                "evidence_claims":
                    " | ".join(
                        txt(
                            row[
                                "evidence_claim"
                            ]
                        )
                        for row
                        in evidence
                    ),
                "evidence_decisions":
                    ";".join(
                        sorted(
                            set(
                                txt(
                                    row[
                                        "author_decision"
                                    ]
                                )
                                for row
                                in evidence
                            )
                        )
                    ),
            }
        )

    validated_condition_fields = (
        ARCH_FIELDS
        + (
            "retained_evidence_count",
            "evidence_ids",
            "pmids",
            "dois",
            "evidence_claims",
            "evidence_decisions",
        )
    )

    write_tsv(
        root
        / "validated_condition_architecture.tsv",
        validated_condition_fields,
        validated_conditions,
    )

    covered_conditions = sum(
        int(
            row[
                "retained_evidence_count"
            ]
        )
        > 0
        for row
        in validated_conditions
    )

    if (
        retained_genes
        and covered_conditions
        == 6
    ):
        status = "COMPLETE"

    elif (
        retained_genes
        or covered_conditions
    ):
        status = "PARTIAL"

    else:
        status = (
            "PENDING_MANUAL_LITERATURE_CURATION"
        )

    summary = [
        "# FSF biological validation",
        "",
        f"Status: **{status}**",
        "",
        (
            "Generated from locked FSF/annotation authorities "
            "plus explicit manual author decisions and "
            "literature-evidence registries. "
            "No literature interpretation is hidden in the builder."
        ),
        "",
        "## Condition-level architecture evidence",
        "",
        (
            "| Condition | Low Stability | Transitional | "
            "Stable | Highly Stable | Evidence | Claim(s) |"
        ),
        (
            "|---|---:|---:|---:|---:|---:|---|"
        ),
    ]

    for row in validated_conditions:

        summary.append(
            f"| {row['condition']} | "
            f"{float(row['low_stability_proportion']) * 100:.2f}% | "
            f"{float(row['transitional_proportion']) * 100:.2f}% | "
            f"{float(row['stable_proportion']) * 100:.2f}% | "
            f"{float(row['highly_stable_proportion']) * 100:.2f}% | "
            f"{row['retained_evidence_count']} | "
            f"{md(row['evidence_claims']) or 'Pending'} |"
        )

    summary.extend(
        [
            "",
            "## Validated gene shortlist",
            "",
        ]
    )

    if retained_genes:

        summary.extend(
            [
                (
                    "| Condition | FSF class | Feature | Gene | SSI | "
                    "Decision | PMID/DOI | Manuscript use |"
                ),
                (
                    "|---|---|---|---|---:|---|---|---|"
                ),
            ]
        )

        for row in retained_genes:

            identifiers = "; ".join(
                value
                for value in (
                    row[
                        "pmids"
                    ],
                    row[
                        "dois"
                    ],
                )
                if value
            )

            summary.append(
                f"| {row['condition']} | "
                f"{row['validation_class']} | "
                f"{row['feature_id']} | "
                f"{md(row['eggnog_preferred_name'])} | "
                f"{row['ssi']} | "
                f"{row['author_decision']} | "
                f"{md(identifiers)} | "
                f"{md(row['manuscript_use'])} |"
            )

    else:

        summary.append(
            "No gene has yet been author-retained "
            "with independent literature evidence."
        )

    summary.extend(
        [
            "",
            "## Guardrails",
            "",
            (
                "- Low Stability remains one "
                "non-directional FSF class."
            ),
            (
                "- Annotation prioritizes literature review "
                "but is not independent validation."
            ),
            (
                "- Retained genes require explicit author "
                "decisions and explicit PMID/DOI evidence."
            ),
            (
                "- Condition relevance does not automatically "
                "prove the exact FSF direction."
            ),
            "",
        ]
    )

    (
        root
        / "biological_validation_summary.md"
    ).write_text(
        "\n".join(
            summary
        ),
        encoding="utf-8",
    )

    return (
        status,
        len(
            retained_genes
        ),
        covered_conditions,
    )


def readme(root):

    (
        root
        / "README.md"
    ).write_text(
        """# Current FSF biological validation

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
""",
        encoding="utf-8",
    )


def manifest(
    root,
    fsf,
    annotation,
    status,
    retained_gene_count,
    covered_condition_count,
):

    items = [
        (
            "status",
            status,
        ),
        (
            "fsf_sha256",
            sha(
                fsf
            ),
        ),
        (
            "annotation_filename",
            annotation.name,
        ),
        (
            "annotation_sha256",
            sha(
                annotation
            ),
        ),
        (
            "fsf_rows",
            N_FSF,
        ),
        (
            "unique_features",
            N_FEATURES,
        ),
        (
            "validation_groups_after_low_stability_collapse",
            N_GROUPS,
        ),
        (
            "candidate_top5_rows",
            N_TOP5,
        ),
        (
            "retained_genes",
            retained_gene_count,
        ),
        (
            "conditions_with_retained_evidence",
            covered_condition_count,
        ),
    ]

    for name in (
        (
            "condition_architecture.tsv",
            "candidate_shortlist_top5.tsv",
        )
        + MANUAL
        + (
            "validated_gene_shortlist.tsv",
            "validated_condition_architecture.tsv",
            "biological_validation_summary.md",
            "README.md",
        )
    ):

        items.append(
            (
                f"sha256:{name}",
                sha(
                    root
                    / name
                ),
            )
        )

    write_tsv(
        root
        / "biological_validation_manifest.tsv",
        (
            "item",
            "value",
        ),
        (
            {
                "item":
                    key,
                "value":
                    value,
            }
            for key, value
            in items
        ),
    )


def generate(
    root,
    fsf,
    annotation,
    initialize=False,
):

    root.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        fsf_rows,
        annotation_map,
    ) = authorities(
        fsf,
        annotation,
    )

    architecture_rows = architecture(
        fsf_rows
    )

    candidate_rows = candidates(
        fsf_rows,
        annotation_map,
    )

    write_tsv(
        root
        / "condition_architecture.tsv",
        ARCH_FIELDS,
        architecture_rows,
    )

    write_tsv(
        root
        / "candidate_shortlist_top5.tsv",
        CAND_FIELDS,
        candidate_rows,
    )

    if initialize:
        init_manual(
            root,
            candidate_rows,
        )

    for name in MANUAL:

        if not (
            root
            / name
        ).is_file():

            raise AssertionError(
                f"Missing manual file {name}; "
                "run --initialize first."
            )

    readme(
        root
    )

    (
        status,
        retained_gene_count,
        covered_condition_count,
    ) = final(
        root,
        candidate_rows,
        architecture_rows,
    )

    manifest(
        root,
        fsf,
        annotation,
        status,
        retained_gene_count,
        covered_condition_count,
    )

    return (
        status,
        retained_gene_count,
        covered_condition_count,
    )


def verify(
    root,
    fsf,
    annotation,
):

    with tempfile.TemporaryDirectory(
        prefix="fsf-bioval-verify."
    ) as temp_directory:

        temporary_root = Path(
            temp_directory
        )

        for name in MANUAL:

            shutil.copy2(
                root
                / name,
                temporary_root
                / name,
            )

        generate(
            temporary_root,
            fsf,
            annotation,
            False,
        )

        for name in GENERATED:

            current = (
                root
                / name
            )

            rebuilt = (
                temporary_root
                / name
            )

            if (
                not current.is_file()
                or not filecmp.cmp(
                    current,
                    rebuilt,
                    shallow=False,
                )
            ):
                raise AssertionError(
                    (
                        "OUTPUT_DRIFT",
                        name,
                    )
                )

    print(
        "BIOLOGICAL_VALIDATION_REPRODUCIBILITY=PASS"
    )


def main():

    csv_limit()

    if (
        len(sys.argv)
        != 2
        or sys.argv[1]
        not in {
            "--initialize",
            "--build",
            "--verify",
        }
    ):
        raise SystemExit(
            "Usage: "
            "build_current_biological_validation.py "
            "--initialize|--build|--verify"
        )

    (
        _,
        fsf,
        annotation,
        root,
    ) = paths()

    mode = sys.argv[1]

    if mode == "--verify":

        verify(
            root,
            fsf,
            annotation,
        )

        return

    (
        status,
        retained_gene_count,
        covered_condition_count,
    ) = generate(
        root,
        fsf,
        annotation,
        mode
        == "--initialize",
    )

    print(
        "BIOLOGICAL_VALIDATION_ROOT="
        + str(
            root
        )
    )

    print(
        "BIOLOGICAL_VALIDATION_STATUS="
        + status
    )

    print(
        "RETAINED_GENE_COUNT="
        + str(
            retained_gene_count
        )
    )

    print(
        "CONDITIONS_WITH_RETAINED_EVIDENCE="
        + str(
            covered_condition_count
        )
    )

    print(
        "LOW_STABILITY_COLLAPSE=PASS"
    )

    print(
        "BIOLOGICAL_VALIDATION_BUILD=PASS"
    )


if __name__ == "__main__":
    main()
