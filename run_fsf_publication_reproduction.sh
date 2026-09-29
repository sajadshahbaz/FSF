#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$ROOT"

ANNOTATION=${FSF_ANNOTATION_MASTER:-$ROOT/governed_artifacts/FSF_v1_annotation_master.tsv}
EXTERNAL=${FSF_EXTERNAL_ARTIFACT_ROOT:-$ROOT/governed_artifacts}
REAL_EFFECTS="$ROOT/data/processed/effect_estimates_real.tsv"
SYNTHETIC_EFFECTS="$ROOT/data/processed/effect_estimates.tsv"
RESULTS="$ROOT/results/current_fsf_v1"
MANUSCRIPT="$RESULTS/manuscript"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/fsf-publication-run.XXXXXX")
FSF_LIB="$WORK/lib"
mkdir -p "$WORK/build" "$FSF_LIB"

export FSF_EXTERNAL_ARTIFACT_ROOT="$EXTERNAL"
export FSF_ANNOTATION_MASTER="$ANNOTATION"
export FSF_REPO_ROOT="$ROOT"

sha() { sha256sum "$1" | awk '{print $1}'; }
require_hash() {
  [ -f "$1" ]
  [ "$(sha "$1")" = "$2" ]
}

require_hash "$REAL_EFFECTS" afed2cdbc7794c62845e32ad81ce673d5e2e8697ef94e5c2f8e52595a322cce4
require_hash "$SYNTHETIC_EFFECTS" 7e561f3ece8e4ae4656579bf82d627b5e836710057a13d9e280539d98d351aa0
require_hash "$ANNOTATION" 25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584

(
  cd "$WORK/build"
  R CMD build --no-build-vignettes --no-manual "$ROOT"
)
PKG=$(find "$WORK/build" -maxdepth 1 -type f -name 'FSF_*.tar.gz' -print | head -n 1)
[ -f "$PKG" ]
R CMD INSTALL --library="$FSF_LIB" "$PKG"

Rscript --vanilla - "$FSF_LIB" <<'RS'
lib <- normalizePath(commandArgs(trailingOnly=TRUE)[1], mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
library(FSF,lib.loc=lib)
p <- normalizePath(find.package("FSF"),mustWork=TRUE)
stopifnot(startsWith(p,paste0(lib,"/")))
cat("INSTALLED_FSF_PATH=",p,"\n",sep="")
cat("INSTALLED_FSF_VERSION=",as.character(packageVersion("FSF")),"\n",sep="")
RS

Rscript --vanilla - "$FSF_LIB" <<'RS'
lib <- normalizePath(commandArgs(trailingOnly=TRUE)[1], mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
library(testthat)
library(FSF,lib.loc=lib)
stopifnot(startsWith(normalizePath(find.package("FSF")),paste0(lib,"/")))
test_dir("tests/testthat",reporter="summary",stop_on_failure=TRUE,stop_on_warning=TRUE)
RS

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/reproduce_current_real_data_fsf.R \
  "$ROOT" "$REAL_EFFECTS" "$WORK/core"
mkdir -p "$RESULTS"
cp -p "$WORK/core"/current_fsf_* "$RESULTS/"

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/build_current_manuscript_tables.R \
  "$ANNOTATION" "$SYNTHETIC_EFFECTS" "$WORK/manuscript"
for d in tables enrichment benchmarks source_data; do
  rm -rf -- "$MANUSCRIPT/$d"
  cp -a -- "$WORK/manuscript/$d" "$MANUSCRIPT/$d"
done

# Audit contains both regenerated outputs and frozen release-governance files.
# Overlay regenerated audit outputs without deleting the governed release guards.
mkdir -p "$MANUSCRIPT/audit"
cp -a -- "$WORK/manuscript/audit/." "$MANUSCRIPT/audit/"

R_LIBS="$FSF_LIB" Rscript --vanilla - "$ANNOTATION" "$RESULTS/representative_features" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
source("scripts/current/build_current_representative_features.R")
out <- build_current_representative_features(args[1])
dest <- args[2]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
external <- "current_fsf_stable_signal_annotated_catalog"
for (nm in setdiff(names(out),external)) {
  readr::write_tsv(out[[nm]],file.path(dest,paste0(nm,".tsv")),na="NA")
}
RS

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/build_current_external_artifact_manifests.R

R_LIBS="$FSF_LIB" Rscript --vanilla - "$WORK/review-template" <<'RS'
root <- commandArgs(trailingOnly=TRUE)[1]
source("scripts/current/build_current_biological_review_package.R")
build_current_biological_review_package(root)
audit_current_biological_review_package(root,review_state="template")
RS
for f in current_go_evidence.tsv current_kegg_evidence.tsv uv_specific_review.tsv review_summary.md; do
  cp -p "$WORK/review-template/$f" "$MANUSCRIPT/review/$f"
done

REVIEW_AUTHORITY="$EXTERNAL/current_fsf_v1/manuscript/review_authority"

test "$(sha256sum "$REVIEW_AUTHORITY/deferred_profile_term_evidence.tsv" | awk '{print $1}')" = \
  "e03c94ef535209fa76a587694c536381eb73a69472d693ef61b9c2ae48b613be"

test "$(sha256sum "$REVIEW_AUTHORITY/final_biological_interpretation_summary.md" | awk '{print $1}')" = \
  "eeeeaa1180bac85e2de5417b43248c37c84280199e6d44a830204dc3979ed95e"

cp -p \
  "$REVIEW_AUTHORITY/deferred_profile_term_evidence.tsv" \
  "$MANUSCRIPT/review/deferred_profile_term_evidence.tsv"

cp -p \
  "$REVIEW_AUTHORITY/final_biological_interpretation_summary.md" \
  "$MANUSCRIPT/review/final_biological_interpretation_summary.md"

cp -p "$MANUSCRIPT/review/author_decision_packet_LOCKED.tsv"   "$MANUSCRIPT/review/author_decision_packet.tsv"
Rscript --vanilla - "$MANUSCRIPT/review/author_decision_packet.tsv"   "$MANUSCRIPT/review/author_decision_lock_summary.md" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
x <- read.delim(args[1],check.names=FALSE,quote="",stringsAsFactors=FALSE)
levels <- c("RETAIN","RETAIN_WITH_REVISED_WORDING","PARTIALLY_RETAIN","DROP","REVIEW_FURTHER")
counts <- table(factor(x$author_decision,levels=levels))
unresolved <- x$decision_id[x$author_decision=="REVIEW_FURTHER"]
rules <- c(
 "FSF method is frozen.",
 "No additional threshold redesign.",
 "No new dataset or analytical branch.",
 "Absence of significant enrichment may be reported descriptively.",
 "Absence of enrichment is not evidence of absence of biology.",
 "Broad/manual catch-all themes are not manuscript biological claims.",
 "Partial support means only the currently supported component is retained.",
 "Causal/adaptive/fitness interpretations are not inferred from architecture or enrichment alone.",
 "UV is Transitional-dominant under the current authority.",
 "Representative-feature discovery is closed; only rationale revision is allowed."
)
text <- c(
 "# Author biological decision lock summary","",
 "## Decision counts","",
 paste0("- ",levels,": ",as.integer(counts)),"",
 "## Intentionally unresolved profiles","",
 paste0("- ",unresolved),"",
 "## Locked interpretation rules","",
 paste0(seq_along(rules),". ",rules)
)
writeLines(text,args[2],useBytes=TRUE)
RS

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/semantic/build_current_semantic_authority.R \
  "$EXTERNAL/current_fsf_v1/manuscript/semantic_authority/raw" \
  "$MANUSCRIPT/semantic_authority"

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/biological_themes/build_current_biological_theme_review.R \
  "$MANUSCRIPT/biological_themes"

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/biological_theme_decisions/build_author_biological_theme_decision_packet.R \
  "$WORK/theme-decisions"

for f in author_biological_theme_decision_packet.tsv          author_biological_theme_decision_hashes.tsv          author_biological_theme_review_summary.md; do
  cp -p "$WORK/theme-decisions/$f" "$MANUSCRIPT/biological_theme_decisions/$f"
done

python3 scripts/current/build_current_biological_validation.py --build

mkdir -p "$MANUSCRIPT/figures/new"
cp -p "$MANUSCRIPT/figures/main/Figure1.png"   "$MANUSCRIPT/figures/new/Figure_1.png"
cp -p "$MANUSCRIPT/figures/main/Figure7.png"   "$MANUSCRIPT/figures/new/Figure_7.png"

R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/figures/build_author_reviewed_biological_theme_panels.R
R_LIBS="$FSF_LIB" Rscript --vanilla \
  scripts/current/figures/build_current_manuscript_figures.R

STAGE12_DEV="$WORK/stage12-dev"
mkdir -p "$STAGE12_DEV/governed_historical_authority"
cp -p \
  "$EXTERNAL/current_fsf_v1/manuscript/historical_authority/Supplementary_Table_S1_metadata.csv" \
  "$STAGE12_DEV/governed_historical_authority/Supplementary_Table_S1_metadata.csv"
rm -rf -- "$MANUSCRIPT/final_tables_candidate"
python3 scripts/current/build_final_table_package.py \
  "$ROOT" "$STAGE12_DEV" "$MANUSCRIPT/final_tables_candidate"

workflow_count=0
while IFS= read -r test_file; do
  [ -n "$test_file" ] || continue
  R_LIBS="$FSF_LIB" Rscript --vanilla "$test_file"
  workflow_count=$((workflow_count + 1))
done < <(find tests/workflow -maxdepth 1 -type f -name '*.R' -printf '%p\n' | LC_ALL=C sort)

printf 'WORKFLOW_TEST_FILE_COUNT=%s\n' "$workflow_count"
printf 'PACKAGE_TESTS=PASS\n'
printf 'WORKFLOW_TESTS=PASS\n'
printf 'PUBLICATION_WORKFLOW=PASS\n'
printf 'FSF_PACKAGE_LIBRARY_PATH_USED=%s\n' "$FSF_LIB"
printf 'RUN_WORK_DIRECTORY=%s\n' "$WORK"
