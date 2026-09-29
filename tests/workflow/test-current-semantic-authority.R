#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out <- file.path(root, "results/current_fsf_v1/manuscript/semantic_authority")
external_root <- Sys.getenv("FSF_EXTERNAL_ARTIFACT_ROOT", unset = "")
stopifnot(nzchar(external_root))
resolver <- file.path(root, "scripts/current/build_current_external_artifact_manifests.R")
source(resolver)
external_root <- .fsf_external_root(external_root)
raw_locator <- "current_fsf_v1/manuscript/semantic_authority/raw"
raw <- .fsf_resolve_external_locator(raw_locator, external_root)
builder <- file.path(root, "scripts/current/semantic/build_current_semantic_authority.R")
go_path <- file.path(root, "results/current_fsf_v1/manuscript/source_data/figure5_go_source.tsv")
ke_path <- file.path(root, "results/current_fsf_v1/manuscript/source_data/figure5_kegg_source.tsv")

sha <- function(p) sub(" .*", "", system2("sha256sum", p, stdout = TRUE)[1])
rd <- function(p) read.delim(p, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
stopifnot(sha(go_path) == "0ded08c5f67abdcf6dbbc2d77bcd47530319ecdeeb7b6ad13cd402a65d5dcdc2")
stopifnot(sha(ke_path) == "780415dbe1a73d9246b16b58b58892b24f20d89eb1d3cf9f6458479d2c2f856b")

go <- rd(go_path); go <- go[go$significant & go$gene_set_family == "main_class", ]
ke <- rd(ke_path); ke <- ke[ke$significant & ke$gene_set_family == "main_class", ]
ga <- rd(file.path(out, "current_significant_go_annotated.tsv"))
ka <- rd(file.path(out, "current_significant_kegg_annotated.tsv"))
stopifnot(nrow(ga) == nrow(go), nrow(ka) == nrow(ke))

gkey <- function(x) paste(x$gene_set_id, x$go_term, sep = "\r")
kkey <- function(x) paste(x$gene_set_id, x$enrichment_type, x$kegg_term, sep = "\r")
stopifnot(identical(gkey(ga), gkey(go)), identical(kkey(ka), kkey(ke)))
stopifnot(isTRUE(all.equal(ga[names(go)], go, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(ka[names(ke)], ke, check.attributes = FALSE)))
go_identity <- function(x) paste(x$condition, x$stability_region, "GO", x$go_term, sep = "\r")
ke_identity <- function(x) paste(x$condition, x$stability_region, "KEGG", x$kegg_term, sep = "\r")
stopifnot(sum(duplicated(go_identity(ga))) == sum(duplicated(go_identity(go))))
stopifnot(sum(duplicated(ke_identity(ka))) == sum(duplicated(ke_identity(ke))))
stopifnot(setequal(go_identity(ga), go_identity(go)))
stopifnot(setequal(ke_identity(ka), ke_identity(ke)))

go_auth <- rd(file.path(out, "go_term_authority.tsv"))
ke_auth <- rd(file.path(out, "kegg_term_authority.tsv"))
stopifnot(all(nzchar(go_auth$go_name[go_auth$mapping_status != "unmapped"])))
stopifnot(all(nzchar(ke_auth$kegg_name[ke_auth$mapping_status == "mapped"])))
go_problem <- rd(file.path(out, "go_unmapped_or_problematic.tsv"))
ke_problem <- rd(file.path(out, "kegg_unmapped_or_problematic.tsv"))
stopifnot(setequal(go_problem$go_id, go_auth$go_id[go_auth$mapping_status != "mapped"]))
stopifnot(setequal(ke_problem$kegg_id, ke_auth$kegg_id[ke_auth$mapping_status != "mapped"]))

all_columns <- c(names(go_auth), names(ke_auth), names(ga), names(ka),
                 names(rd(file.path(out, "biological_semantic_review_packet.tsv"))))
forbidden <- c("theme", "interpretation", "adaptive", "causal", "mechanism",
               "fitness", "author_decision", "historical")
stopifnot(!any(vapply(forbidden, function(z)
  any(grepl(z, all_columns, ignore.case = TRUE)), logical(1))))

coverage <- rd(file.path(out, "semantic_coverage.tsv"))
stopifnot(coverage$total_unique_ids[coverage$source == "GO"] == length(unique(go$go_term)))
stopifnot(sum(coverage$total_unique_ids[coverage$source == "KEGG"]) ==
          nrow(unique(ke[c("enrichment_type", "kegg_term")])))

manifest <- rd(file.path(out, "semantic_source_manifest.tsv"))
stopifnot(all(vapply(manifest$local_artifact, .fsf_validate_external_locator, logical(1L))))
stopifnot(all(file.exists(vapply(manifest$local_artifact, .fsf_resolve_external_locator,
  character(1L), root = external_root))))
rejects <- c("/absolute/path", "C:/absolute/path", "\\\\server\\share\\file",
  "current/../escape", "..", "folder\\file")
stopifnot(all(vapply(rejects, function(locator)
  inherits(try(.fsf_validate_external_locator(locator), silent = TRUE), "try-error"), logical(1L))))
outside <- tempfile("fsf-outside-"); dir.create(outside)
stopifnot(inherits(try(.fsf_external_locator_from_path(outside, external_root), silent = TRUE), "try-error"))
unlink(outside, recursive = TRUE)

tmp <- tempfile("semantic-authority-rebuild-")
dir.create(tmp)
status <- system2("Rscript", c("--vanilla", builder, raw, tmp),
                  stdout = TRUE, stderr = TRUE,
                  env = paste0("FSF_EXTERNAL_ARTIFACT_ROOT=", shQuote(external_root)))
stopifnot(is.character(status), is.null(attr(status, "status")))
stopifnot(file.exists(file.path(tmp, "semantic_output_hashes.tsv")))
expected <- rd(file.path(out, "semantic_output_hashes.tsv"))
rebuilt <- rd(file.path(tmp, "semantic_output_hashes.tsv"))
stopifnot(identical(expected, rebuilt))
unlink(tmp, recursive = TRUE)

cat("current semantic authority tests: PASS\n")

