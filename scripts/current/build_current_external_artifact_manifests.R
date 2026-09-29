#!/usr/bin/env Rscript

# Portable external-artifact locator contract and deterministic manifest writer.
.fsf_external_root <- function(value = Sys.getenv("FSF_EXTERNAL_ARTIFACT_ROOT", unset = "")) {
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("FSF_EXTERNAL_ARTIFACT_ROOT is required for repository-external current artifacts.", call. = FALSE)
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

.fsf_validate_external_locator <- function(locator) {
  if (length(locator) != 1L || is.na(locator) || !nzchar(locator)) {
    stop("External artifact locator must be one nonempty relative path.", call. = FALSE)
  }
  if (grepl("^/", locator) || grepl("^[A-Za-z]:[/\\\\]", locator) ||
      grepl("^\\\\\\\\", locator) || grepl("\\\\", locator)) {
    stop("External artifact locator must use a relative forward-slash path.", call. = FALSE)
  }
  parts <- strsplit(locator, "/", fixed = TRUE)[[1L]]
  if (any(!nzchar(parts)) || any(parts %in% c(".", ".."))) {
    stop("External artifact locator contains an empty, '.' or '..' path component.", call. = FALSE)
  }
  TRUE
}

.fsf_path_inside_root <- function(path, root) {
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

.fsf_resolve_external_locator <- function(locator, root = .fsf_external_root(), must_work = TRUE) {
  .fsf_validate_external_locator(locator)
  root <- .fsf_external_root(root)
  resolved <- normalizePath(file.path(root, locator), winslash = "/", mustWork = must_work)
  if (!.fsf_path_inside_root(resolved, root)) {
    stop("Resolved external artifact escapes FSF_EXTERNAL_ARTIFACT_ROOT.", call. = FALSE)
  }
  resolved
}

.fsf_external_locator_from_path <- function(path, root = .fsf_external_root()) {
  root <- .fsf_external_root(root)
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!.fsf_path_inside_root(resolved, root) || identical(resolved, root)) {
    stop("External artifact path is not a file or directory beneath FSF_EXTERNAL_ARTIFACT_ROOT.", call. = FALSE)
  }
  locator <- substring(resolved, nchar(root) + 2L)
  .fsf_validate_external_locator(locator)
  locator
}

.fsf_derive_external_root <- function(path, required_suffix) {
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  .fsf_validate_external_locator(required_suffix)
  suffix <- paste0("/", required_suffix)
  if (!endsWith(resolved, suffix)) {
    stop("RAW_DIR is not under the canonical semantic raw suffix; set FSF_EXTERNAL_ARTIFACT_ROOT explicitly.", call. = FALSE)
  }
  root <- substr(resolved, 1L, nchar(resolved) - nchar(suffix))
  if (!nzchar(root)) root <- "/"
  root <- .fsf_external_root(root)
  if (!identical(.fsf_resolve_external_locator(required_suffix, root), resolved)) {
    stop("Unable to derive truthful portable RAW_DIR provenance.", call. = FALSE)
  }
  root
}

.fsf_sha256 <- function(path) {
  tool <- Sys.which("sha256sum")
  if (!nzchar(tool)) stop("sha256sum is required.", call. = FALSE)
  output <- system2(tool, path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 calculation failed.", call. = FALSE)
  sub("[[:space:]].*$", "", output[[1L]])
}

.fsf_read_header <- function(path) {
  strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
}

write_current_external_artifact_manifests <- function(
    root = normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    external_root = .fsf_external_root()) {
  external_root <- .fsf_external_root(external_root)
  artifacts <- data.frame(
    key = c("go", "kegg", "representative"),
    locator = c(
      "current_fsf_v1/manuscript/enrichment/go/go_enrichment_all.tsv",
      "current_fsf_v1/manuscript/enrichment/kegg/kegg_enrichment_all.tsv",
      "current_fsf_v1/representative_features/current_fsf_stable_signal_annotated_catalog.tsv"
    ),
    bytes = c(306249266, 35469610, 965459032),
    sha256 = c(
      "1b991db9a178bead03adf3fd96246d9b8e4b08e40828af300d29aca75b260411",
      "bcf428a508975b0ffb81324852243adcf8cc24f136de5ff63bb266a10150fa32",
      "2002148dcaa38a1b4e84c173d78e8aa34495a2dbb98f0fe5738b1ef0445fabdc"
    ), stringsAsFactors = FALSE
  )
  paths <- vapply(artifacts$locator, .fsf_resolve_external_locator, character(1L), root = external_root)
  expected_schema <- list(
    go = c("gene_set_id","gene_set_family","condition","stability_region","dominant_state","signal_class","go_term","query_size","background_size","query_overlap_count","background_term_count","query_fraction","background_fraction","enrichment_ratio","p_value","adjusted_p_value","significant","contributing_feature_ids","all_feature_count","annotated_feature_count"),
    kegg = c("gene_set_id","gene_set_family","condition","stability_region","dominant_state","signal_class","enrichment_type","kegg_term","query_size","background_size","query_overlap_count","background_term_count","query_fraction","background_fraction","enrichment_ratio","p_value","adjusted_p_value","significant","contributing_feature_ids","all_feature_count","annotated_feature_count"),
    representative = c("condition","feature_id","n_perturbations","n_up","n_down","n_const","p_up","p_down","p_const","dominant_state","ssi","stability_region","signal_class","stability_deviation","locus_tag","protein_id","protein_length","translation_status","eggnog_annotated","eggnog_seed_ortholog","eggnog_evalue","eggnog_score","eggnog_ogs","eggnog_max_annot_lvl","eggnog_cog_category","eggnog_description","eggnog_preferred_name","eggnog_go","eggnog_ec","eggnog_kegg_ko","eggnog_kegg_pathway","eggnog_kegg_module","eggnog_pfam","interpro_annotated","interpro_analyses","interpro_signature_accessions","interpro_signature_descriptions","interpro_ids","interpro_descriptions","interpro_go","interpro_pathways","n_interpro_rows","pfam_annotated","pfam_domains","pfam_accessions","pfam_descriptions","n_pfam_hits","annotation_source_count","annotation_status")
  )
  observed_bytes <- as.numeric(file.info(paths)$size)
  observed_hashes <- vapply(paths, .fsf_sha256, character(1L))
  observed_schema <- lapply(paths, .fsf_read_header)
  if (!identical(observed_bytes, artifacts$bytes) ||
      !identical(unname(observed_hashes), artifacts$sha256) ||
      !all(vapply(seq_along(paths), function(i) identical(observed_schema[[i]], expected_schema[[artifacts$key[[i]]]]), logical(1L)))) {
    stop("Locked external artifact size, SHA-256, or schema validation failed.", call. = FALSE)
  }
  external_manifest <- data.frame(
    relative_repository_path = c("enrichment/go/go_enrichment_all.tsv", "enrichment/kegg/kegg_enrichment_all.tsv"),
    external_artifact_path = artifacts$locator[1:2], bytes = artifacts$bytes[1:2], sha256 = artifacts$sha256[1:2],
    artifact_role = c("complete GO enrichment result authority", "complete combined KO/PATHWAY enrichment result authority"),
    generator_script = "scripts/current/build_current_manuscript_enrichment_outputs.R",
    package_tag = "fsf-r-v0.1.0", package_commit = "ddb347a19e9877032a49aec2ebfeacc551b0b315",
    scientific_lock_sha256 = "0c7a1462d210144c726041952cb4dfd7a59b386a7136d161090fa443d1d3b9d1",
    stringsAsFactors = FALSE
  )
  large_manifest <- data.frame(
    artifact_name = "current_fsf_stable_signal_annotated_catalog.tsv",
    artifact_role = "current Stable/Highly Stable FSF candidate catalog joined to authoritative biological annotation",
    bytes = artifacts$bytes[3], sha256 = artifacts$sha256[3], git_policy = "EXTERNAL_LARGE_ARTIFACT_NOT_GIT_BLOB",
    scientific_identity = "current_fsf_stable_signal_annotated_catalog.tsv + SHA-256 + validated 49-column schema and frozen-authority provenance",
    local_persistent_locator = artifacts$locator[3],
    publication_status = "LOCAL_PERSISTENT_COPY_VALIDATED; PUBLIC_ARCHIVE_NOT_YET_ASSIGNED",
    stringsAsFactors = FALSE
  )
  external_out <- file.path(root, "results/current_fsf_v1/manuscript/audit/external_artifacts.tsv")
  large_out <- file.path(root, "results/current_fsf_v1/representative_features/LARGE_ARTIFACT_MANIFEST.tsv")
  write.table(external_manifest, external_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  write.table(large_manifest, large_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(c(external_out, large_out))
}

if (sys.nframe() == 0L) {
  write_current_external_artifact_manifests()
  cat("current external artifact manifests generated: PASS\n")
}
