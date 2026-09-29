clean_release_sha256 <- function(path) {
  if (!file.exists(path)) {
    stop("Protected clean-release path is missing: ", path)
  }

  value <- system2(
    "sha256sum",
    path,
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(value, "status")

  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for: ", path)
  }

  sub(" .*", "", value[1])
}

assert_paths_match_clean_manifest <- function(paths) {
  root <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )

  manifest_path <- file.path(
    root,
    "results",
    "current_fsf_v1",
    "manuscript",
    "audit",
    "current_release_test_guard.sha256.tsv"
  )

  if (!file.exists(manifest_path)) {
    stop("Clean-release test guard manifest is missing")
  }

  manifest <- read.delim(
    manifest_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!identical(
    names(manifest),
    c("sha256", "relative_path")
  )) {
    stop("Unexpected clean-release test guard schema")
  }

  expand_path <- function(rel) {
    p <- file.path(root, rel)

    if (!file.exists(p)) {
      stop("Protected path is missing: ", rel)
    }

    if (dir.exists(p)) {
      files <- list.files(
        p,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
      )

      files[
        file.info(files)$isdir %in% FALSE
      ]
    } else {
      p
    }
  }

  files <- unique(
    unlist(
      lapply(paths, expand_path),
      use.names = FALSE
    )
  )

  if (!length(files)) {
    stop("Protected path set resolved to zero files")
  }

  canonical_root <- paste0(root, "/")

  relative <- vapply(
    normalizePath(
      files,
      winslash = "/",
      mustWork = TRUE
    ),
    function(p) {
      if (!startsWith(p, canonical_root)) {
        stop("Protected path escapes clean root: ", p)
      }

      substring(
        p,
        nchar(canonical_root) + 1L
      )
    },
    character(1)
  )

  idx <- match(
    relative,
    manifest$relative_path
  )

  if (anyNA(idx)) {
    stop(
      "Protected file absent from clean-release test guard: ",
      paste(
        relative[is.na(idx)],
        collapse = "; "
      )
    )
  }

  actual <- vapply(
    files,
    clean_release_sha256,
    character(1)
  )

  expected <- manifest$sha256[idx]

  mismatch <- which(
    actual != expected
  )

  if (length(mismatch)) {
    stop(
      "Protected clean-release file changed: ",
      paste(
        relative[mismatch],
        collapse = "; "
      )
    )
  }

  invisible(TRUE)
}
