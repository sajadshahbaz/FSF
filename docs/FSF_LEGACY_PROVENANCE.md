# FSF legacy provenance

## Authority boundary

The sole current scientific definition authority is
`docs/FSF_V1_CURRENT_SCIENTIFIC_LOCK.md`. Historical and legacy payloads are
not active FSF authority and cannot supply current terminology, thresholds,
classifications, numerical results, figures, tables, or conclusions.

## Preserved payload identities

- The historical archive payload contains exactly 101 files. Its public
  checksum manifest is `docs/provenance/historical-archive-files.sha256`.
- The legacy bundle contains exactly 57 files. Its public checksum manifest is
  `docs/provenance/legacy-bundle-files.sha256`.

The manifests preserve repository-relative path identities and SHA-256 values.
References use `historical-archive:<relative-path>` and
`legacy-bundle:<relative-path>` as provenance identifiers, not active filesystem
dependencies.

Current FSF package execution does not depend on either payload. Current
manuscript regeneration does not depend on either payload. Current workflows
must use `scripts/current/`, current governed inputs, and verified current
artifacts rather than resolving these historical identifiers.

## Deposit status

Durable archive locator: `UNASSIGNED_PENDING_PUBLIC_DEPOSIT`

No DOI, URL, release ID, or public archive identifier is currently claimed.
The checksum manifests and Git history preserve auditable identity and change
provenance while a durable public deposit is pending. The historical payloads
must remain recoverable until that locator is assigned and independently
verified.
