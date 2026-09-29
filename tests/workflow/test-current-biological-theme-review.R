#!/usr/bin/env Rscript
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
out<-file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
sem<-file.path(root,"results/current_fsf_v1/manuscript/semantic_authority")
builder<-file.path(root,"scripts/current/biological_themes/build_current_biological_theme_review.R")
source(file.path(root,"scripts/current/build_current_external_artifact_manifests.R"))
external_root<-.fsf_external_root()
rd<-function(p)read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
historical_authority<-.fsf_resolve_external_locator(
 "current_fsf_v1/manuscript/historical_authority/Table_SignalClass_Biology_CURATED.tsv",
 external_root)
stopifnot(file.exists(historical_authority),
 sha(historical_authority)=="c4704fb6d59e0185919ea75ab9b3ee42a004f64da39a264ac0f5c169febf2529")
blank<-function(x)all(is.na(x)|!nzchar(x))
snapshot<-function(paths){
 files<-unlist(lapply(paths,function(path){
  if(!dir.exists(path))return(character())
  list.files(path,all.files=TRUE,recursive=TRUE,full.names=TRUE,no..=TRUE)
 }),use.names=FALSE)
 if(!length(files))return(setNames(character(),character()))
 info<-file.info(files)
 files<-sort(normalizePath(files[!is.na(info[["isdir"]])&!info[["isdir"]]],
  winslash="/",mustWork=TRUE))
 relative<-substring(files,nchar(root)+2L)
 setNames(vapply(files,function(path)sha(shQuote(path)),""),relative)
}
protected_scope<-file.path(root,c("manuscript","docs"))
protected_before<-snapshot(protected_scope)
expected<-c(go_term_authority.tsv="608020c5f4c95263a34749b54b23a797c621184a27979dd59652506fcee1a218",
 kegg_term_authority.tsv="9652ad9bf7a89cf6cbb7f73da5cf99e8808241be284531bd47f85d21eef50671",
 current_significant_go_annotated.tsv="6bcf373e531b9dbfd8e1e89dc25ba77e2b2407c9c5dc3741408da281571786dd",
 current_significant_kegg_annotated.tsv="6bf9666965f84c7f183e75331f90cb878f2eabe0b8caf7adff6d520dfcb39735",
 biological_semantic_review_packet.tsv="c3b720079f405a8f9b20a4704435ae09e3a4e0887b559106bd4b85f996c42403")
stopifnot(all(vapply(file.path(sem,names(expected)),sha,"")==expected))
semantic_manifest<-rd(file.path(sem,"semantic_source_manifest.tsv"))
stopifnot(all(vapply(semantic_manifest$local_artifact,.fsf_validate_external_locator,logical(1L))))
obo_row<-semantic_manifest$resource_type=="go-basic OBO"
stopifnot(sum(obo_row)==1L)
obo_path<-.fsf_resolve_external_locator(semantic_manifest$local_artifact[obo_row],external_root)
stopifnot(sha(obo_path)==semantic_manifest$sha256[obo_row])

source_packet<-rd(file.path(sem,"biological_semantic_review_packet.tsv"))
trace<-rd(file.path(out,"biological_theme_evidence_trace.tsv"))
stopifnot(nrow(trace)==nrow(source_packet))
stopifnot(isTRUE(all.equal(trace[names(source_packet)],source_packet,check.attributes=FALSE)))
stopifnot(!anyDuplicated(trace$evidence_row_id))

mapping<-rd(file.path(out,"candidate_term_theme_mapping.tsv")
)
ga<-rd(file.path(sem,"go_term_authority.tsv"));ka<-rd(file.path(sem,"kegg_term_authority.tsv"))
stopifnot(all(mapping$term_id[mapping$source=="GO"]%in%ga$go_id))
km<-mapping[mapping$source=="KEGG",]
stopifnot(all(paste(km$term_id,ifelse(grepl("^K[0-9]",km$term_id),"KO","PATHWAY"))%in%
 paste(ka$kegg_id,ka$kegg_identifier_type)))
allowed_basis<-c("ONTOLOGY_RELATION","OFFICIAL_NAME_MATCH","KEGG_PATHWAY_IDENTITY",
 "HISTORICAL_HYPOTHESIS_SUPPORTED","CURRENT_EVIDENCE_DISCOVERY")
stopifnot(all(mapping$mapping_basis%in%allowed_basis))

files_with_author<-c("candidate_term_theme_mapping.tsv","candidate_theme_evidence.tsv",
 "within_condition_region_contrast.tsv","historical_theme_reconciliation.tsv",
 "biological_theme_author_review_packet.tsv")
for(f in files_with_author){x<-rd(file.path(out,f));stopifnot(blank(x$author_decision),blank(x$author_notes))}
ev<-rd(file.path(out,"candidate_theme_evidence.tsv"))
allowed_status<-c("STRONG_MULTI_SOURCE_SUPPORT","GO_SUPPORTED","KEGG_SUPPORTED","PARTIAL_SUPPORT",
 "WEAK_OR_BROAD_SUPPORT","NO_CURRENT_SUPPORT","REQUIRES_AUTHOR_REVIEW")
stopifnot(all(ev$evidence_status%in%allowed_status))
review<-rd(file.path(out,"biological_theme_author_review_packet.tsv"))
stopifnot(all(review$recommended_action_options==
 "RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER"))

regions<-c("Low Stability","Transitional","Stable","Highly Stable")
conditions<-c("DES","GAM","HT","LT","OSM","UV")
stopifnot(setequal(unique(trace$condition),conditions))
stopifnot(all(unique(trace$FSF_region)%in%regions))
for(f in c("candidate_term_theme_mapping.tsv","candidate_theme_evidence.tsv",
 "condition_region_biological_profile.tsv","within_condition_region_contrast.tsv")){
 x<-rd(file.path(out,f));cols<-intersect(names(x),c("FSF_region","region_A","region_B"))
 stopifnot(all(unlist(x[cols])%in%regions))
}
stopifnot(!any(tolower(unlist(ev[c("FSF_region")]))%in%c("instable","weakly stable","weakly_stable")))

defs<-rd(file.path(out,"candidate_theme_definitions.tsv"))
profile<-rd(file.path(out,"condition_region_biological_profile.tsv"))
summary_text<-readLines(file.path(out,"biological_theme_review_summary.md"),warn=FALSE)
authored<-tolower(c(defs$candidate_theme,profile$candidate_themes,profile$strongest_themes,
 profile$GO_only_themes,profile$KEGG_only_themes,profile$cross_source_themes,summary_text))
prohibited<-c("adaptive program","protective","fitness","causal","mechanistic","functionally required","robust program")
stopifnot(!any(vapply(prohibited,function(z)any(grepl(z,authored,fixed=TRUE)),logical(1))))

groups<-rd(file.path(out,"go_redundancy_groups.tsv"))
go<-rd(file.path(sem,"current_significant_go_annotated.tsv"))
members<-unlist(strsplit(groups$member_term_ids,";",fixed=TRUE))
stopifnot(setequal(members,unique(go$go_term)))

figure5_paths<-file.path(root,"results/current_fsf_v1/manuscript/figures/main",
 c("Figure5.pdf","Figure5.png"))
figure7_paths<-file.path(root,"results/current_fsf_v1/manuscript/figures/main",
 c("Figure7.pdf","Figure7.png"))
figure7_source<-file.path(root,
 "results/current_fsf_v1/manuscript/figures/new/Figure_7.png")
stopifnot(all(file.exists(c(figure5_paths,figure7_paths,figure7_source))))
figure7_png<-figure7_paths[basename(figure7_paths)=="Figure7.png"]
protected_figure7_hash<-"6ab8be9958a8073303661814fa7def0de9e1a5d9ece745adb23dc1bf377771a0"
stopifnot(system2("cmp",c("--silent",figure7_png,figure7_source))==0L,
 sha(figure7_png)==sha(figure7_source),
 sha(figure7_png)==protected_figure7_hash,
 sha(figure7_source)==protected_figure7_hash)

tmp<-tempfile("biological-theme-review-");dir.create(tmp)
log<-system2("Rscript",c("--vanilla",builder,tmp),stdout=TRUE,stderr=TRUE)
protected_after<-snapshot(protected_scope)
stopifnot(identical(protected_before,protected_after))
rebuilt<-rd(file.path(tmp,"biological_theme_output_hashes.tsv"))
frozen<-rd(file.path(out,"biological_theme_output_hashes.tsv"))
stopifnot(identical(frozen,rebuilt))
unlink(tmp,recursive=TRUE)
cat("current biological theme review tests: PASS\n")

