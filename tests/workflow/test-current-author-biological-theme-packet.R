#!/usr/bin/env Rscript
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
src<-file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
out<-file.path(root,"results/current_fsf_v1/manuscript/biological_theme_decisions")
builder<-file.path(root,"scripts/current/biological_theme_decisions/build_author_biological_theme_decision_packet.R")
rd<-function(p)read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
blank<-function(x)all(is.na(x)|!nzchar(x))
bhash<-c(candidate_theme_definitions.tsv="9dc65d51b39bc45d5885e000ae0e92ad7371e21732c17376527a5e1a2a450dfa",
 current_biological_evidence_inventory.tsv="bef472609c2790f1c6985493aa24149a65b7cbeb2b5f66ee7328bffefc600c58",
 biological_theme_evidence_trace.tsv="51be21df2d9484f83a4428ffe02a0d62e484415667eb102f5b785d9b5fe9ceee",
 candidate_term_theme_mapping.tsv="abab8afc65f4189d4c1b7220c495ba26ef793f84d57b1c67e67dfc76958e967b",
 go_redundancy_groups.tsv="6059c1e761f42c213222c5ea9540dc5a1ddee2155ba0c37905eee72b6faca639",
 candidate_theme_evidence.tsv="9c0a354ba0f98956ad844b35ee746e3dc54159cf39fb8d23c9d707816787c2a6",
 condition_region_biological_profile.tsv="62900e152b6f4b45529a5babe40396687dfa8a9e1fd6bcdc9531a9301d82a76c",
 within_condition_region_contrast.tsv="dad4cc1e90f58d7d731dac1d0535dd0eea6099776fe5e95429cfc7833be42f08",
 historical_theme_reconciliation.tsv="4a70702ddb57ae231247b4c0facd8f0ea31ee9bfc315410e0aaa1e7f09321958",
 biological_theme_author_review_packet.tsv="cb15514e28a556781f2a681b14f2c6f2ad56d48fd1a60fb0ead52b745c4e3f5f",
 biological_theme_review_summary.md="42335d4503e5d08e73f236fae208d941b6505403daa177fb6e76fb44090e0c9e")
stopifnot(all(vapply(file.path(src,names(bhash)),sha,"")==bhash))

defs<-rd(file.path(src,"candidate_theme_definitions.tsv"))
ev<-rd(file.path(src,"candidate_theme_evidence.tsv"))
mapping<-rd(file.path(src,"candidate_term_theme_mapping.tsv"))
trace<-rd(file.path(src,"biological_theme_evidence_trace.tsv"))
packet<-rd(file.path(out,"author_biological_theme_decision_packet.tsv"))
within<-rd(file.path(out,"within_condition_biological_review.tsv"))
global<-rd(file.path(out,"candidate_theme_global_review.tsv"))
stopifnot(nrow(global)==18L,setequal(global$candidate_theme,defs$candidate_theme))
required<-ev[ev$evidence_status=="REQUIRES_AUTHOR_REVIEW",]
required_key<-paste(required$condition,required$FSF_region,required$candidate_theme,sep="|")
stopifnot(all(required_key%in%packet$trace_key))
unit_key<-paste(ev$condition,ev$FSF_region,ev$candidate_theme,sep="|")
stopifnot(all(packet$trace_key%in%unit_key))
map_key<-paste(mapping$condition,mapping$FSF_region,mapping$candidate_theme,sep="|")
stopifnot(all(packet$trace_key%in%map_key))
trace_theme_key<-unique(unlist(lapply(seq_len(nrow(trace)),function(i){
 if(is.na(trace$candidate_themes[i])||!nzchar(trace$candidate_themes[i]))return(character())
 paste(trace$condition[i],trace$FSF_region[i],strsplit(trace$candidate_themes[i],"; ",fixed=TRUE)[[1]],sep="|")
})))
stopifnot(all(packet$trace_key%in%trace_theme_key))

for(x in list(packet,within,global))stopifnot(blank(x$author_decision),blank(x$author_notes))
options<-"RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER"
stopifnot(all(packet$recommended_action_options==options),all(global$recommended_action_options==options))
stopifnot(all(packet$review_priority%in%c("HIGH","MEDIUM","LOW")))
regions<-c("Low Stability","Transitional","Stable","Highly Stable")
stopifnot(all(packet$FSF_region%in%regions),all(within$region_A%in%regions),all(within$region_B%in%regions))

summary<-readLines(file.path(out,"author_biological_theme_review_summary.md"),warn=FALSE)
authored<-tolower(c(packet$candidate_theme,packet$within_condition_distinctiveness,
 packet$conflicting_or_ambiguous_evidence,within$review_interpretation_status,
 global$candidate_theme,global$region_specificity_summary,global$ambiguity_summary,summary))
prohibited<-c("adaptive program","protective","fitness","causal","mechanistic","functionally required","robust program")
stopifnot(!any(vapply(prohibited,function(z)any(grepl(z,authored,fixed=TRUE)),logical(1))))

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
protected<-c("results/current_fsf_v1/manuscript/semantic_authority",
 "results/current_fsf_v1/manuscript/biological_themes",
 "scripts/current/biological_themes")
source("tests/workflow/clean_release_hash_guard.R")
assert_paths_match_clean_manifest(protected)
review_test<-file.path(root,"tests","workflow","test-current-biological-theme-review.R")
stopifnot(file.exists(review_test),is.expression(parse(file=review_test)))

tmp<-tempfile("author-biological-theme-packet-");dir.create(tmp)
log<-system2("Rscript",c("--vanilla",builder,tmp),stdout=TRUE,stderr=TRUE)
rebuilt<-rd(file.path(tmp,"author_biological_theme_decision_hashes.tsv"))
frozen<-rd(file.path(out,"author_biological_theme_decision_hashes.tsv"))
stopifnot(identical(frozen,rebuilt))
unlink(tmp,recursive=TRUE)
cat("current author biological theme packet tests: PASS\n")

