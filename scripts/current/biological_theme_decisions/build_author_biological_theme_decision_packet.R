#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
a<-commandArgs(trailingOnly=TRUE)
out<-if(length(a))a[1] else file.path(root,"results/current_fsf_v1/manuscript/biological_theme_decisions")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src<-file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
rd<-function(f)read.delim(file.path(src,f),check.names=FALSE,quote="",stringsAsFactors=FALSE)
wr<-function(x,p)write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="")
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
conditions<-c("DES","GAM","HT","LT","OSM","UV");regions<-c("Low Stability","Transitional","Stable","Highly Stable")
defs<-rd("candidate_theme_definitions.tsv");ev<-rd("candidate_theme_evidence.tsv")
profile<-rd("condition_region_biological_profile.tsv");contrast<-rd("within_condition_region_contrast.tsv")
recon<-rd("historical_theme_reconciliation.tsv");base_packet<-rd("biological_theme_author_review_packet.tsv")
trace<-rd("biological_theme_evidence_trace.tsv");mapping<-rd("candidate_term_theme_mapping.tsv")
units<-ev[ev$evidence_status!="NO_CURRENT_SUPPORT",]
units<-units[order(match(units$condition,conditions),match(units$FSF_region,regions),units$candidate_theme),]
collapse_unique<-function(x,sep="; ")paste(sort(unique(x[!is.na(x)&nzchar(x)])),collapse=sep)
strongest<-function(g,k){z<-c(g,k);z<-z[is.finite(z)];if(length(z))min(z)else NA_real_}
fmt<-function(x)ifelse(is.na(x),"",format(x,digits=4,scientific=TRUE))
trace_key<-paste(units$condition,units$FSF_region,units$candidate_theme,sep="|")
support_regions<-lapply(seq_len(nrow(units)),function(i){
 z<-units[units$condition==units$condition[i]&units$candidate_theme==units$candidate_theme[i],]
 z$FSF_region
})
distinct<-vapply(seq_len(nrow(units)),function(i){
 rs<-unique(support_regions[[i]]);pop<-unique(ev$FSF_region[ev$condition==units$condition[i]&ev$evidence_status!="NO_CURRENT_SUPPORT"])
 if(length(pop)==1)"SINGLE_POPULATED_REGION" else if(length(rs)==1)"REGION_UNIQUE_CANDIDATE" else paste0("SHARED_ACROSS_",length(rs),"_REGIONS")
},"")
examples<-vapply(seq_len(nrow(units)),function(i){
 z<-mapping[mapping$condition==units$condition[i]&mapping$FSF_region==units$FSF_region[i]&mapping$candidate_theme==units$candidate_theme[i],]
 z<-z[order(z$adjusted_p,z$source,z$term_id),]
 paste(head(paste(z$source,z$term_id,z$official_term_name,sep=":"),5),collapse="; ")
},"")
ambig<-vapply(seq_len(nrow(units)),function(i){
 z<-mapping[mapping$condition==units$condition[i]&mapping$FSF_region==units$FSF_region[i]&mapping$candidate_theme==units$candidate_theme[i],]
 q<-character()
 if(any(z$evidence_status=="REQUIRES_AUTHOR_REVIEW"))q<-c(q,"Obsolete GO term evidence present")
 if(units$candidate_theme[i]=="Broad regulatory and signaling processes")q<-c(q,"Broad candidate definition")
 if(units$historical_relation[i]=="HISTORICAL_CANDIDATE")q<-c(q,"Historical hypothesis; author review pending")
 if(!length(q))"None identified by deterministic audit" else paste(q,collapse="; ")
},"")
pminv<-mapply(strongest,units$strongest_GO_adjusted_p,units$strongest_KEGG_adjusted_p)
priority<-ifelse(units$cross_source_support=="YES","HIGH",
 ifelse(distinct=="REGION_UNIQUE_CANDIDATE"&pminv<=0.01,"HIGH",
 ifelse(units$evidence_status%in%c("GO_SUPPORTED","KEGG_SUPPORTED")&pminv<=0.05&
 units$candidate_theme!="Broad regulatory and signaling processes","MEDIUM","LOW")))
packet<-data.frame(decision_id=sprintf("AUTHOR_BIO_THEME_%04d",seq_len(nrow(units))),
 condition=units$condition,FSF_region=units$FSF_region,candidate_theme=units$candidate_theme,
 GO_support=paste0(units$GO_term_count_raw," terms / ",units$GO_redundancy_group_count," redundancy groups"),
 KEGG_support=paste0(units$KEGG_term_count," terms"),cross_source_support=units$cross_source_support,
 GO_redundancy_group_count=units$GO_redundancy_group_count,
 strongest_GO_adjusted_p=units$strongest_GO_adjusted_p,strongest_KEGG_adjusted_p=units$strongest_KEGG_adjusted_p,
 evidence_status=units$evidence_status,historical_relation=units$historical_relation,
 within_condition_distinctiveness=distinct,supporting_term_examples=examples,
 conflicting_or_ambiguous_evidence=ambig,
 recommended_action_options="RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER",
 review_priority=priority,trace_key=trace_key,author_decision="",author_notes="",stringsAsFactors=FALSE)
wr(packet,file.path(out,"author_biological_theme_decision_packet.tsv"))

distinguish<-function(c,r,themes){
 if(is.na(themes)||!nzchar(themes))return("")
 th<-strsplit(themes,"; ",fixed=TRUE)[[1]]
 z<-units[units$condition==c&units$FSF_region==r&units$candidate_theme%in%th,]
 if(!nrow(z))return("")
 p<-mapply(strongest,z$strongest_GO_adjusted_p,z$strongest_KEGG_adjusted_p);o<-order(p,z$candidate_theme)
 paste(head(paste0(z$candidate_theme[o]," [p_adj=",fmt(p[o]),"]"),5),collapse="; ")
}
wc<-contrast[contrast$condition%in%c("DES","GAM","LT","UV"),]
within<-data.frame(condition=wc$condition,region_A=wc$region_A,region_B=wc$region_B,
 Jaccard_similarity=wc$profile_similarity_measure,shared_candidate_themes=wc$themes_shared,
 unique_themes_A=wc$themes_unique_A,unique_themes_B=wc$themes_unique_B,
 strongest_distinguishing_evidence_A=mapply(distinguish,wc$condition,wc$region_A,wc$themes_unique_A),
 strongest_distinguishing_evidence_B=mapply(distinguish,wc$condition,wc$region_B,wc$themes_unique_B),
 review_interpretation_status="DESCRIPTIVE_REGION_COMPARISON_REQUIRES_AUTHOR_REVIEW",
 author_decision="",author_notes="",stringsAsFactors=FALSE)
wr(within,file.path(out,"within_condition_biological_review.tsv"))

split_values<-function(x)unique(unlist(strsplit(x[nzchar(x)],"; ",fixed=TRUE)))
global<-do.call(rbind,lapply(seq_len(nrow(defs)),function(i){
 t<-defs$candidate_theme[i];z<-units[units$candidate_theme==t,]
 rr<-recon[grepl(t,recon$current_matching_evidence,fixed=TRUE),]
 p<-mapply(strongest,z$strongest_GO_adjusted_p,z$strongest_KEGG_adjusted_p)
 data.frame(candidate_theme=t,historical_or_new=ifelse(defs$route[i]=="HISTORICAL_HYPOTHESIS","HISTORICAL","NEW_CURRENT"),
  conditions_supported=collapse_unique(z$condition),regions_supported=collapse_unique(z$FSF_region),
  GO_support_total=sum(z$GO_term_count_raw),KEGG_support_total=sum(z$KEGG_term_count),
  cross_source_support_count=sum(z$cross_source_support=="YES"),
  region_specificity_summary=paste0(sum(vapply(seq_len(nrow(z)),function(j){
   sum(units$condition==z$condition[j]&units$candidate_theme==t)==1},logical(1)))," region-unique units; ",nrow(z)," supported units"),
  historical_status=if(defs$route[i]=="HISTORICAL_HYPOTHESIS")collapse_unique(rr$current_status)else"NOT_HISTORICAL",
  evidence_strength_summary=paste0("minimum adjusted p=",fmt(min(p,na.rm=TRUE))),
  ambiguity_summary=paste0(sum(z$evidence_status=="REQUIRES_AUTHOR_REVIEW")," units require review"),
  recommended_action_options="RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER",
  author_decision="",author_notes="",stringsAsFactors=FALSE)
}))
wr(global,file.path(out,"candidate_theme_global_review.tsv"))

cond_section<-function(c,i){
 z<-packet[packet$condition==c,]
 c(paste0("## ",i,". ",c),"",paste0("- Decision units: ",nrow(z),". Priorities: HIGH ",sum(z$review_priority=="HIGH"),", MEDIUM ",sum(z$review_priority=="MEDIUM"),", LOW ",sum(z$review_priority=="LOW"),"."),paste0("- Regions: ",collapse_unique(z$FSF_region)),"")
}
broad<-global[global$cross_source_support_count>=8|grepl("Broad",global$candidate_theme),]
unsupported<-unique(recon$historical_theme[recon$current_status=="NOT_SUPPORTED_CURRENTLY"])
summary<-c("# Author biological theme review summary","","## 1. Overview","",
 sprintf("- %d decision units from 18 frozen candidate themes.",nrow(packet)),
 sprintf("- Review priorities: HIGH %d; MEDIUM %d; LOW %d.",sum(priority=="HIGH"),sum(priority=="MEDIUM"),sum(priority=="LOW")),
 "- Priority orders review only and does not select an action.","- Author decision and notes fields are blank.","",
 unlist(Map(cond_section,conditions,2:7)),
 "## 8. Cross-condition themes","",paste0("- ",global$candidate_theme,": ",global$conditions_supported),"",
 "## 9. Strongest region-specific biological differences","",paste0("- ",within$condition," ",within$region_A," vs ",within$region_B,": Jaccard ",sprintf("%.3f",within$Jaccard_similarity)),"",
 "## 10. Themes with broad/non-specific support","",paste0("- ",broad$candidate_theme,": ",broad$region_specificity_summary),"",
 "## 11. Historical themes not supported by current evidence","",paste0("- ",unsupported),"",
 "## 12. Ambiguous cases requiring author judgment","",paste0("- ",packet$decision_id[packet$review_priority=="LOW"],": ",packet$conflicting_or_ambiguous_evidence[packet$review_priority=="LOW"]),"",
 "## 13. Suggested review order","","- Review HIGH units first, followed by MEDIUM and LOW units.","- Suggested order is based only on deterministic evidence criteria.","",
 "Review document only; no action option has been selected.")
writeLines(summary,file.path(out,"author_biological_theme_review_summary.md"),useBytes=TRUE)
arts<-c("author_biological_theme_decision_packet.tsv","within_condition_biological_review.tsv","candidate_theme_global_review.tsv","author_biological_theme_review_summary.md")
wr(data.frame(artifact=arts,sha256=vapply(file.path(out,arts),sha,"")),file.path(out,"author_biological_theme_decision_hashes.tsv"))
cat("author biological theme decision packet build: PASS\n")
cat("decision units:",nrow(packet),"HIGH",sum(priority=="HIGH"),"MEDIUM",sum(priority=="MEDIUM"),"LOW",sum(priority=="LOW"),"\n")

