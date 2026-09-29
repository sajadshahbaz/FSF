#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root<-normalizePath(getwd(),winslash="/",mustWork=TRUE)
source(file.path(root,"scripts/current/build_current_external_artifact_manifests.R"))
a<-commandArgs(trailingOnly=TRUE)
out<-if(length(a))a[1] else file.path(root,"results/current_fsf_v1/manuscript/biological_themes")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sem<-file.path(root,"results/current_fsf_v1/manuscript/semantic_authority")
review<-file.path(root,"results/current_fsf_v1/manuscript/review")
rd<-function(p)read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
wr<-function(x,p)write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="")
sha<-function(p)sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
regions<-c("Low Stability","Transitional","Stable","Highly Stable")
conditions<-c("DES","GAM","HT","LT","OSM","UV")
expected<-c(go_term_authority.tsv="608020c5f4c95263a34749b54b23a797c621184a27979dd59652506fcee1a218",
 kegg_term_authority.tsv="9652ad9bf7a89cf6cbb7f73da5cf99e8808241be284531bd47f85d21eef50671",
 current_significant_go_annotated.tsv="6bcf373e531b9dbfd8e1e89dc25ba77e2b2407c9c5dc3741408da281571786dd",
 current_significant_kegg_annotated.tsv="6bf9666965f84c7f183e75331f90cb878f2eabe0b8caf7adff6d520dfcb39735",
 biological_semantic_review_packet.tsv="c3b720079f405a8f9b20a4704435ae09e3a4e0887b559106bd4b85f996c42403")
stopifnot(all(vapply(file.path(sem,names(expected)),sha,"")==expected))
go<-rd(file.path(sem,"current_significant_go_annotated.tsv"))
ke<-rd(file.path(sem,"current_significant_kegg_annotated.tsv"))
ga<-rd(file.path(sem,"go_term_authority.tsv"));ka<-rd(file.path(sem,"kegg_term_authority.tsv"))
packet<-rd(file.path(sem,"biological_semantic_review_packet.tsv"))

defs<-data.frame(
 theme_id=c(sprintf("H%02d",1:13),sprintf("C%02d",1:5)),
 candidate_theme=c(
 "Translation, ribosome and protein targeting","DNA repair and genome maintenance",
 "Chromosome and nuclear organization","Cell cycle and division",
 "Metabolism and redox processes","Energy and core metabolism",
 "Protein turnover and proteostasis","Cell adhesion and extracellular organization",
 "Development and morphogenesis","Neural, projection and behavior processes",
 "Stress response and signaling","Endocrine and physiological regulation",
 "Broad regulatory and signaling processes","Membrane transport and ion homeostasis",
 "RNA processing and expression","Vesicle trafficking and organelle organization",
 "Immune and defense processes","Cytoskeleton and cell motility"),
 pattern=c(
 "ribosom|translation|translational|protein targeting|signal recognition particle|aminoacyl|trna synthet|peptide biosynth",
 "dna repair|dna recombin|dna damage|genome maintenance|double.strand break|nucleotide.excision|mismatch repair|homologous recombination|base.excision",
 "chromosom|chromatin|nucleosome|nuclear organization|nucleus organization|kinetochore|centromer|telomer|histone",
 "cell cycle|cell division|mitotic|mitosis|meiotic|meiosis|cytokinesis|spindle|checkpoint",
 "metabol|oxidoreduct|redox|oxidation.reduction|biosynth|catabol|dehydrogenase|transferase activity",
 "energy|oxidative phosphorylation|electron transport|respiratory|atp synth|carbon metabolism|glycolysis|citrate cycle|tca cycle",
 "proteasom|ubiquitin|protein degrad|protein catabol|proteolysis|peptidase|proteostasis|protein folding|chaperon|autophag",
 "cell adhesion|extracellular|cell junction|focal adhesion|matrix|cadherin|integrin",
 "develop|morphogen|differentiation|organ formation|anatomical structure|growth",
 "neural|neuron|axon|dendrit|synap|behavior|neurotrans|projection",
 "stress|response to stimulus|signaling pathway|signal transduction|mapk|heat shock|oxidative stress",
 "endocrine|hormone|physiological|insulin|thyroid|steroid|secretion",
 "regulation|regulatory|signaling|signal transduction|kinase activity|transcription factor|protein binding",
 "transport|transmembrane|channel activity|ion homeostasis|ion transport|porter activity|pump activity",
 "rna process|rna splicing|rna modification|rna metabol|transcription|mrna|rrna|ncrna|gene expression",
 "vesicle|endocyt|exocyt|golgi|endoplasmic reticulum|organelle organization|intracellular transport|protein localization",
 "immune|defense response|inflammatory|cytokine|antigen|infection|host defense|innate immunity",
 "cytoskeleton|actin|microtubule|cell motility|cell migration|cilium|flagell"),
 route=c(rep("HISTORICAL_HYPOTHESIS",13),rep("CURRENT_EVIDENCE",5)),
 stringsAsFactors=FALSE)
wr(defs,file.path(out,"candidate_theme_definitions.tsv"))

map_terms<-function(x,source,id,name,extra=""){
 txt<-tolower(paste(x[[name]],if(nzchar(extra))x[[extra]] else ""))
 z<-lapply(seq_len(nrow(defs)),function(i){
  hit<-grepl(defs$pattern[i],txt,perl=TRUE)
  if(!any(hit))return(NULL)
  d<-x[hit,]
  basis<-if(defs$route[i]=="HISTORICAL_HYPOTHESIS")"HISTORICAL_HYPOTHESIS_SUPPORTED" else
   if(source=="KEGG"&&"enrichment_type"%in%names(d)&&any(d$enrichment_type=="PATHWAY"))"KEGG_PATHWAY_IDENTITY" else "CURRENT_EVIDENCE_DISCOVERY"
  data.frame(theme_id=defs$theme_id[i],candidate_theme=defs$candidate_theme[i],source=source,
   term_id=d[[id]],official_term_name=d[[name]],condition=d$condition,FSF_region=d$stability_region,
   adjusted_p=d$adjusted_p_value,enrichment_measure=d$enrichment_ratio,mapping_basis=basis,
   historical_theme_relation=ifelse(defs$route[i]=="HISTORICAL_HYPOTHESIS","HISTORICAL_HYPOTHESIS","NEW_CURRENT_CANDIDATE"),
   evidence_status=if(source=="GO")ifelse(d$obsolete_status=="obsolete","REQUIRES_AUTHOR_REVIEW","CURRENT_TERM_SUPPORT") else "CURRENT_TERM_SUPPORT",
   author_decision="",author_notes="",stringsAsFactors=FALSE)
 })
 do.call(rbind,z)
}
gm<-map_terms(go,"GO","go_term","go_name")
km<-map_terms(ke,"KEGG","kegg_term","kegg_name","kegg_description_if_available")
mapping<-rbind(gm,km)
mapping<-unique(mapping)
mapping<-mapping[order(match(mapping$condition,conditions),match(mapping$FSF_region,regions),mapping$theme_id,mapping$source,mapping$adjusted_p,mapping$term_id),]
rownames(mapping)<-NULL
wr(mapping,file.path(out,"candidate_term_theme_mapping.tsv"))
trace_key<-function(condition,region,source,term)paste(condition,region,source,term,sep="\r")
mk<-trace_key(mapping$condition,mapping$FSF_region,mapping$source,mapping$term_id)
theme_ids_by_term<-tapply(mapping$theme_id,mk,function(z)paste(sort(unique(z)),collapse=";"))
theme_names_by_term<-tapply(mapping$candidate_theme,mk,function(z)paste(sort(unique(z)),collapse="; "))
pk<-trace_key(packet$condition,packet$FSF_region,packet$source,packet$term_id)
trace<-packet
trace$evidence_row_id<-sprintf("SEMANTIC_EVIDENCE_%05d",seq_len(nrow(trace)))
trace$candidate_theme_ids<-unname(theme_ids_by_term[pk])
trace$candidate_themes<-unname(theme_names_by_term[pk])
trace$candidate_theme_ids[is.na(trace$candidate_theme_ids)]<-""
trace$candidate_themes[is.na(trace$candidate_themes)]<-""
trace$theme_assignment_status<-ifelse(nzchar(trace$candidate_theme_ids),
 "CANDIDATE_MAPPING_AVAILABLE","UNASSIGNED_EVIDENCE_RETAINED")
trace<-trace[c("evidence_row_id",names(packet),"candidate_theme_ids",
 "candidate_themes","theme_assignment_status")]
wr(trace,file.path(out,"biological_theme_evidence_trace.tsv"))

manifest<-rd(file.path(sem,"semantic_source_manifest.tsv"))
obo_locator<-manifest$local_artifact[manifest$resource_type=="go-basic OBO"]
stopifnot(length(obo_locator)==1)
external_root<-.fsf_external_root()
obo<-.fsf_resolve_external_locator(obo_locator,external_root)
stopifnot(sha(obo)==manifest$sha256[manifest$resource_type=="go-basic OBO"])
parse_parents<-function(p){
 x<-readLines(p,warn=FALSE);s<-which(x=="[Term]");e<-c(s[-1]-1,length(x));ans<-list()
 for(i in seq_along(s)){b<-x[(s[i]+1):e[i]];id<-sub("^id: ","",b[startsWith(b,"id: ")][1]);pa<-sub(" !.*$","",sub("^is_a: ","",b[startsWith(b,"is_a: ")]));if(length(id)&&!is.na(id))ans[[id]]<-pa}
 ans
}
parents<-parse_parents(obo)
make_groups<-function(d){
 ids<-sort(unique(d$go_term));can<-setNames(d$canonical_go_id[match(ids,d$go_term)],ids)
 n<-length(ids);par<-seq_len(n)
 find<-function(i){while(par[i]!=i){par[i]<<-par[par[i]];i<-par[i]};i}
 unite<-function(i,j){a<-find(i);b<-find(j);if(a!=b)par[b]<<-a}
 cids<-unname(can);names(cids)<-ids;ix<-setNames(seq_len(n),cids)
 for(i in seq_len(n)){ps<-parents[[cids[i]]];for(p in ps)if(p%in%names(ix))unite(i,ix[[p]])}
 inv<-list()
 for(i in seq_len(n))for(p in parents[[cids[i]]])inv[[p]]<-c(inv[[p]],i)
 for(v in inv)if(length(v)>1)for(j in v[-1])unite(v[1],j)
 comp<-split(seq_len(n),vapply(seq_len(n),find,integer(1)))
 do.call(rbind,lapply(seq_along(comp),function(k){
  mem<-ids[comp[[k]]];q<-d[d$go_term%in%mem,];best<-q[which.min(q$adjusted_p_value),]
  data.frame(redundancy_group_id=sprintf("GR_%s_%s_%04d",unique(d$condition),gsub(" ","_",unique(d$stability_region)),k),
   representative_term_id=best$go_term,representative_term_name=best$go_name,
   member_term_ids=paste(mem,collapse=";"),member_count=length(mem),condition=unique(d$condition),
   FSF_region=unique(d$stability_region),selection_basis="DIRECT_IS_A_OR_SHARED_PARENT; representative has minimum adjusted p",stringsAsFactors=FALSE)
 }))
}
gsplit<-split(go,paste(go$condition,go$stability_region,sep="\r"))
groups<-do.call(rbind,lapply(gsplit,make_groups));rownames(groups)<-NULL
wr(groups,file.path(out,"go_redundancy_groups.tsv"))

sets<-unique(rbind(go[c("condition","stability_region")],ke[c("condition","stability_region")]))
sets<-sets[order(match(sets$condition,conditions),match(sets$stability_region,regions)),]
collapse_unique<-function(x)paste(sort(unique(x[nzchar(x)])),collapse="; ")
inventory<-do.call(rbind,lapply(seq_len(nrow(sets)),function(i){
 c<-sets$condition[i];r<-sets$stability_region[i];g<-go[go$condition==c&go$stability_region==r,];k<-ke[ke$condition==c&ke$stability_region==r,];gr<-groups[groups$condition==c&groups$FSF_region==r,]
 data.frame(condition=c,FSF_region=r,significant_GO_terms=length(unique(g$go_term)),significant_KEGG_terms=length(unique(k$kegg_term)),
  strongest_GO_adjusted_p=if(nrow(g))min(g$adjusted_p_value)else NA,strongest_KEGG_adjusted_p=if(nrow(k))min(k$adjusted_p_value)else NA,
  GO_namespaces=collapse_unique(g$go_namespace),official_GO_names=collapse_unique(g$go_name),official_KEGG_names=collapse_unique(k$kegg_name),
  GO_redundancy_groups=nrow(gr),GO_terms_in_multiterm_groups=sum(gr$member_count[gr$member_count>1]),
  obsolete_GO_terms=length(unique(g$go_term[g$obsolete_status=="obsolete"])),
  unresolved_KEGG_identifiers=collapse_unique(k$kegg_term[k$mapping_status=="unmapped"]),stringsAsFactors=FALSE)
}))
wr(inventory,file.path(out,"current_biological_evidence_inventory.tsv"))

all_evidence<-expand.grid(condition=sets$condition[match(unique(paste(sets$condition,sets$stability_region)),paste(sets$condition,sets$stability_region))],
 FSF_region=sets$stability_region[match(unique(paste(sets$condition,sets$stability_region)),paste(sets$condition,sets$stability_region))],
 theme_id=defs$theme_id,stringsAsFactors=FALSE)
all_evidence<-merge(sets,defs,by=NULL)
summ<-lapply(seq_len(nrow(all_evidence)),function(i){
 z<-all_evidence[i,];m<-mapping[mapping$condition==z$condition&mapping$FSF_region==z$stability_region&mapping$theme_id==z$theme_id,]
 g<-m[m$source=="GO",];k<-m[m$source=="KEGG",]
 gr<-groups[groups$condition==z$condition&groups$FSF_region==z$stability_region,]
 gc<-if(nrow(g))length(unique(g$term_id))else 0;kc<-if(nrow(k))length(unique(k$term_id))else 0
 rg<-if(gc)sum(vapply(strsplit(gr$member_term_ids,";",fixed=TRUE),function(q)any(q%in%g$term_id),logical(1)))else 0
 es<-if(gc>0&&kc>0)"STRONG_MULTI_SOURCE_SUPPORT" else if(gc>0)"GO_SUPPORTED" else if(kc>0)"KEGG_SUPPORTED" else "NO_CURRENT_SUPPORT"
 if(any(m$evidence_status=="REQUIRES_AUTHOR_REVIEW"))es<-"REQUIRES_AUTHOR_REVIEW"
 data.frame(condition=z$condition,FSF_region=z$stability_region,candidate_theme=z$candidate_theme,
  GO_term_count_raw=gc,GO_redundancy_group_count=rg,KEGG_term_count=kc,
  strongest_GO_adjusted_p=if(nrow(g))min(g$adjusted_p)else NA,strongest_KEGG_adjusted_p=if(nrow(k))min(k$adjusted_p)else NA,
  GO_support=collapse_unique(paste(g$term_id,g$official_term_name,sep=": ")),
  KEGG_support=collapse_unique(paste(k$term_id,k$official_term_name,sep=": ")),
  cross_source_support=ifelse(gc>0&&kc>0,"YES","NO"),evidence_status=es,
  historical_relation=ifelse(z$route=="HISTORICAL_HYPOTHESIS","HISTORICAL_CANDIDATE","NEW_CURRENT_CANDIDATE"),
  author_decision="",author_notes="",stringsAsFactors=FALSE)
})
evidence<-do.call(rbind,summ)
evidence<-evidence[order(match(evidence$condition,conditions),match(evidence$FSF_region,regions),evidence$candidate_theme),];rownames(evidence)<-NULL
wr(evidence,file.path(out,"candidate_theme_evidence.tsv"))

supported<-evidence[evidence$evidence_status!="NO_CURRENT_SUPPORT",]
profile<-do.call(rbind,lapply(seq_len(nrow(sets)),function(i){
 c<-sets$condition[i];r<-sets$stability_region[i];z<-supported[supported$condition==c&supported$FSF_region==r,]
 score<-pmin(z$strongest_GO_adjusted_p,z$strongest_KEGG_adjusted_p,na.rm=TRUE);score[!is.finite(score)]<-ifelse(is.na(z$strongest_GO_adjusted_p),z$strongest_KEGG_adjusted_p,z$strongest_GO_adjusted_p)
 ord<-order(score)
 data.frame(condition=c,FSF_region=r,candidate_themes=collapse_unique(z$candidate_theme),supported_theme_count=nrow(z),
  strongest_themes=paste(head(z$candidate_theme[ord],5),collapse="; "),
  GO_only_themes=collapse_unique(z$candidate_theme[z$GO_term_count_raw>0&z$KEGG_term_count==0]),
  KEGG_only_themes=collapse_unique(z$candidate_theme[z$GO_term_count_raw==0&z$KEGG_term_count>0]),
  cross_source_themes=collapse_unique(z$candidate_theme[z$GO_term_count_raw>0&z$KEGG_term_count>0]),
  unresolved_or_ambiguous_evidence=collapse_unique(c(inventory$unresolved_KEGG_identifiers[inventory$condition==c&inventory$FSF_region==r],
   z$candidate_theme[z$evidence_status=="REQUIRES_AUTHOR_REVIEW"])),stringsAsFactors=FALSE)
}))
wr(profile,file.path(out,"condition_region_biological_profile.tsv"))

jac<-function(a,b){u<-union(a,b);if(!length(u))NA_real_ else length(intersect(a,b))/length(u)}
multi<-c("DES","GAM","LT","UV")
contr<-do.call(rbind,lapply(multi,function(c){
 rs<-sets$stability_region[sets$condition==c]
 if(length(rs)<2)return(NULL)
 cmb<-combn(rs,2,simplify=FALSE)
 do.call(rbind,lapply(cmb,function(v){
  a<-supported$candidate_theme[supported$condition==c&supported$FSF_region==v[1]]
  b<-supported$candidate_theme[supported$condition==c&supported$FSF_region==v[2]]
  ga<-unique(go$go_term[go$condition==c&go$stability_region==v[1]]);gb<-unique(go$go_term[go$condition==c&go$stability_region==v[2]])
  kaa<-unique(ke$kegg_term[ke$condition==c&ke$stability_region==v[1]]);kb<-unique(ke$kegg_term[ke$condition==c&ke$stability_region==v[2]])
  data.frame(condition=c,region_A=v[1],region_B=v[2],themes_shared=collapse_unique(intersect(a,b)),
   themes_unique_A=collapse_unique(setdiff(a,b)),themes_unique_B=collapse_unique(setdiff(b,a)),
   GO_overlap=jac(ga,gb),KEGG_overlap=jac(kaa,kb),profile_similarity_measure=jac(a,b),
   interpretation_status="DESCRIPTIVE_COMPARISON_REQUIRES_AUTHOR_REVIEW",author_decision="",author_notes="",stringsAsFactors=FALSE)
 }))
}))
wr(contr,file.path(out,"within_condition_region_contrast.tsv"))

hist<-rd(file.path(review,"combined_profile_review.tsv"))
stopifnot(nrow(hist)==40L)
historical_locator<-"current_fsf_v1/manuscript/historical_authority/Table_SignalClass_Biology_CURATED.tsv"
historical_path<-.fsf_resolve_external_locator(historical_locator,external_root)
historical_source_sha256<-"c4704fb6d59e0185919ea75ab9b3ee42a004f64da39a264ac0f5c169febf2529"
stopifnot(sha(historical_path)==historical_source_sha256)
historical<-rd(historical_path)
historical$key<-paste(toupper(historical$condition),
 ifelse(historical$stability_level=="instable","Low Stability",
  ifelse(historical$stability_level=="weakly_stable_transitional","Transitional",
   ifelse(historical$stability_level=="highly_stable","Highly Stable","Stable"))),
 ifelse(historical$dominant_state=="mixed","tied",historical$dominant_state),sep="\r")
hist$key<-paste(hist$condition,hist$stability_region,hist$dominant_state,sep="\r")
historical_match_count<-vapply(hist$key,function(k)sum(historical$key==k),integer(1))
stopifnot(length(historical_match_count)==40L,
 sum(historical_match_count==0L)==0L,
 sum(historical_match_count>1L)==0L)
historical_index<-match(hist$key,historical$key)
hist$historical_profile<-historical$final_biological_theme[historical_index]
hist$historical_GO_theme<-historical$GO_theme[historical_index]
hist$historical_KEGG_theme<-historical$KEGG_theme[historical_index]
hist_rows<-unique(rbind(
 data.frame(historical_theme=hist$historical_profile,historical_condition=hist$condition,historical_class_or_region=hist$region_or_class),
 data.frame(historical_theme=hist$historical_GO_theme,historical_condition=hist$condition,historical_class_or_region=hist$region_or_class),
 data.frame(historical_theme=hist$historical_KEGG_theme,historical_condition=hist$condition,historical_class_or_region=hist$region_or_class)))
hist_rows<-hist_rows[!is.na(hist_rows$historical_theme)&nzchar(hist_rows$historical_theme),]
theme_alias<-function(x){
 x<-tolower(x)
 hit<-vapply(seq_len(nrow(defs)),function(i)grepl(defs$pattern[i],x,perl=TRUE),logical(1))
 defs$candidate_theme[hit]
}
recon<-do.call(rbind,lapply(seq_len(nrow(hist_rows)),function(i){
 h<-hist_rows[i,];r<-if(h$historical_class_or_region%in%regions)h$historical_class_or_region else ""
 cand<-theme_alias(h$historical_theme)
 z<-supported[supported$condition==h$historical_condition&(r==""|supported$FSF_region==r)&supported$candidate_theme%in%cand,]
 gs<-sum(z$GO_term_count_raw);ks<-sum(z$KEGG_term_count)
 status<-if(!length(cand))"AMBIGUOUS_REQUIRES_AUTHOR_REVIEW" else if(nrow(z)&&gs>0&&ks>0)"SUPPORTED_CURRENTLY" else if(nrow(z))"PARTIALLY_SUPPORTED" else "NOT_SUPPORTED_CURRENTLY"
 data.frame(historical_theme=h$historical_theme,historical_condition=h$historical_condition,
  historical_class_or_region=h$historical_class_or_region,current_matching_evidence=collapse_unique(z$candidate_theme),
  current_GO_support=gs,current_KEGG_support=ks,current_status=status,
  reason=ifelse(!length(cand),"No transparent match to a current candidate definition",ifelse(nrow(z),"Current official-name evidence matched","No current official-name evidence matched")),
  author_decision="",author_notes="",stringsAsFactors=FALSE)
}))
wr(recon,file.path(out,"historical_theme_reconciliation.tsv"))

review_units<-supported
review_packet<-data.frame(decision_id=sprintf("BIO_THEME_%04d",seq_len(nrow(review_units))),
 condition=review_units$condition,FSF_region=review_units$FSF_region,candidate_theme=review_units$candidate_theme,
 GO_support_summary=review_units$GO_support,KEGG_support_summary=review_units$KEGG_support,
 strongest_evidence=pmin(review_units$strongest_GO_adjusted_p,review_units$strongest_KEGG_adjusted_p,na.rm=TRUE),
 historical_relation=review_units$historical_relation,evidence_status=review_units$evidence_status,
 recommended_action_options="RETAIN;RETAIN_WITH_REVISED_WORDING;MERGE_WITH_OTHER_THEME;DROP;REVIEW_FURTHER",
 author_decision="",author_notes="",stringsAsFactors=FALSE)
review_packet$strongest_evidence[!is.finite(review_packet$strongest_evidence)]<-ifelse(is.na(review_units$strongest_GO_adjusted_p),review_units$strongest_KEGG_adjusted_p,review_units$strongest_GO_adjusted_p)[!is.finite(review_packet$strongest_evidence)]
wr(review_packet,file.path(out,"biological_theme_author_review_packet.tsv"))

recurring<-aggregate(condition~candidate_theme,supported,function(x)length(unique(x)));names(recurring)[2]<-"conditions_with_support";recurring<-recurring[order(-recurring$conditions_with_support,recurring$candidate_theme),]
sec<-function(c,i){z<-profile[profile$condition==c,];c(paste0("## ",i,". ",c),"",if(!nrow(z))"- No evidence-bearing region." else paste0("- ",z$FSF_region,": ",z$supported_theme_count," candidate themes; strongest: ",z$strongest_themes),"")}
summary<-c("# Biological theme review summary","","## 1. Overview","",
 sprintf("- %d candidate definitions; %d evidence-bearing condition × region sets; %d review decision units.",nrow(defs),nrow(sets),nrow(review_packet)),
 sprintf("- %d semantic evidence rows are retained; %d have no candidate mapping.",nrow(trace),sum(trace$theme_assignment_status=="UNASSIGNED_EVIDENCE_RETAINED")),
 "- Evidence is based only on official names and frozen enrichment statistics.","- Author decisions remain blank.","",
 unlist(Map(sec,conditions,2:7)),
 "## 8. Cross-condition recurring themes","",paste0("- ",head(recurring$candidate_theme,15),": support in ",head(recurring$conditions_with_support,15)," conditions"),"",
 "## 9. Region-specific functional differences","",paste0("- ",contr$condition," ",contr$region_A," vs ",contr$region_B,": theme Jaccard ",format(round(contr$profile_similarity_measure,3),nsmall=3)),"",
 "## 10. Historical-theme reconciliation","",paste0("- ",names(table(recon$current_status)),": ",as.integer(table(recon$current_status))),"",
 "## 11. Ambiguous evidence","",paste0("- GO terms flagged for review in mappings: ",length(unique(mapping$term_id[mapping$evidence_status=="REQUIRES_AUTHOR_REVIEW"]))),paste0("- Unresolved KEGG identifiers in evidence inventory: ",sum(nzchar(inventory$unresolved_KEGG_identifiers))),"",
 "## 12. Unsupported historical themes","",paste0("- ",unique(recon$historical_theme[recon$current_status=="NOT_SUPPORTED_CURRENTLY"])),"",
 "## 13. Candidate themes suitable for eventual visualization","",paste0("- ",defs$candidate_theme),"",
 "This document reports review evidence only and contains no selected themes or final biological conclusions.")
writeLines(summary,file.path(out,"biological_theme_review_summary.md"),useBytes=TRUE)

arts<-c("candidate_theme_definitions.tsv","current_biological_evidence_inventory.tsv","biological_theme_evidence_trace.tsv","candidate_term_theme_mapping.tsv","go_redundancy_groups.tsv","candidate_theme_evidence.tsv","condition_region_biological_profile.tsv","within_condition_region_contrast.tsv","historical_theme_reconciliation.tsv","biological_theme_author_review_packet.tsv","biological_theme_review_summary.md")
hashes<-data.frame(artifact=arts,sha256=vapply(file.path(out,arts),sha,""))
wr(hashes,file.path(out,"biological_theme_output_hashes.tsv"))
cat("biological theme review build: PASS\n")
cat("candidate themes:",nrow(defs),"mapping rows:",nrow(mapping),"review units:",nrow(review_packet),"\n")

