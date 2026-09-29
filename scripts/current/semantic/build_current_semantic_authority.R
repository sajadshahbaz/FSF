#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts/current/build_current_external_artifact_manifests.R"))
a <- commandArgs(trailingOnly = TRUE)
raw_suffix <- "current_fsf_v1/manuscript/semantic_authority/raw"
declared_root <- Sys.getenv("FSF_EXTERNAL_ARTIFACT_ROOT", unset = "")
if (length(a)) {
  raw_dir <- normalizePath(a[1], winslash = "/", mustWork = TRUE)
  external_root <- if (nzchar(declared_root)) {
    .fsf_external_root(declared_root)
  } else {
    .fsf_derive_external_root(raw_dir, raw_suffix)
  }
  raw_locator <- .fsf_external_locator_from_path(raw_dir, external_root)
} else {
  external_root <- .fsf_external_root(declared_root)
  raw_locator <- raw_suffix
  raw_dir <- .fsf_resolve_external_locator(raw_locator, external_root)
}
out_dir <- if (length(a) > 1) a[2] else file.path(root, "results/current_fsf_v1/manuscript/semantic_authority")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
date <- "2026-09-19"; gov <- "releases/2026-07-26"
kov <- "KEGG KO 2026-09-18"; kpv <- "KEGG PATHWAY 2026-09-17"
goh <- "0ded08c5f67abdcf6dbbc2d77bcd47530319ecdeeb7b6ad13cd402a65d5dcdc2"
keh <- "780415dbe1a73d9246b16b58b58892b24f20d89eb1d3cf9f6458479d2c2f856b"
gp <- file.path(root,"results/current_fsf_v1/manuscript/source_data/figure5_go_source.tsv")
kp <- file.path(root,"results/current_fsf_v1/manuscript/source_data/figure5_kegg_source.tsv")
sha <- function(p) sub(" .*","",system2("sha256sum",p,stdout=TRUE)[1])
rd <- function(p) read.delim(p,check.names=FALSE,quote="",stringsAsFactors=FALSE)
wr <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="")
stopifnot(sha(gp)==goh,sha(kp)==keh)
rf <- c(go="go-basic.obo",ko="kegg_list_ko.tsv",map="kegg_list_pathway.tsv",
 ko_path="kegg_list_pathway_ko.tsv",info="kegg_info_kegg.txt",
 info_ko="kegg_info_ko.txt",info_path="kegg_info_pathway.txt")
rp <- file.path(raw_dir,rf); names(rp)<-names(rf); stopifnot(all(file.exists(rp)))
rh <- c(go="b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc",
 ko="bf37ae22839f8805d3034a9e6cc76b55332f901aacaf10a5670749d43b5b0a92",
 map="46c36040530f1de249615192d7b054cefd366fac0d48243cf82bec22755cfe61",
 ko_path="f0905f26762b75cd2ca9337352f249a11adf939a3cd714c7b940954a428be03d",
 info="ce7f0636cb90da1fb1a92c86339fbd74a06f497be3143c709e64e767b61b8955",
 info_ko="9fdbd543cb2971901f0be21c87339a28db26c87bff9469333534c7a240f466ac",
 info_path="b59412b70632a0b1cae3e89d6b28bcb39919478f70b28117f05eda2b564421bd")
oh<-vapply(rp,sha,""); stopifnot(identical(unname(oh),unname(rh)))
go<-rd(gp); ke<-rd(kp)
go<-go[go$significant & go$gene_set_family=="main_class",]
ke<-ke[ke$significant & ke$gene_set_family=="main_class",]
gids<-sort(unique(go$go_term)); ki<-unique(ke[c("kegg_term","enrichment_type")])
ki<-ki[order(ki$enrichment_type,ki$kegg_term),]

parse_obo <- function(p) {
 x<-readLines(p,warn=FALSE); s<-which(x=="[Term]"); e<-c(s[-1]-1,length(x))
 one<-function(b,k,m=FALSE){z<-substring(b[startsWith(b,paste0(k,": "))],nchar(k)+3);if(!length(z))"" else if(m)paste(z,collapse=";") else z[1]}
 z<-lapply(seq_along(s),function(i){b<-x[(s[i]+1):e[i]];data.frame(
  canonical_go_id=one(b,"id"),go_name=one(b,"name"),go_namespace=one(b,"namespace"),
  obsolete_status=ifelse(tolower(one(b,"is_obsolete"))=="true","obsolete","current"),
  replaced_by=one(b,"replaced_by",TRUE),consider=one(b,"consider",TRUE),
  alt_id=one(b,"alt_id",TRUE),stringsAsFactors=FALSE)})
 do.call(rbind,z)
}
ont<-parse_obo(rp["go"]); stopifnot(any(ont$canonical_go_id=="GO:0008150"))
ci<-match(gids,ont$canonical_go_id); ai<-which(nzchar(ont$alt_id))
am<-do.call(rbind,lapply(ai,function(i)data.frame(go_id=strsplit(ont$alt_id[i],";",fixed=TRUE)[[1]],row=i)))
mi<-match(gids,am$go_id); oi<-ci; ua<-is.na(oi)&!is.na(mi); oi[ua]<-am$row[mi[ua]]; mg<-!is.na(oi)
ga<-data.frame(go_id=gids,canonical_go_id=ifelse(mg,ont$canonical_go_id[oi],""),
 go_name=ifelse(mg,ont$go_name[oi],""),go_namespace=ifelse(mg,ont$go_namespace[oi],""),
 obsolete_status=ifelse(mg,ont$obsolete_status[oi],"unmapped"),
 replacement_or_consider=ifelse(mg,ifelse(nzchar(ont$replaced_by[oi]),paste0("replaced_by:",ont$replaced_by[oi]),ifelse(nzchar(ont$consider[oi]),paste0("consider:",ont$consider[oi]),"")),""),
 source_version=gov,mapping_status=ifelse(!mg,"unmapped",ifelse(ua,"alternate_id",ifelse(ont$obsolete_status[oi]=="obsolete","obsolete","mapped"))))

kl<-function(p){x<-read.delim(p,header=FALSE,sep="\t",quote="",fill=TRUE,col.names=c("id","text"));x$id<-sub("^[^:]+:","",x$id);x}
kol<-kl(rp["ko"]); mpl<-kl(rp["map"]); kpl<-kl(rp["ko_path"])
ka<-do.call(rbind,lapply(seq_len(nrow(ki)),function(i){
 id<-ki$kegg_term[i]; ty<-ki$enrichment_type[i]; lu<-if(ty=="KO")kol else if(startsWith(id,"ko"))kpl else mpl
 j<-match(id,lu$id); ok<-!is.na(j); tx<-if(ok)lu$text[j] else ""
 if(ty=="KO"&&ok&&grepl("; ",tx,fixed=TRUE)){q<-strsplit(tx,"; ",fixed=TRUE)[[1]];nm<-q[1];ds<-paste(q[-1],collapse="; ")}else{nm<-tx;ds<-""}
 data.frame(kegg_id=id,kegg_identifier_type=ty,kegg_name=nm,kegg_description_if_available=ds,
 source_version=ifelse(ty=="KO",kov,kpv),mapping_status=ifelse(ok,"mapped","unmapped"))
}))

gj<-merge(go,ga,by.x="go_term",by.y="go_id",all.x=TRUE,sort=FALSE)
key<-function(x,c)do.call(paste,c(x[c],sep="\r"))
gj<-gj[match(key(go,c("gene_set_id","go_term")),key(gj,c("gene_set_id","go_term"))),];rownames(gj)<-NULL
kj<-merge(ke,ka,by.x=c("kegg_term","enrichment_type"),by.y=c("kegg_id","kegg_identifier_type"),all.x=TRUE,sort=FALSE)
kj<-kj[match(key(ke,c("gene_set_id","enrichment_type","kegg_term")),key(kj,c("gene_set_id","enrichment_type","kegg_term"))),];rownames(kj)<-NULL
stopifnot(nrow(gj)==nrow(go),nrow(kj)==nrow(ke),identical(gj$go_term,go$go_term),identical(kj$kegg_term,ke$kegg_term))

cov<-rbind(data.frame(source="GO",identifier_type="GO",total_unique_ids=nrow(ga),
 mapped_ids=sum(ga$mapping_status!="unmapped"),unmapped_ids=sum(ga$mapping_status=="unmapped"),
 obsolete_ids=sum(ga$obsolete_status=="obsolete"),ambiguous_or_replaced_ids=sum(ga$mapping_status=="alternate_id"|nzchar(ga$replacement_or_consider))),
 do.call(rbind,lapply(c("KO","PATHWAY"),function(t){z<-ka[ka$kegg_identifier_type==t,];data.frame(source="KEGG",identifier_type=t,total_unique_ids=nrow(z),mapped_ids=sum(z$mapping_status=="mapped"),unmapped_ids=sum(z$mapping_status=="unmapped"),obsolete_ids=NA_integer_,ambiguous_or_replaced_ids=NA_integer_)})))
ctx<-function(x,t){k<-paste(x$condition,x$stability_region,sep="\r");n<-tapply(x[[t]],k,function(z)length(unique(z)));as.integer(n[k])}
gpack<-data.frame(condition=gj$condition,FSF_region=gj$stability_region,source="GO",term_id=gj$go_term,
 official_term_name=gj$go_name,official_namespace_or_identifier_type=gj$go_namespace,
 adjusted_p=gj$adjusted_p_value,enrichment_measure=gj$enrichment_ratio,
 significant_term_count_context=ctx(gj,"go_term"),mapping_status=gj$mapping_status)
kpack<-data.frame(condition=kj$condition,FSF_region=kj$stability_region,source="KEGG",term_id=kj$kegg_term,
 official_term_name=kj$kegg_name,official_namespace_or_identifier_type=kj$enrichment_type,
 adjusted_p=kj$adjusted_p_value,enrichment_measure=kj$enrichment_ratio,
 significant_term_count_context=ctx(kj,"kegg_term"),mapping_status=kj$mapping_status)
packet<-rbind(gpack,kpack); co<-c("DES","GAM","HT","LT","OSM","UV");ro<-c("Low Stability","Transitional","Stable","Highly Stable")
packet<-packet[order(match(packet$condition,co),match(packet$FSF_region,ro),packet$source,packet$adjusted_p,packet$term_id),];rownames(packet)<-NULL

manifest<-data.frame(source=c("Gene Ontology",rep("KEGG",6)),
 resource_type=c("go-basic OBO","KO descriptions","PATHWAY map descriptions","PATHWAY ko descriptions","KEGG database metadata","KO metadata","PATHWAY metadata"),
 retrieval_date=date,resource_version_or_release=c(gov,kov,kpv,kpv,"KEGG 2026-09-18",kov,kpv),
 source_location=c("https://purl.obolibrary.org/obo/go/go-basic.obo","https://rest.kegg.jp/list/ko","https://rest.kegg.jp/list/pathway","https://rest.kegg.jp/list/pathway/ko","https://rest.kegg.jp/info/kegg","https://rest.kegg.jp/info/ko","https://rest.kegg.jp/info/pathway"),
 local_artifact=unname(vapply(rp,.fsf_external_locator_from_path,character(1L),root=external_root)),file_size=as.numeric(file.info(rp)$size),sha256=unname(oh),
 identifier_type=c("GO","KO","PATHWAY map","PATHWAY ko","metadata","metadata","metadata"),
 notes=c("Official GO basic ontology; version from OBO header","Official KEGG REST list; split at first semicolon","Official KEGG REST map pathway list","Official KEGG REST ko pathway list","Official KEGG REST release metadata","Official KEGG REST KO metadata","Official KEGG REST PATHWAY metadata"))
sets<-unique(rbind(go[c("condition","stability_region")],ke[c("condition","stability_region")]))
rich<-aggregate(term_id~condition+FSF_region+source,packet,function(z)length(unique(z)));names(rich)[4]<-"n";rich<-rich[order(-rich$n),]
sl<-c("# Biological semantic review summary","","Official names are attached to frozen current significant Figure 5 evidence. No theme or interpretation is assigned.","","## Coverage","",
 sprintf("- GO: %d unique; %d mapped; %d unmapped; %d obsolete; %d alternate/replaced.",cov$total_unique_ids[1],cov$mapped_ids[1],cov$unmapped_ids[1],cov$obsolete_ids[1],cov$ambiguous_or_replaced_ids[1]),
 sprintf("- KEGG KO: %d unique; %d mapped; %d unmapped.",cov$total_unique_ids[2],cov$mapped_ids[2],cov$unmapped_ids[2]),
 sprintf("- KEGG PATHWAY: %d unique; %d mapped; %d unmapped.",cov$total_unique_ids[3],cov$mapped_ids[3],cov$unmapped_ids[3]),
 sprintf("- Annotated evidence: %d GO rows and %d KEGG rows across %d condition × region sets.",nrow(gj),nrow(kj),nrow(sets)),"","## Evidence-rich condition/region combinations","",
 paste0("- ",head(rich$condition,12)," / ",head(rich$FSF_region,12)," / ",head(rich$source,12),": ",head(rich$n,12)," significant IDs"),"",
 "## Limitations","","- GO semantics are limited to the frozen go-basic snapshot.","- KEGG list responses provide KO names/descriptions and pathway names; pathway descriptions are unavailable in the list response.","- Semantic availability is not an interpretation or author decision.")
al<-c("# Semantic authority audit","",paste("- Build date:",date),paste("- Frozen Figure 5 GO SHA-256:",goh),paste("- Frozen Figure 5 KEGG SHA-256:",keh),
 "- Scope: significant TRUE and gene_set_family main_class.",sprintf("- GO: %d rows; %d IDs.",nrow(go),length(gids)),sprintf("- KEGG: %d rows; %d IDs.",nrow(ke),nrow(ki)),sprintf("- Evidence sets: %d.",nrow(sets)),
 "- Existing authoritative descriptions in frozen current evidence: none.","- Sources: official GO go-basic OBO and official KEGG REST lists.","- Original quantitative columns are preserved.","- Historical manual theme tables were not used.","- No themes or biological interpretations were created.","","## Coverage","",paste(capture.output(print(cov,row.names=FALSE)),collapse="\n"))
wr(manifest,file.path(out_dir,"semantic_source_manifest.tsv"));wr(ga,file.path(out_dir,"go_term_authority.tsv"));wr(ka,file.path(out_dir,"kegg_term_authority.tsv"));wr(cov,file.path(out_dir,"semantic_coverage.tsv"))
wr(ga[ga$mapping_status!="mapped",],file.path(out_dir,"go_unmapped_or_problematic.tsv"));wr(ka[ka$mapping_status!="mapped",],file.path(out_dir,"kegg_unmapped_or_problematic.tsv"))
wr(gj,file.path(out_dir,"current_significant_go_annotated.tsv"));wr(kj,file.path(out_dir,"current_significant_kegg_annotated.tsv"));wr(packet,file.path(out_dir,"biological_semantic_review_packet.tsv"))
writeLines(sl,file.path(out_dir,"biological_semantic_review_summary.md"),useBytes=TRUE);writeLines(al,file.path(out_dir,"semantic_authority_audit.md"),useBytes=TRUE)
arts<-c("semantic_source_manifest.tsv","go_term_authority.tsv","kegg_term_authority.tsv","semantic_coverage.tsv","go_unmapped_or_problematic.tsv","kegg_unmapped_or_problematic.tsv","current_significant_go_annotated.tsv","current_significant_kegg_annotated.tsv","biological_semantic_review_packet.tsv","biological_semantic_review_summary.md","semantic_authority_audit.md")
wr(data.frame(artifact=arts,sha256=vapply(file.path(out_dir,arts),sha,"")),file.path(out_dir,"semantic_output_hashes.tsv"))
cat("semantic authority build: PASS\n");print(cov,row.names=FALSE)

