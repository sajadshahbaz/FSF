#!/usr/bin/env python3
import csv,gzip,hashlib,re,shutil,sys
from datetime import datetime
from pathlib import Path
from openpyxl import Workbook
from openpyxl.styles import Alignment,Font,PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import Table,TableStyleInfo
csv.field_size_limit(sys.maxsize)
ROOT,DEV,OUT=map(lambda x:Path(x).resolve(),sys.argv[1:4])
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rr(p,d='\t'):
 with open(p,newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter=d))
def wr(p,fields,data,d='\t'):
 with open(p,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter=d,lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(data)
def cp(src,name):shutil.copyfile(src,OUT/name)
def gzbytes(src,dst):
 with open(src,'rb') as i,open(dst,'wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=6,mtime=0) as o:shutil.copyfileobj(i,o)
def gzrows(dst,fields,data):
 payload=('\t'.join(fields)+'\n'+''.join('\t'.join(r[x] for x in fields)+'\n' for r in data)).encode()
 with open(dst,'wb') as raw:
  with gzip.GzipFile(filename='',mode='wb',fileobj=raw,compresslevel=6,mtime=0) as o:o.write(payload)
if OUT.exists():raise SystemExit(f'OUTPUT_ALREADY_EXISTS={OUT}')
OUT.mkdir(parents=False)
regions=rr(ROOT/'results/current_fsf_v1/current_fsf_region_architecture.tsv');metrics=rr(ROOT/'results/current_fsf_v1/current_fsf_feature_metrics.tsv')
conds=['DES','GAM','HT','LT','OSM','UV'];regs=['Low Stability','Transitional','Stable','Highly Stable'];tf=['condition','n_features','n_perturbations','low_stability_n','low_stability_percent','transitional_n','transitional_percent','stable_n','stable_percent','highly_stable_n','highly_stable_percent','design_note'];note='Single perturbation contrast; SSI = 1 is structurally determined and does not demonstrate cross-perturbation reproducibility.';t=[]
for c in conds:
 m=[x for x in metrics if x['condition']==c];ids={x['feature_id'] for x in m};np={x['n_perturbations'] for x in m};r={x['stability_region']:x for x in regions if x['condition']==c}
 if len(m)!=len(ids) or len(np)!=1 or set(r)!=set(regs):raise SystemExit(f'TABLE1_INPUT_INCONSISTENT={c}')
 v=[]
 for q in regs:n=int(r[q]['region_count']);v.extend([str(n),f'{100*n/len(ids):.1f}'])
 t.append(dict(zip(tf,[c,str(len(ids)),next(iter(np)),*v,note if c in ('HT','OSM') else ''])))
wr(OUT/'Table1_Condition_Level_FSF_Stability_Region_Architecture.tsv',tf,t)
s1=DEV/'governed_historical_authority/Supplementary_Table_S1_metadata.csv'
if sha(s1)!='6da9e049d318a88a37a80d10605580971a2ee6c8e45350518e990ceeb362bbf2':raise SystemExit('S1_HASH')
cp(s1,'Supplementary_Table_S1_metadata.csv');ann=ROOT/'governed_artifacts/FSF_v1_annotation_master.tsv'
if sha(ann)!='25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584':raise SystemExit('S2_HASH')
gzbytes(ann,OUT/'Supplementary_Table_S2_Functional_Annotation_Atlas.tsv.gz')
bio=ROOT/'results/current_fsf_v1/manuscript/biological_validation';short={(x['condition'],x['feature_id']):x for x in rr(bio/'validated_gene_shortlist.tsv')};dec={(x['condition'],x['feature_id']):x for x in rr(bio/'gene_author_decisions.tsv')};mi={(x['condition'],x['feature_id']):x for x in metrics};sf=['condition','feature_id','preferred_gene_name','annotation_description','signal_class','stability_region','dominant_state','n_perturbations','p_up','p_down','p_const','ssi','candidate_rank','author_decision','pmid','doi','literature_citation','evidence_role','supports_condition_relevance','supports_directionality','interpretation_note'];s3=[]
for e in rr(bio/'gene_literature_evidence.tsv'):
 k=(e['condition'],e['feature_id'])
 if k not in short:continue
 s,d,m=short[k],dec[k],mi[k];s3.append(dict(zip(sf,[k[0],k[1],s['eggnog_preferred_name'],s['eggnog_description'],m['signal_class'],m['stability_region'],m['dominant_state'],m['n_perturbations'],m['p_up'],m['p_down'],m['p_const'],m['ssi'],s['candidate_rank'],d['author_decision'],e['pmid'],e['doi'],e['citation'],e['evidence_role'],e['supports_condition_relevance'],e['supports_directionality'],d['author_notes']])))
if len(s3)!=8:raise SystemExit('S3_ROWS')
wr(OUT/'Supplementary_Table_S3_Literature_Supported_Representative_Genes.tsv',sf,s3)
def enrich(src,dst,n):
 with open(src,newline='',encoding='utf-8') as f:r=csv.DictReader(f,delimiter='\t');fields=r.fieldnames;data=[x for x in r if x['gene_set_family']=='main_class' and x['significant']=='TRUE']
 if len(data)!=n:raise SystemExit(f'{dst}_ROWS')
 gzrows(OUT/dst,fields,data)
td=ROOT/'results/current_fsf_v1/manuscript/tables';enrich(td/'current_go_statistical_summary.tsv','Supplementary_Table_S4_GO_Enrichment_Main_Class.tsv.gz',6011);enrich(td/'current_kegg_statistical_summary.tsv','Supplementary_Table_S5_KEGG_Enrichment_Main_Class.tsv.gz',986)
bd=ROOT/'results/current_fsf_v1/manuscript/benchmarks';bs=rr(bd/'baseline_summary.tsv');bc=rr(bd/'baseline_truth_by_class.tsv');s6f=['truth_scenario','n_features','mean_ssi','median_ssi','mean_stability_deviation','median_stability_deviation','predicted_signal_class','class_count','class_proportion'];s6=[]
for b in bs:
 for c in sorted((x for x in bc if x['truth_scenario']==b['truth_scenario']),key=lambda x:x['signal_class']):s6.append({**b,'predicted_signal_class':c['signal_class'],'class_count':c['n_features'],'class_proportion':c['proportion']})
wr(OUT/'Supplementary_Table_S6_Synthetic_Baseline_Benchmark.tsv',s6f,s6)
direct={'noise_summary.tsv':'Supplementary_Table_S7A_Noise_Gradient_SSI_Summary.tsv','noise_dominant_recovery.tsv':'Supplementary_Table_S7B_Noise_Gradient_Dominant_Class_Recovery.tsv','noise_truth_by_class.tsv':'Supplementary_Table_S7C_Noise_Gradient_Class_Assignments.tsv','tau_truth_by_class.tsv':'Supplementary_Table_S8A_Tau_Sensitivity_Class_Composition.tsv','tau_summary.tsv':'Supplementary_Table_S8B_Tau_Sensitivity_SSI_Summary.tsv'}
for a,b in direct.items():cp(bd/a,b)
meta=[('Table 1','Table1_Condition_Level_FSF_Stability_Region_Architecture.tsv','Table 1. Condition-level FSF stability-region architecture across six environmental stress conditions','results/current_fsf_v1/current_fsf_region_architecture.tsv; results/current_fsf_v1/current_fsf_feature_metrics.tsv; results/current_fsf_v1/manuscript/tables/current_region_architecture.tsv','Results 3.1','Figure 3','Exact region counts; percentages displayed to one decimal place; single-contrast HT/OSM guardrail retained.'),('S1','Supplementary_Table_S1_metadata.csv','Supplementary Table S1. RNA-seq sample metadata used for FSF analysis','Git object d13ab123b021672d20485f14b1381214175ac0e9:results/manuscript/tables/supplementary/locked_submission_tables/Supplementary_Table_S1_metadata.csv','Methods 2.1','','Preserved byte-identical shared metadata table.'),('S2','Supplementary_Table_S2_Functional_Annotation_Atlas.tsv.gz','Supplementary Table S2. Functional annotation atlas of the FSF feature universe','governed_artifacts/FSF_v1_annotation_master.tsv (SHA-256 25627cfcaf544dd2793d9fd2462772d9cb76f4728bbcab8c9d29599520f5d584)','Methods 2.6 and annotation-coverage Results subsection','Supplementary Figure S1','One frozen annotation-master row per analyzed feature; no external enrichment or inferred names.')]
# Remaining descriptive rows share the verified schema and do not encode scientific values.
extra=[('S3','Supplementary_Table_S3_Literature_Supported_Representative_Genes.tsv','Literature-supported representative genes','validated_gene_shortlist.tsv; gene_literature_evidence.tsv; gene_author_decisions.tsv; current_fsf_feature_metrics.tsv','Results 3.4','','Eight author-retained literature-supported representative genes; directionality support is not claimed.'),('S4','Supplementary_Table_S4_GO_Enrichment_Main_Class.tsv.gz','Significant GO enrichment results','results/current_fsf_v1/manuscript/tables/current_go_statistical_summary.tsv','Results 3.3','Figure 5','Raw unique GO term IDs; main_class significant rows only; adjusted P <= 0.05.'),('S5','Supplementary_Table_S5_KEGG_Enrichment_Main_Class.tsv.gz','Significant KEGG enrichment results','results/current_fsf_v1/manuscript/tables/current_kegg_statistical_summary.tsv','Results 3.3','Figure 5','Raw unique KEGG IDs retained; ko/map-normalized pathway counts are not substituted.'),('S6','Supplementary_Table_S6_Synthetic_Baseline_Benchmark.tsv','Baseline synthetic benchmark','baseline_summary.tsv; baseline_truth_by_class.tsv','Synthetic-baseline Results subsection','Figure 6','Generator truth-scenario labels are provenance labels.'),('S7A','Supplementary_Table_S7A_Noise_Gradient_SSI_Summary.tsv','Noise-gradient SSI summary','noise_summary.tsv','Noise-robustness Results subsection','','Authority preserved exactly.'),('S7B','Supplementary_Table_S7B_Noise_Gradient_Dominant_Class_Recovery.tsv','Noise-gradient dominant-class recovery','noise_dominant_recovery.tsv','Noise-robustness Results subsection','Figure 6','Direct Figure 6 authority.'),('S7C','Supplementary_Table_S7C_Noise_Gradient_Class_Assignments.tsv','Noise-gradient class assignments','noise_truth_by_class.tsv','Noise-robustness Results subsection','','Authority preserved exactly.'),('S8A','Supplementary_Table_S8A_Tau_Sensitivity_Class_Composition.tsv','Tau-sensitivity class composition','tau_truth_by_class.tsv','Tau-sensitivity Results subsection','','Governed tau authority.'),('S8B','Supplementary_Table_S8B_Tau_Sensitivity_SSI_Summary.tsv','Tau-sensitivity SSI summary','tau_summary.tsv','Tau-sensitivity Results subsection','Supplementary Figure S2','Mean SSI is the plotted endpoint.')];meta+=extra
readme=['# Final current FSF manuscript table package','','This directory contains manuscript-facing derivatives of frozen/current FSF authorities. It is not a new scientific authority and does not regenerate FSF results.','','| Table | Filename | Title/contents | Manuscript location | Related figure | Authority and guardrails |','|---|---|---|---|---|---|']+[f'| {a} | `{b}` | {c} | {e} | {f or "—"} | {d}. {g} |' for a,b,c,d,e,f,g in meta]+['','## Interpretation guardrails','','- Table 1 reports current stability-region architecture.','- S1 is preserved byte-for-byte.','- S2 is a deterministic gzip representation of the frozen annotation master. Annotation coverage and enrichment-universe sizes are distinct governed concepts: 13,467 features have at least one annotation source; 9,321 are GO-mapped; the KO background contains 7,178 features; and the KEGG pathway background contains 4,680 features.','- S5 retains raw KEGG identifiers.','- S7D is retired as redundant.','- The workbook is a human-review aid only.']
(OUT/'FINAL_TABLE_PACKAGE_README.md').write_text('\n'.join(readme)+'\n',encoding='utf-8')
mf=['table_id','filename','title','source_authority','row_count','column_count','sha256','manuscript_section','related_figure','counting_convention','status'];mr=[]
for a,b,c,d,e,f,g in meta:
 op=gzip.open if b.endswith('.gz') else open;dl=',' if b.endswith('.csv') else '\t'
 with op(OUT/b,'rt',encoding='utf-8',newline='') as h:x=list(csv.reader(h,delimiter=dl))
 mr.append(dict(zip(mf,[a,b,c,d,str(len(x)-1),str(len(x[0])),sha(OUT/b),e,f,g,'PASS'])))
wr(OUT/'SUPPLEMENTARY_TABLES_MANIFEST.tsv',mf,mr)
def load(p):
 op=gzip.open if str(p).endswith('.gz') else open;dl=',' if str(p).endswith('.csv') else '\t'
 with op(p,'rt',encoding='utf-8',newline='') as f:return list(csv.reader(f,delimiter=dl))
def val(x):
 if x=='':return None
 if re.fullmatch(r'-?(?:0|[1-9]\d*)',x):return int(x)
 if re.fullmatch(r'-?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?',x):return float(x)
 return x
embed={'Table 1':'Table1','S1':'S1_Metadata','S2':'Not embedded; use governed file','S3':'S3_Representative_Genes','S4':'Not embedded; use governed file','S5':'Not embedded; use governed file','S6':'S6_Baseline','S7A':'S7A_Noise_SSI','S7B':'S7B_Noise_Recovery','S7C':'Not embedded; use governed file','S8A':'S8A_Tau_Class','S8B':'S8B_Tau_SSI'};idx=[['table_id','filename','title','row_count','column_count','manuscript_section','workbook_sheet']]+[[r['table_id'],r['filename'],r['title'],int(r['row_count']),int(r['column_count']),r['manuscript_section'],embed[r['table_id']]] for r in mr]
sheets=[('Table1',meta[0][1]),('S1_Metadata',meta[1][1]),('S3_Representative_Genes',meta[3][1]),('S6_Baseline',meta[6][1]),('S7A_Noise_SSI',meta[7][1]),('S7B_Noise_Recovery',meta[8][1]),('S8A_Tau_Class',meta[10][1]),('S8B_Tau_SSI',meta[11][1]),('Manifest','SUPPLEMENTARY_TABLES_MANIFEST.tsv')]
wb=Workbook();wb.remove(wb.active);wb.properties.creator='openpyxl';wb.properties.created=wb.properties.modified=datetime(2026,9,27,2,0,39)
for name,data in [('Index',idx)]+[(n,load(OUT/f)) for n,f in sheets]:
 ws=wb.create_sheet(name)
 for row in data:ws.append([val(x) if isinstance(x,str) else x for x in row])
 for c in ws[1]:c.font=Font(b=True,color='00FFFFFF');c.fill=PatternFill('solid',fgColor='001F4E78');c.alignment=Alignment(vertical='top',wrap_text=True)
 for row in ws.iter_rows(min_row=2):
  for c in row:c.alignment=Alignment(vertical='top',wrap_text=True)
 ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;tab=Table(displayName=f'{name}_Table',ref=ws.dimensions);tab.tableStyleInfo=TableStyleInfo(name='TableStyleMedium2',showFirstColumn=False,showLastColumn=False,showRowStripes=True,showColumnStripes=False);ws.add_table(tab)
wb.save(OUT/'FSF_Final_Manuscript_Table_Index.xlsx')
expected={x[1] for x in meta}|{'FINAL_TABLE_PACKAGE_README.md','SUPPLEMENTARY_TABLES_MANIFEST.tsv','FSF_Final_Manuscript_Table_Index.xlsx'};actual={x.name for x in OUT.iterdir() if x.is_file()}
if actual!=expected:raise SystemExit(f'FILESET missing={expected-actual} extra={actual-expected}')
print('STAGE12_PRODUCER_EXECUTION=PASS')
