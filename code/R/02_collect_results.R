# Collect per-analysis MR results into the reported tables, recompute instrument strength (Supplementary Table 8 rows)
# and Steiger filtering from the harmonized data, and prepare figure inputs.
# Usage: Rscript code/R/02_collect_results.R
suppressPackageStartupMessages({library(data.table);library(TwoSampleMR)})
pkg <- normalizePath(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'..','..'))
man <- fread(file.path(pkg,'data/manifest.tsv'))
runs <- file.path(pkg,'results/runs/published_seed'); tabs <- file.path(pkg,'results/tables'); figd <- file.path(pkg,'results/figure_data')
dir.create(tabs,recursive=TRUE,showWarnings=FALSE); dir.create(figd,recursive=TRUE,showWarnings=FALSE)
stopifnot(all(file.exists(file.path(runs,man$analysis,'status.txt'))))
rd <- function(f) rbindlist(lapply(man$analysis,function(a) fread(file.path(runs,a,f),colClasses=list(character=intersect(c('global_p','distortion_p','Pvalue'),names(fread(file.path(runs,a,f),nrows=0)))))),fill=TRUE)
ce <- rd('causal_effects.tsv'); he <- rd('heterogeneity.tsv'); eg <- rd('egger_intercept.tsv'); pt <- rd('PRESSO_tests.tsv')
ce <- merge(ce,man[,.(analysis,analysis_type,published_table)],by='analysis')

## BH within method across the analyses of each family (main: 4; negative control: 8; leave-many-out: 2).
## As in the published tables, analyses without an estimate (MR-PRESSO outlier-corrected when no outliers) are counted in m.
ce[,BH_p:=p.adjust(p_value,'BH',n=.N),by=.(analysis_type,method)]
ce[status!='estimated',BH_p:=NA_real_]

## Selected (highlighted) method: Q P<.05 & intercept P<.05 -> MR-Egger; Q P<.05 -> weighted median;
## otherwise MR-PRESSO outlier-corrected if estimated, else IVW.
sel <- merge(he[method=='IVW',.(analysis,Q_p=Q_pval)],eg[,.(analysis,intercept_p=p_value)],by='analysis')
sel <- merge(sel,ce[method=='MR-PRESSO outlier-corrected',.(analysis,presso_corrected=status)],by='analysis')
sel[,selected_method:=fcase(Q_p<.05&intercept_p<.05,'MR-Egger',Q_p<.05,'Weighted median',
  intercept_p>=.05&presso_corrected=='estimated','MR-PRESSO outlier-corrected',intercept_p>=.05,'IVW',default=NA_character_)]
ce[,selected:=method==sel$selected_method[match(analysis,sel$analysis)]]
morder <- c('IVW','MR-Egger','Weighted median','MR-PRESSO raw','MR-PRESSO outlier-corrected','MR-RAPS')
ce[,method:=factor(method,morder)];setorder(ce,published_table,analysis,method);ce[,method:=as.character(method)]
ce[,analysis:=factor(analysis,man$analysis)];setorder(ce,analysis,method);ce[,analysis:=as.character(analysis)]
est <- ce[,.(published_table,analysis_type,analysis,method,nsnp,OR,SE_OR_delta,OR_lower95,OR_upper95,p_value,BH_p,log_OR=beta,SE_log_OR=se,ci_distribution,ci_df,status,selected,seed)]
fwrite(est,file.path(tabs,'MR_estimates_ST9_ST11_ST12.tsv'),sep='\t',na='NA')
fwrite(he,file.path(tabs,'Heterogeneity_Q_ST9_ST11_ST12.tsv'),sep='\t')
fwrite(eg,file.path(tabs,'Egger_intercept_ST9_ST11_ST12.tsv'),sep='\t')
fwrite(pt,file.path(tabs,'MR_PRESSO_tests.tsv'),sep='\t')
fwrite(sel,file.path(tabs,'Selected_method.tsv'),sep='\t')

## Instrument strength (Supplementary Table 8) recomputed from beta/SE/EAF with the published formulas.
TR <- list(short=list(s='lmm',P=.258,K=.258),long=list(s='lmm',P=.101,K=.101),MDD=list(s='logit',P=.218,K=.15))
r2_normal <- function(b,af,t){z<-dnorm(qnorm(1-t$K))
 if(t$s=='lmm'){bl<-b*t$P*(1-t$P);2*af*(1-af)*bl^2/(t$P*(1-t$P))*t$K^2*(1-t$K)^2/(t$P*(1-t$P)*z^2)}
 else{fp<-get_population_allele_frequency(af,rep(t$P,length(af)),exp(b),t$K);2*fp*(1-fp)*(b*t$K*(1-t$K)/z)^2}}
strength <- list(); steiger <- list()
for(i in seq_len(nrow(man))){a<-man[i];h<-fread(file.path(pkg,a$file),na.strings='NA');d<-h[included_in_MR==TRUE];if(a$analysis_type=='leave_many_out') d<-d[LMO_excluded==FALSE]
 F<-(d$beta.exposure/d$se.exposure)^2;r2<-r2_normal(d$beta.exposure,d$eaf.exposure,TR[[a$exposure]])
 stopifnot(max(abs(r2-d$rsq.exposure)/r2)<1e-10,max(abs(F-d$F_exposure)/F)<1e-10)
 strength[[i]]<-data.table(published_table='ST8',analysis_type=a$analysis_type,analysis=a$analysis,N_IV=nrow(d),Total_R2=sum(r2),R2_scale='Normal liability',
  Mean_F=mean(F),Median_F=median(F),Min_F=min(F),Max_F=max(F),N_F_lt10=sum(F<10))
 if(a$analysis_type=='main'){  # Steiger re-check on all candidate instruments
  rx<-r2_normal(h$beta.exposure,h$eaf.exposure,TR[[a$exposure]]);ry<-r2_normal(h$beta.outcome,h$eaf.outcome,TR[[a$outcome]])
  z<-(atanh(sqrt(rx))-atanh(sqrt(ry)))/sqrt(1/(h$effective_n.exposure-3)+1/(h$effective_n.outcome-3))
  steiger[[i]]<-data.table(analysis=a$analysis,SNP=h$SNP,R2_exposure=rx,R2_outcome=ry,steiger_dir=rx>ry,steiger_p=2*pnorm(-abs(z)),
   steiger_dir_in_data=h$steiger_dir,steiger_p_in_data=h$steiger_pval,included_in_MR=h$included_in_MR)}
}
strength<-rbindlist(strength);steiger<-rbindlist(steiger)
strength[analysis_type=='leave_many_out',published_table:='not in ST8 (leave-many-out set)']
steiger[,agree:=steiger_dir==steiger_dir_in_data & abs(steiger_p-steiger_p_in_data)<=1e-10*pmax(steiger_p,1e-300) & included_in_MR==steiger_dir]
stopifnot(all(steiger$agree))
fwrite(strength,file.path(tabs,'IV_strength_ST8.tsv'),sep='\t');fwrite(steiger,file.path(tabs,'Steiger_recheck_main.tsv'),sep='\t')

## Leave-many-out screening: excluded variants = variants with any predictor association P <= 5e-8 (complete lookup table)
lk <- fread(file.path(pkg,'data/leave_many_out/predictor_lookup_all_pairs.tsv'),na.strings='NA')
lsum <- lk[,.(SNPs_screened=uniqueN(SNP),predictors=uniqueN(predictor_id),pairs=.N,pairs_found=sum(lookup_status=='FOUND'),
  pairs_not_found=sum(lookup_status!='FOUND'),gws_associations=sum(genome_wide_significant),SNPs_excluded=uniqueN(SNP[genome_wide_significant])),by=direction]
for(x in lsum$direction){h<-fread(file.path(pkg,man[analysis==x]$file),na.strings='NA')
 stopifnot(setequal(h[LMO_excluded==TRUE]$SNP,unique(lk[direction==x & genome_wide_significant==TRUE]$SNP)),
           all(lk[direction==x & genome_wide_significant==TRUE]$p_value<=5e-8),all(lk[direction==x & lookup_status=='FOUND' & genome_wide_significant==FALSE]$p_value>5e-8))}
lsum[,SNPs_retained:=SNPs_screened-SNPs_excluded]
fwrite(lsum,file.path(tabs,'Leave_many_out_screening_summary.tsv'),sep='\t')

## Figure inputs for the four main diagnostic figures: single-SNP Wald ratios and IVW leave-one-out.
for(a in man[analysis_type=='main']$analysis){d<-fread(file.path(runs,a,'input_IVs.tsv'))
 dd<-as.data.frame(d[,.(SNP,beta.exposure,beta.outcome,se.exposure,se.outcome,exposure=a,outcome=a,id.exposure=a,id.outcome=a,mr_keep=TRUE)])
 z<-as.data.table(mr_singlesnp(dd,all_method=character()));loo<-as.data.table(mr_leaveoneout(dd))
 z[,`:=`(lower=b-qnorm(.975)*se,upper=b+qnorm(.975)*se)];loo[,`:=`(lower=b-qnorm(.975)*se,upper=b+qnorm(.975)*se)]
 stopifnot(nrow(z)==nrow(d),nrow(loo)==nrow(d)+1,abs(loo[SNP=='All']$b-ce[analysis==a & method=='IVW']$beta)<1e-12,
           max(abs(z$b-d$beta.outcome/d$beta.exposure))<1e-12)
 fwrite(z[,.(SNP,b,se,p,lower,upper)],file.path(figd,paste0(a,'_single_SNP.tsv')),sep='\t')
 fwrite(loo[,.(SNP,b,se,p,lower,upper)],file.path(figd,paste0(a,'_leave_one_out.tsv')),sep='\t')}
message('Tables written to results/tables; figure inputs to results/figure_data')
