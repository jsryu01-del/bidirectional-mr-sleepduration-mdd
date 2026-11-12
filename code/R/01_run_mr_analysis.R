# Run the MR analyses for one analysis from the submitted harmonized dataset.
# Usage: Rscript code/R/01_run_mr_analysis.R <analysis> [seed_offset]
#   <analysis>  : an `analysis` value from data/manifest.tsv (main, negative_control or leave_many_out)
#   seed_offset : optional integer added to the published seed (0 = published results; used only for Monte Carlo sensitivity)
# Methods and settings:
#   IVW (multiplicative random effects; TwoSampleMR::mr_ivw), MR-Egger (TwoSampleMR::mr_egger_regression; t(n-2) CI),
#   weighted median (TwoSampleMR::mr_weighted_median; 10,000 bootstraps), MR-PRESSO (MRPRESSO::mr_presso; 10,000 simulations;
#   alpha 0.05; t(n_used-1) CI) for all analyses; MR-RAPS (mr.raps::mr.raps.mle; over-dispersion, Huber loss, no pruning,
#   sandwich SE) additionally for the main and leave-many-out analyses.
# Monte Carlo methods (weighted median, MR-PRESSO) depend on the seed AND on SNP row order; rows are analysed in `mr_row_order`.
suppressPackageStartupMessages({library(data.table);library(TwoSampleMR)})
args <- commandArgs(TRUE); target <- args[1]; offset <- if(length(args)>1) as.integer(args[2]) else 0L
pkg <- normalizePath(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'..','..'))
man <- fread(file.path(pkg,'data/manifest.tsv')); a <- man[analysis==target]; stopifnot(nrow(a)==1)
tag <- if(offset==0L) 'published_seed' else paste0('seed_offset_',offset)
out <- file.path(pkg,'results/runs',tag,target); dir.create(out,recursive=TRUE,showWarnings=FALSE)
h <- fread(file.path(pkg,a$file),na.strings='NA')
# leave-many-out analyses use the main instrument set minus variants associated (P <= 5e-8) with candidate predictors
d <- if(a$analysis_type=='leave_many_out') h[included_in_MR==TRUE & LMO_excluded==FALSE][order(mr_row_order)] else h[included_in_MR==TRUE][order(mr_row_order)]
stopifnot(nrow(d)==a$n_included,!anyDuplicated(d$SNP),all(d$harmonisation_mr_keep),!is.unsorted(d$mr_row_order),
          a$analysis_type=='leave_many_out' || identical(d$mr_row_order,seq_len(nrow(d))),
          all(is.finite(d$beta.exposure)),all(is.finite(d$beta.outcome)),all(d$se.exposure>0),all(d$se.outcome>0),
          identical(d$effect_allele.exposure,d$effect_allele.outcome))
n <- nrow(d); seed <- a$seed+offset; set.seed(seed)
fwrite(d,file.path(out,'input_IVs.tsv'),sep='\t')
warn <- list()
run <- function(label,expr) withCallingHandlers(expr,warning=function(w){warn[[length(warn)+1]]<<-data.table(analysis=target,method=label,message=conditionMessage(w));invokeRestart('muffleWarning')})
x<-d$beta.exposure;y<-d$beta.outcome;sx<-d$se.exposure;sy<-d$se.outcome
par<-default_parameters();par$nboot<-10000L
ivw<-run('IVW',mr_ivw(x,y,sx,sy,par));eg<-run('MR-Egger',mr_egger_regression(x,y,sx,sy,par))
wm<-run('Weighted median',mr_weighted_median(x,y,sx,sy,par))
use_raps<-a$analysis_type %in% c('main','leave_many_out')
ra<-if(use_raps) run('MR-RAPS',mr.raps::mr.raps.mle(x,y,sx,sy,over.dispersion=TRUE,loss.function='huber',diagnostics=FALSE,pruning=FALSE,se.method='sandwich')) else NULL
row<-function(method,b,se,p,nsnp=n,df=Inf,status='estimated'){crit<-if(is.finite(df)) qt(.975,df) else qnorm(.975)
 data.table(analysis=target,seed,method,nsnp,beta=b,se,ci_lower=b-crit*se,ci_upper=b+crit*se,OR=exp(b),OR_lower95=exp(b-crit*se),OR_upper95=exp(b+crit*se),
  SE_OR_delta=exp(b)*se,p_value=as.numeric(p),ci_distribution=if(is.finite(df)) 't' else 'normal',ci_df=df,status)}
res<-list(row('IVW',ivw$b,ivw$se,ivw$pval),row('MR-Egger',eg$b,eg$se,eg$pval,df=n-2),row('Weighted median',wm$b,wm$se,wm$pval))
if(use_raps) res[[length(res)+1]]<-row('MR-RAPS',ra$beta.hat,ra$beta.se,ra$beta.p.value)
message(target,': standard methods done; MR-PRESSO (10,000 simulations, seed ',seed,')')
pr<-run('MR-PRESSO',MRPRESSO::mr_presso(BetaOutcome='beta.outcome',BetaExposure='beta.exposure',SdOutcome='se.outcome',SdExposure='se.exposure',
 data=as.data.frame(d),OUTLIERtest=TRUE,DISTORTIONtest=TRUE,NbDistribution=10000,SignifThreshold=.05,seed=seed))
tests<-pr[['MR-PRESSO results']];global<-tests[['Global Test']];dist<-tests[['Distortion Test']]
idx<-if(!is.null(dist) && is.numeric(dist[['Outliers Indices']])) dist[['Outliers Indices']] else integer()
main<-pr[['Main MR results']]
for(j in seq_len(nrow(main))){corr<-grepl('correct',as.character(main[j,2]),ignore.case=TRUE)
 b<-as.numeric(main[j,3]);se<-as.numeric(main[j,4]);p<-as.numeric(main[j,6]);cnt<-if(corr) n-length(idx) else n
 res[[length(res)+1]]<-row(if(corr)'MR-PRESSO outlier-corrected' else 'MR-PRESSO raw',b,se,p,cnt,df=cnt-1,
  status=if(is.finite(b)) 'estimated' else if(!length(idx)) 'not_applicable_no_outliers' else 'not_estimable')}
fwrite(rbindlist(res),file.path(out,'causal_effects.tsv'),sep='\t')
fwrite(rbindlist(lapply(c('IVW','MR-Egger'),function(m){z<-if(m=='IVW') ivw else eg;data.table(analysis=target,method=m,Q=z$Q,Q_df=z$Q_df,Q_pval=z$Q_pval)})),file.path(out,'heterogeneity.tsv'),sep='\t')
fwrite(data.table(analysis=target,nsnp=n,intercept=eg$b_i,se=eg$se_i,ci_lower=eg$b_i-qt(.975,n-2)*eg$se_i,ci_upper=eg$b_i+qt(.975,n-2)*eg$se_i,p_value=eg$pval_i,
 scale='outcome log odds; not exponentiated'),file.path(out,'egger_intercept.tsv'),sep='\t')
fwrite(data.table(analysis=target,seed,nsnp=n,NbDistribution=10000,global_RSS=as.numeric(global[['RSSobs']]),global_p=as.character(global[['Pvalue']]),
 distortion_percent=if(is.null(dist[['Distortion Coefficient']])) NA_real_ else as.numeric(dist[['Distortion Coefficient']]),
 distortion_p=if(is.null(dist[['Pvalue']])) NA_character_ else as.character(dist[['Pvalue']]),n_outliers=length(idx),outlier_SNPs=paste(d$SNP[idx],collapse=';')),file.path(out,'PRESSO_tests.tsv'),sep='\t')
ot<-if(!is.null(tests[['Outlier Test']])) as.data.table(tests[['Outlier Test']]) else data.table(RSSobs=NA_real_,Pvalue=NA_character_)[rep(1,n)]
ot[,`:=`(analysis=target,SNP=d$SNP,is_outlier=seq_len(n) %in% idx)];fwrite(ot,file.path(out,'PRESSO_outlier_tests.tsv'),sep='\t')
fwrite(if(length(warn)) rbindlist(warn) else data.table(analysis=character(),method=character(),message=character()),file.path(out,'warnings.tsv'),sep='\t')
saveRDS(Filter(Negate(is.null),list(ivw=ivw,egger=eg,weighted_median=wm,raps=ra,presso=pr)),file.path(out,'raw_fits.rds'))
writeLines('SUCCESS',file.path(out,'status.txt'));message(target,' (',tag,'): complete')
