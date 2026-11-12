# Compare reproduced results with the values printed in Supplementary Tables 8, 9 and 12 (main and negative-control
# analyses) and 10-11 (leave-many-out).
# Usage: Rscript code/R/03_compare_with_published.R      (exit status 1 if any value differs)
suppressPackageStartupMessages(library(data.table))
pkg <- normalizePath(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'..','..'))
TOL <- 1e-8   # relative tolerance for ST8, ST9 and ST11, whose values are stored at full double precision
# ST12 is printed rounded (e.g. "1.014 (0.956, 1.076)", "5.39×10−6"): a value matches if it rounds to the printed value,
# i.e. |reproduced - printed| <= half a unit of the last printed digit.
printed_tol <- function(txt) vapply(txt,function(s){ if(is.na(s)) return(NA_real_)
 s <- gsub(intToUtf8(0x2212),'-',enc2utf8(s),fixed=TRUE); tx <- paste0(intToUtf8(0xD7),'10')
 if(grepl(tx,s,fixed=TRUE)){p<-strsplit(s,tx,fixed=TRUE)[[1]];m<-p[1];e<-as.numeric(p[2])} else {m<-s;e<-0}
 d <- if(grepl('.',m,fixed=TRUE)) nchar(sub('^[^.]*\\.','',m)) else 0
 0.5*10^(e-d)*(1+1e-9)},numeric(1),USE.NAMES=FALSE)
pub <- fread(file.path(pkg,'data/published_values/published_values.tsv'),na.strings='NA',encoding='UTF-8')
tabs <- file.path(pkg,'results/tables')
est <- fread(file.path(tabs,'MR_estimates_ST9_ST11_ST12.tsv'),na.strings='NA')
he <- fread(file.path(tabs,'Heterogeneity_Q_ST9_ST11_ST12.tsv')); eg <- fread(file.path(tabs,'Egger_intercept_ST9_ST11_ST12.tsv'))
st <- fread(file.path(tabs,'IV_strength_ST8.tsv')); sel <- fread(file.path(tabs,'Selected_method.tsv'))
# not-estimated rows are printed as empty cells
for(cc in c('nsnp','OR','SE_OR_delta','OR_lower95','OR_upper95','p_value','BH_p','log_OR','SE_log_OR')) est[status!='estimated',(cc):=NA]
long <- rbind(
 melt(est,id.vars=c('published_table','analysis','method'),measure.vars=c('nsnp','OR','SE_OR_delta','OR_lower95','OR_upper95','p_value','BH_p','log_OR','SE_log_OR'),variable.name='field',value.name='reproduced')[,table:=published_table][,published_table:=NULL],
 melt(he,id.vars=c('analysis','method'),measure.vars=c('Q','Q_df','Q_pval'),variable.name='field',value.name='reproduced')[,table:=NA_character_],
 melt(eg[,.(analysis,method='Egger intercept',intercept,se,ci_lower,ci_upper,p_value)],id.vars=c('analysis','method'),variable.name='field',value.name='reproduced')[,table:=NA_character_],
 melt(st[,.(analysis,method='IV strength',N_IV,Total_R2,Mean_F,Median_F,Min_F,Max_F,N_F_lt10)],id.vars=c('analysis','method'),variable.name='field',value.name='reproduced')[,table:='ST8'],
 use.names=TRUE)
long[,field:=as.character(field)]
cmp <- merge(pub[field!='selected_method' & table!='ST10'],long[,!'table'],by=c('analysis','method','field'),all.x=TRUE)
cmp[,rel_diff:=fifelse(is.na(value)&is.na(reproduced),0,fifelse(is.na(value)|is.na(reproduced),Inf,abs(value-reproduced)/pmax(abs(value),1e-300)))]
cmp[,abs_diff:=fifelse(is.na(value)&is.na(reproduced),0,fifelse(is.na(value)|is.na(reproduced),Inf,abs(value-reproduced)))]
cmp[,printed_tolerance:=fifelse(table=='ST12',printed_tol(text),NA_real_)]
cmp[,status:=fifelse(fifelse(table=='ST12',abs_diff<=printed_tolerance,rel_diff<=TOL),'Match','DIFFERENT')]
cmp[table=='ST12' & is.na(value) & is.na(reproduced),status:='Match']
s <- merge(pub[field=='selected_method',.(table,analysis,published=text)],sel[,.(analysis,reproduced=selected_method)],by='analysis',all=TRUE)
# ST10: printed excluded SNP-predictor associations vs complete lookup table
lk <- fread(file.path(pkg,'data/leave_many_out/predictor_lookup_all_pairs.tsv'),na.strings='NA')[genome_wide_significant==TRUE]
p10 <- pub[table=='ST10',.(direction=analysis,SNP=method,predictor_id=text,published=value)]
x10 <- merge(p10,lk[,.(direction,SNP,predictor_id,reproduced=p_value)],by=c('direction','SNP','predictor_id'),all=TRUE)
x10[,status:=fifelse(!is.na(published)&!is.na(reproduced)&abs(published-reproduced)<=1e-6*abs(reproduced),'Match','DIFFERENT')]
fwrite(x10,file.path(pkg,'results/comparison_ST10_excluded_associations.tsv'),sep='\t',na='NA')
s[,status:=fifelse(!is.na(published)&published==reproduced,'Match','DIFFERENT')]
fwrite(cmp[,.(table,analysis,method,field,published=value,printed=text,reproduced,rel_diff,abs_diff,printed_tolerance,status)],file.path(pkg,'results/comparison_with_published.tsv'),sep='\t',na='NA')
fwrite(s,file.path(pkg,'results/comparison_selected_method.tsv'),sep='\t')
# max_rel_diff: full-precision tables; max_diff_over_printed_tol: ST12 (<= 1 means the reproduced value rounds to the printed one)
summ <- rbind(cmp[,.(comparison=fifelse(table[1]=='ST12','printed precision (rounded table)',paste('relative tolerance',TOL)),n_values=.N,n_match=sum(status=='Match'),n_different=sum(status!='Match'),
                   max_rel_diff=if(table[1]=='ST12') NA_real_ else max(rel_diff),
                   max_diff_over_printed_tol=if(table[1]=='ST12') max(fifelse(is.na(value)&is.na(reproduced),0,abs_diff/printed_tolerance)) else NA_real_),by=table],
              s[,.(table='Selected method (highlighted/bold row)',comparison='exact',n_values=.N,n_match=sum(status=='Match'),n_different=sum(status!='Match'),max_rel_diff=NA_real_,max_diff_over_printed_tol=NA_real_)],
              x10[,.(table='ST10 excluded SNP-predictor associations',comparison='relative tolerance 1e-6',n_values=.N,n_match=sum(status=='Match'),n_different=sum(status!='Match'),max_rel_diff=NA_real_,max_diff_over_printed_tol=NA_real_)])
fwrite(summ,file.path(pkg,'results/comparison_summary.tsv'),sep='\t')
print(summ)
if(any(summ$n_different>0)){print(cmp[status!='Match']);print(s[status!='Match']);print(x10[status!='Match']);quit(status=1)}
message('All reproduced values match the published Supplementary Tables 8-12 (ST8, ST9, ST11: relative tolerance ',TOL,
        '; ST12: printed precision; ST10: relative tolerance 1e-6).')
