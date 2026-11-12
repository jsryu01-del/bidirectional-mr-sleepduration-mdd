# Check that required packages are installed and record the versions used for this run.
# Usage: Rscript code/R/00_check_environment.R
pkg <- normalizePath(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'..','..'))
need <- c('data.table','TwoSampleMR','MRPRESSO','mr.raps')
miss <- need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop('Missing R packages: ',paste(miss,collapse=', '),'. See README (Installation).')
used <- read.delim(file.path(pkg,'environment/package_versions_original_analysis.tsv'))
now <- data.frame(package=c('R',need),version=c(paste(R.version$major,R.version$minor,sep='.'),vapply(need,function(p) as.character(packageVersion(p)),'')))
m <- merge(used,now,by='package',suffixes=c('_original','_this_run'))
m$same <- m$version_original==m$version_this_run
dir.create(file.path(pkg,'results'),showWarnings=FALSE)
write.table(m,file.path(pkg,'results/package_versions_this_run.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(pkg,'results/sessionInfo_this_run.txt'))
print(m)
if(!all(m$same)) message('NOTE: package versions differ from the original analysis; deterministic estimates should agree, ',
 'but Monte Carlo methods (weighted median bootstrap, MR-PRESSO) may differ slightly.')
