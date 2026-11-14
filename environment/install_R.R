# Install required versions observed in the captured analysis environment.
args <- commandArgs(trailingOnly=TRUE)
tab <- read.delim('environment/R_packages.tsv',stringsAsFactors=FALSE)
required <- c('data.table','jsonlite','metafor','sandwich','digest')
if(!requireNamespace('remotes',quietly=TRUE)) install.packages('remotes',repos='https://cloud.r-project.org')
for(pkg in required){
 row <- tab[tab$Package==pkg,,drop=FALSE]
 if(nrow(row)!=1) stop(paste('No unique recorded version for',pkg))
 if(!requireNamespace(pkg,quietly=TRUE) || as.character(packageVersion(pkg))!=row$Version)
  remotes::install_version(pkg,version=row$Version,repos='https://cloud.r-project.org',upgrade='never')
}
cat('Required R packages match the captured version list.\n')
