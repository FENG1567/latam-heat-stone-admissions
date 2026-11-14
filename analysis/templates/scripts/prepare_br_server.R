#!/usr/bin/env Rscript
# Complete-calendar implementation of prepare_panels.py using the same weather
# bridge. It avoids transferring a large CSV over the slow SSH connection.
.libPaths(c('/data/admin/fengshuo/latam_stone_repair_20261112/Rlib',.libPaths()))
library(data.table)
setDTthreads(4)
base <- '/data/admin/fengshuo/latam_stone_repair_20261112'
source_path <- '/data/admin/fengshuo/latam_stone_20261109/data/processed/brazil_analysis.csv'
d <- fread(source_path,colClasses=list(character='hrid'))
setnames(d,c('hrid','K80','K81','N20_N23','sum3','pop_region'),
         c('unit','y_K80','y_K81','y_urinary','y_sum3','population'))
keep <- c('unit','date','tmean','tmax','population','y_K80','y_K81','y_urinary','y_sum3')
d <- d[,..keep]
d[,date:=as.IDate(date)]
bridge <- fread(file.path(base,'br_weather_bridge.csv'),colClasses=list(character='unit'))
bridge[,date:=as.IDate(date)]
stopifnot(!anyDuplicated(d[,.(unit,date)]),!anyDuplicated(bridge[,.(unit,date)]))
d <- rbindlist(list(d,bridge),fill=TRUE)
stopifnot(!anyDuplicated(d[,.(unit,date)]))
grid <- CJ(unit=sort(unique(d$unit)),date=seq(as.IDate('2005-12-12'),as.IDate('2020-12-31'),by='day'))
d <- merge(grid,d,by=c('unit','date'),all.x=TRUE,sort=TRUE)
stopifnot(d[,all(diff(as.integer(date))==1L),by=unit]$V1)
for(k in c(3L,7L,14L,21L)) d[,(paste0('E',k)):=frollmean(tmean,n=k,align='right',na.rm=FALSE),by=unit]
d[,Emax3:=frollmean(tmax,n=3L,align='right',na.rm=FALSE),by=unit]
d <- d[date>=as.IDate('2006-01-01')]
d[,`:=`(dow=as.integer(format(date,'%u')),month=as.integer(format(date,'%m')),year=as.integer(format(date,'%Y')),country='BR')]
d[,population:=nafill(nafill(population,type='locf'),type='nocb'),by=.(unit,year)]
stopifnot(nrow(d)==2465550L,uniqueN(d$unit)==450L,sum(is.na(d$y_K80))==1135L,
          sum(is.na(d$tmean))==0L,sum(d$y_K80,na.rm=TRUE)==2173811L,
          sum(d$y_urinary,na.rm=TRUE)==1108174L)
destination <- file.path(base,'data','br_calendar.csv.gz')
if(file.exists(destination)) stopifnot(file.rename(destination,paste0(destination,'.partial_transfer')))
fwrite(d,destination,compress='gzip')
saveRDS(d,file.path(base,'data','br_calendar_server.rds'),compress='gzip')
cat('prepared BR full calendar',nrow(d),'rows, no missing weather; outcome gaps retained\n')
