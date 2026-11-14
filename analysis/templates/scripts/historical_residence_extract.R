#!/usr/bin/env Rscript
# Count admissions by ADMISSION month (DT_INTER), reading the target year's
# processing files plus early files of the following year (late submissions).
# Usage: Rscript extract_year_by_admission.R YEAR [NEXT_MONTHS] [NCORES]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: extract_year_by_admission.R YEAR [NEXT_MONTHS] [NCORES]")
year <- as.integer(args[1])
next_months <- if (length(args) >= 2) as.integer(args[2]) else 6L
ncores <- if (length(args) >= 3) as.integer(args[3]) else 16L
.libPaths(c("/data/admin/fengshuo/climate_biliary_20261102/env/Rlib", .libPaths()))
suppressPackageStartupMessages(library(read.dbc))
suppressPackageStartupMessages(library(parallel))
root <- "/data/admin/fengshuo/climate_biliary_20261102"
dirs <- c(file.path(root, "data/raw/SIHSUS/RD", year),
          file.path(root, "data/raw/SIHSUS/RD", year + 1))
files <- unlist(lapply(dirs, function(d) list.files(d, pattern="\\.dbc$", recursive=TRUE, full.names=TRUE)))
if (year < 2024) {
  pat <- sprintf("RD[A-Z]{2}%02d(0[1-9]|1[0-%d])\\.dbc$", (year + 1) %% 100, next_months)
  # simpler: keep next-year files whose YYMM month <= next_months
  keep <- grepl(sprintf("RD[A-Z]{2}%02d", (year + 1) %% 100), files)
  mm <- as.integer(sub(".*RD[A-Z]{2}[0-9]{2}([0-9]{2})\\.dbc$", "\\1", files))
  files <- files[!keep | (keep & mm <= next_months)]
}
cat("year", year, "files", length(files), "cores", ncores, "\n")
uf_map <- c("11"="RO","12"="AC","13"="AM","14"="RR","15"="PA","16"="AP","17"="TO","21"="MA","22"="PI","23"="CE","24"="RN","25"="PB","26"="PE","27"="AL","28"="SE","29"="BA","31"="MG","32"="ES","33"="RJ","35"="SP","41"="PR","42"="SC","43"="RS","50"="MS","51"="MT","52"="GO","53"="DF")
one <- function(f) {
  d <- tryCatch(read.dbc(f), error=function(e) NULL)
  if (is.null(d) || nrow(d)==0) return(NULL)
  diag <- substr(as.character(d$DIAG_PRINC),1,3)
  dt <- as.Date(as.character(d$DT_INTER), format="%Y%m%d")
  muni <- suppressWarnings(sprintf("%06d", as.integer(as.character(d$MUNIC_RES))))
  uf <- substr(muni,1,2)
  ok <- !is.na(dt) & format(dt,"%Y")==as.character(year) & !is.na(diag) & uf %in% names(uf_map)
  if (!any(ok)) return(NULL)
  d2 <- data.frame(uf=uf[ok], muni=muni[ok], month=format(dt[ok],"%Y-%m"), date=as.character(dt[ok]), diag=diag[ok],
                   total=1L,
                   K80=as.integer(diag[ok]=="K80"), K81=as.integer(diag[ok]=="K81"),
                   K85=as.integer(diag[ok]=="K85"),
                   N20_N23=as.integer(diag[ok] %in% c("N20","N21","N22","N23")),
                   stringsAsFactors=FALSE)
  aggregate(cbind(total,K80,K81,K85,N20_N23) ~ uf+muni+month+date, data=d2, FUN=sum)
}
res <- mclapply(files, one, mc.cores=ncores, mc.preschedule=FALSE)
res <- res[!vapply(res, is.null, logical(1))]
all <- do.call(rbind, res)
all$uf_abbr <- unname(uf_map[all$uf])
cat("municipalities", length(unique(all$muni)), "rows", nrow(all), "\n")
daily <- aggregate(cbind(total,K80,K81,K85,N20_N23) ~ uf+uf_abbr+muni+date, data=all, FUN=sum)
write.csv(daily, gzfile(file.path(root, "data/interim", paste0("municipality_day_counts_", year, ".csv.gz"))), row.names=FALSE)
cat("daily rows", nrow(daily), "\n")
muni <- aggregate(cbind(total,K80,K81,K85,N20_N23) ~ uf+uf_abbr+muni+month, data=all, FUN=sum)
write.csv(muni, file.path(root, "data/interim", paste0("municipality_month_byadmit_", year, ".csv")), row.names=FALSE)
state <- aggregate(cbind(total,K80,K81,K85,N20_N23) ~ uf+uf_abbr+month, data=all, FUN=sum)
grid <- expand.grid(uf_abbr=sort(unique(unname(uf_map))), month=sprintf("%d-%02d",year,1:12), stringsAsFactors=FALSE)
grid <- merge(grid, state[,c("uf_abbr","month","total","K80","K81","K85","N20_N23")], by=c("uf_abbr","month"), all.x=TRUE)
for (v in c("total","K80","K81","K85","N20_N23")) grid[[v]][is.na(grid[[v]])] <- 0L
grid <- grid[order(grid$uf_abbr, grid$month),]
write.csv(grid, file.path(root, "data/interim", paste0("state_month_byadmit_", year, ".csv")), row.names=FALSE)
cat("state rows", nrow(grid), "total", sum(grid$total), "N20_N23", sum(grid$N20_N23),
    "K80", sum(grid$K80), "K81", sum(grid$K81), "K85", sum(grid$K85), "\n")
