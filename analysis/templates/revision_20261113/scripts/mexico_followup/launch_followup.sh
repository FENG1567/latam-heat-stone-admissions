#!/usr/bin/env bash
set -euo pipefail

mode="${1:?usage: launch_followup.sh A|B}"
case "$mode" in A|B) ;; *) echo "mode must be A or B" >&2; exit 2 ;; esac

root=/data/admin/fengshuo/latam_stone_repair_20261112/revision_20261113
rlib=/data/admin/fengshuo/latam_stone_repair_20261112/Rlib
script="$root/scripts/mexico_followup/refit_mx2023_followup.R"
logdir="$root/logs/mexico_followup"
mkdir -p "$logdir"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
stem="mexico_followup_${mode}_${stamp}"
log="$logdir/${stem}.log"
exitfile="$logdir/${stem}.exitcode"
pidfile="$logdir/${stem}.pid"

nohup bash -c '
  mode="$1"; exitfile="$2"; root="$3"; script="$4"; rlib="$5"
  export R_LIBS_USER="$rlib"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLAS_NUM_THREADS=1
  /usr/bin/Rscript --vanilla "$script" "--root=$root" "--mode=$mode"
  rc=$?
  printf "%s\n" "$rc" > "$exitfile"
  exit "$rc"
' _ "$mode" "$exitfile" "$root" "$script" "$rlib" >"$log" 2>&1 </dev/null &
pid=$!
printf "%s\n" "$pid" > "$pidfile"
printf "PID=%s\nLOG=%s\nEXIT=%s\n" "$pid" "$log" "$exitfile"

