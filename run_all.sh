#!/bin/bash
# run_all.sh: refits all data analyses in order (udder, bcos, Figure 1).
# The simulation (sim/) is run separately.
#
# Start from code/, with the R environment active:
#   setsid nohup bash run_all.sh > run_all_log.txt 2>&1 < /dev/null &
# Progress: cat run_all_log.txt
#
# If WAIT_SIM=1, waits until no sim_table1.R / sim_shapes.R process is
# running before starting, so the fits do not compete for cores.

cd "$(dirname "$0")"
ROOT=$(pwd)

if [ "${WAIT_SIM:-0}" = "1" ]; then
  while pgrep -f "sim_table1.R|sim_shapes.R" > /dev/null; do
    echo "$(date '+%F %T') waiting for the simulation to finish"; sleep 600
  done
fi

run() {   # run <dir> <script>
  cd "$ROOT/$1"
  local log="${2%.R}_log.txt"
  echo "$(date '+%F %T') START $1/$2"
  local t0=$(date +%s)
  if Rscript "$2" > "$log" 2>&1; then
    echo "$(date '+%F %T') DONE  $1/$2 ($(( $(date +%s) - t0 )) s)"
  else
    echo "$(date '+%F %T') FAIL  $1/$2 -- see $1/$log"; exit 1
  fi
}

run udder 02_fit_frailty.R
run udder fig2_frailty.R
run udder sens_prior_theta.R
run udder sens_prior_sigma.R
run udder 04_fit_baselines.R
run udder 05_tables_95.R
run udder sens_resolution.R
run bcos  02_fit_baselines.R
run bcos  03_sens_prior.R
run figures fig1_baselines.R

cd "$ROOT"; Rscript -e 'suppressMessages({library(cmdstanr); library(posterior); library(loo); library(splines2); library(survivalMPL); library(UdderQuarterInfectionData)}); print(sessionInfo()); cat("CmdStan", cmdstanr::cmdstan_version(), "(scripts set 2.37.0)\n")' > session_info.txt 2>&1
echo "$(date '+%F %T') ALL DONE"
