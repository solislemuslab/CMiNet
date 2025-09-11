## Directory structure

```text
Evaluation/
├─ Git_Gut/
│  ├─ evaluation_all/                  # consensus (all methods) – gut data
│  │  ├─ 01_bootstrap_amgut_counts.R
│  │  ├─ 02_run_cminet_on_bootstraps.R
│  │  ├─ 03_eval_bootstrap_thresholds.R
│  │  ├─ 04_gridsearch_m0_mstar_cstar.R
│  │  ├─ 05_plot_bootstrap_consensus_heatmap.R
│  │  ├─ 06_Jeffreys_CIs.R
│  │  ├─ Figure_boxplot/               # scripts to render boxplots
│  │  │  ├─ 01_bootstrap_CI.R
│  │  │  ├─ 02_boxplot_CI.R
│  │  │  ├─ 03_boxplot_fscore.R
│  │  │  └─ read_sample_data/          # CSVs used by plotting scripts
│  │  └─ sample_result_allmethod/      # precomputed outputs (RData/CSV)
│  ├─ evaluation_conditional_gut/      # conditional-dependence variant (gut)
│  │  ├─ 01_run_cminet_on_bootstraps_conditional.R
│  │  ├─ 02_eval_bootstrap_thresholds_conditional.R
│  │  └─ sample_result_conditional_gut/   # precomputed outputs (RData/CSV)
│  └─ evaluation_correlation_gut/      # correlation-only variant (gut)
│     ├─ 01_run_cminet_on_bootstraps_correlation.R
│     ├─ 02_eval_bootstrap_thresholds_correlation.R
│     └─ sample_result_correlation_gut/   # precomputed outputs (RData/CSV)
│
├─ Git_Soil/
│  ├─ evaluation_all/                  # consensus (all methods) – soil data
│  │  ├─ 01_bootstrap_soil_counts.R
│  │  ├─ 02_run_cminet_on_bootstraps.R
│  │  ├─ 03_eval_bootstrap_thresholds.R
│  │  ├─ 05_plot_bootstrap_consensus_heatmap.R
│  │  ├─ 06_Jeffreys_CIs.R
│  │  ├─ 4-Bootstrap-based edge confidence.R
│  │  ├─ Figure_boxplot/               # same structure as gut
│  │  └─ sample_result_allmethod/
│  ├─ evaluation_conditional_soil/
│  │  ├─ 01_run_cminet_on_bootstraps_conditional.R
│  │  └─ 02_eval_bootstrap_thresholds_conditional.R
│  └─ evaluation_correlation_soil/
│     ├─ 01_run_cminet_on_bootstraps_correlation.R
│     └─ 02_eval_bootstrap_thresholds_correlation.R
│
└─ Git_simulation/
   └─ simulatedresult.R                # helper for simulated results
