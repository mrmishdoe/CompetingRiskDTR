# Setup 
library(tidyverse)
library(rms)
library(parallel)
source("helper_functions_analysis.R")

# Load data (not provided in repository)
data = readRDS("Data/Temp/data.rds")

################# Weighted and Greedy ITRs

# Model specification
treat.mod = DonorHCV ~ age_n + gender_num + RECAA + RECDM + RecHCV + immuno_group + log_size
cens.mod = delta~ age_n + gender_num + RECAA + RECDM +RecHCV + immuno_group+ DonorHCV + DON_TY + donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

out.mod = log(timeto)~ gender_num + RECAA + immuno_group + DonorHCV*RecHCV + DonorHCV*DON_TY + DonorHCV*age_n+  donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

cause.mod = ind~ age_n + gender_num + RECAA + RECDM +RecHCV + immuno_group+ DonorHCV + DON_TY + donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

models = list(treat.mod = treat.mod, cens.mod = cens.mod, out.mod = out.mod,cause.mod  =cause.mod)

set.seed(2024)
n_boot = 1000
boot.results = AFT.boot(n_boot, data, models)

#summary
summary_boot = summaryAFT(data, list(est = boot.results$est, boot.est = boot.results$boot.est))

#Plotting
dat_w = summary_boot$data_w
dat_g = summary_boot$data_g

plots_w = plotAFT(dat_w, "(a) Weighted regime")
lots_g = plotAFT(dat_g, "(a) Greedy regime")
cause_plot=ggarrange(plots_w$p1, plots_g$p1 ,plots_w$p2, plots_g$p2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
overall_plot = ggarrange(plots_w$p3, plots_g$p3, ncol = 2, nrow = 1, common.legend = T, legend = "right")


################# Cause-specific and composite outcome ITRs
treat.mod = DonorHCV ~ age_n + gender_num + RECAA + RECDM + RecHCV + immuno_group + log_size
cens.mod = delta~ age_n + gender_num + RECAA + RECDM +RecHCV + immuno_group+ DonorHCV + DON_TY + donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

out.mod = log(timeto)~ gender_num + RECAA + immuno_group + DonorHCV*RecHCV + DonorHCV*DON_TY + DonorHCV*age_n+  donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

cause.mod = ind~ age_n + gender_num + RECAA + RECDM +RecHCV + immuno_group+ DonorHCV + DON_TY + donorage + DONORAA + genderd_num +
  DONCIG + DONORCVA + COLDISCH + ORGANDIST + log_size

models = list(treat.mod = treat.mod, cens.mod = cens.mod, out.mod = out.mod,cause.mod  =cause.mod)

#Censored and composite outcome bootstrap
set.seed(2024)
n_boot = 1000
boot.results.cens = AFT.boot(n_boot, data, models, regime = "censor")
summary_boot_cens = summaryAFT(data, list(est = boot.results.cens$est, boot.est = boot.results.cens$boot.est, regime = "censor"))

boot.results.comp = AFT.boot(n_boot, data, models, regime = "composite")
summary_boot_comp = summaryAFT(data, list(est = boot.results.comp$est, boot.est = boot.results.comp$boot.est, regime = "composite"))

dat_cens = summary_boot_cens$data_comp
dat_comp = summary_boot_comp$data_comp

plots_cens = plotAFT(dat_cens, "(a) Cause-specific regime")
plots_comp = plotAFT(dat_comp, "(a) Composite regime")

cause_plot_comparator =ggarrange(plots_cens$p1, plots_comp$p1 ,plots_cens$p2, plots_comp$p2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
overall_plot_comparator = ggarrange(plots_cens$p3, plots_comp$p3, ncol = 2, nrow = 1, common.legend = T, legend = "right")

#Compute ITR metrics for all the regimes
cens_mod = glm(cens.mod, data = data, family = binomial)
w = 1/predict(cens_mod, newdata = data %>% filter(delta == 1), type = "response")

ITR_metrics(dat_w, dat_w$oracle, w)
ITR_metrics(dat_g, dat_w$oracle, w)
ITR_metrics(dat_cens, dat_w$oracle, w)
ITR_metrics(dat_comp, dat_w$oracle, w)