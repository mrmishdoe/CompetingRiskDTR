########## Setup
library(tidyverse)
library(ggpubr)
library(parallel)

source("helper_functions_sims.R")

########## Scenario 1 with high POT, high value for both regimes
psi1 = c(0.2, -0.2)
psi2 = c(0.2, 0.2)

models = list(both = list(treat = a~x, cens = delta~x, out = log(Y)~ x+ a+a:x), treat = list(treat = a~x, cens = delta~x, out = log(Y)~ x+z+ a+a:x), 
              out = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+ a+a:x), none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_b = simAFT(n_rep, n_train, n_test, psi1, psi2, models)
sum_b = summaryAFT(raw_b$test, raw_b$results)


########## Scenario 1 but with larger training sample size
psi1 = c(0.2, -0.2)
psi2 = c(0.2, 0.2)

models = list(both = list(treat = a~x, cens = delta~x, out = log(Y)~ x+ a+a:x), treat = list(treat = a~x, cens = delta~x, out = log(Y)~ x+z+ a+a:x), 
              out = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+ a+a:x), none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

#5x sample size so approx half the SE
n_rep = 1000
n_train = 5000
n_test = 10000

set.seed(2024)
raw_b2 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, nclust = 250)
sum_b2 = summaryAFT(raw_b2$test, raw_b2$results)


########## Scenario 2, weighted performs best
#Opposite strategies for both causes, always treat for one cause, never treat for the other
#Smaller POT but much higher value

psi1 = c(3, -0.5)
psi2 = c(-1, 0.2)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2023)
raw_s1 = simAFT(n_rep, n_train, n_test, psi1, psi2, models)
sum_s1 = summaryAFT(raw_s1$test, raw_s1$results)

########## Scenario 3: Greedy performs the best
#When one cause is more likely and the less likely cause does not have very negative blips
#Much higher POT but value a bit smaller

psi1= c(-0.5, -0.7)
psi2= c(0.1, 0.08)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s2 = simAFT(n_rep, n_train, n_test, psi1, psi2, models)
sum_s2 = summaryAFT(raw_s2$test, raw_s2$results)


########## Scenario 4: Same slopes and similar performance for both regimes
psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s3 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1)
sum_s3 = summaryAFT(raw_s3$test, raw_s3$results)


########## Scenario 5: Different ICC
psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s4_1 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1 , ICC = 0.9)
sum_s4_1 = summaryAFT(raw_s4_1$test, raw_s4_1$results)

set.seed(2024)
raw_s4_2 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1 , ICC = 0.1)
sum_s4_2 = summaryAFT(raw_s4_2$test, raw_s4_2$results)

########## Scenario 6: Different error distribution
#Centered gamma distribution for error term

psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s5 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, re_dist = "gamma")
sum_s5 = summaryAFT(raw_s5$test, raw_s5$results)


########## Scenario 7: 50% censoring instead of 20%
psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s6 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, cens_low = F)
sum_s6 = summaryAFT(raw_s6$test, raw_s6$results)


########## Scenario 8: ignoring clustering i.e. use independence structure in GEE over exchangeable
psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s7 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, corstr = "independence")
sum_s7 = summaryAFT(raw_s7$test, raw_s7$results)


########## Scenario 9: ignoring treatment clustering
#3 different levels of clustering

psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s8_1 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, treat_clust = T)
sum_s8_1 = summaryAFT(raw_s8_1$test, raw_s8_1$results)

#Large random effect
set.seed(2024)
raw_s8_2 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, treat_clust = T, treat_re = "high")
sum_s8_2 = summaryAFT(raw_s8_2$test, raw_s8_2$results)

#Small random effect
set.seed(2024)
raw_s8_3 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, total_var = 1, treat_clust = T, treat_re = "low")
sum_s8_3 = summaryAFT(raw_s8_3$test, raw_s8_3$results)

########## Scenario 10: comparison with cause-specific regime
psi1 = c(1, -0.5)
psi2 = c(-3, 0.2)
models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))
#
n_rep = 1000
n_train = 1000
n_test = 10000
#
set.seed(2023)
raw_s10 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, epsilon_0=2)
sum_s10 = summaryAFT(raw_s10$test, raw_s10$results)

#Example of plotting
plots_w = plotAFT(sum_s10$none$data_w, "(a) Weighted regime")
plots_g = plotAFT(sum_s10$none$data_g, "(b) Greedy regime")
plots_ce = plotAFT(sum_s10$none$data_ce, "Cause-specific regime")

cause_plot=ggarrange(plots_w$p1, plots_g$p1 ,plots_w$p2, plots_g$p2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
ggsave("plots/s10_none.png", cause_plot,width=12, height = 12)
plot_ce = ggarrange(plots_ce$p1, plots_ce$p2, ncol = 2, common.legend = T, legend = "right")
ggsave("plots/s10_cs.png", plot_ce,width=12, height = 6)

########## Scenario 11: comparison with composite outcome regime
psi1= c(-5, -5)
psi2= c(5, 5)

models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+a+a:x))#z+ 

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s11 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, epsilon_0 = -1)
sum_s11 = summaryAFT(raw_s11$test, raw_s11$results)

#Plotting
plots_w = plotAFT(sum_s11$none$data_w, "(a) Weighted regime")
plots_g = plotAFT(sum_s11$none$data_g, "(b) Greedy regime")
plots_co = plotAFT(sum_s11$none$data_co, "Composite outcome regime")

plot_co = ggarrange(plots_co$p1, plots_co$p2, ncol = 2, common.legend = T, legend = "right")
ggsave("plots/s11_comp.png", plot_co,width=12, height = 6)
plot_wg =ggarrange(plots_w$p1, plots_g$p1 ,plots_w$p2, plots_g$p2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
ggsave("plots/s11_none.png", plot_wg ,width=12, height = 12)

########## Scenario 12: misspecified cause model
psi1= c(0.6,-0.6)
psi2= c(-0.6, -0.6)
models = list(none = list(treat = a~x+z, cens = delta~x+z, out = log(Y)~ x+z+ a+a:x))

n_rep = 1000
n_train = 1000
n_test = 10000

set.seed(2024)
raw_s12 = simAFT(n_rep, n_train, n_test, psi1, psi2, models, cause_mod = F)
sum_s12 = summaryAFT(raw_s12$test, raw_s12$results)
