#One stage competing risks AFT model
library(lme4)
library(tidyverse)
library(ggpubr)
library(geeM)

expit = function(x) exp(x)/(1 + exp(x))

#aft function with linear mixed model approach
aft_lmm = function(data, treat.mod, cens.mod, outcome.mod, cause){
  
  #Treatment model and censoring model 
  
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")

  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data$a - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  weights = weights[index]
  
  #Fitting linear mixed model
  
  model = lmer(formula = outcome.mod, data = data, weights = weights)
  fixef(model)
  
}

#aft function with GEE's
aft = function(data, treat.mod, cens.mod, outcome.mod, cause){
  
  #Treatment model and censoring model 
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")
  
  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data$a - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  w = weights[index]
  
  data = data %>% mutate(w = w)
  
  #Fitting GEE with exchangeable correlation structure, geeM package requires data.frame
  
  data = data %>% arrange(group) %>% as.data.frame()
  model = geem(formula = outcome.mod, data = data, weights = w, id = group, corstr = "exchangeable")
  
  #Get main effect and interactions with treatment (all named elements of vector containing a)
  coefs = coef(model)
  
  coefs[grepl("a", names(coefs), fixed = T)]
}

sim_data = function(n, psi1, psi2, cens = T, treat_clust = F, re_dist = "norm"){
  #Covariates
  #x <- rnorm(n, 0, 1) #was = 2 initially
  x<- rbinom(n,1, 0.5) #can put this as an option
  z = rnorm(n, 0, 2)
  
  #Clustering in treatment and outcome (we assume that clustering variable is the same)
  nclust = 10
  group = sample(nclust, n, replace = T)
  
  #Clustering in the treatment allocation
  if(treat_clust){
    treat_effects = rnorm(nclust, mean = 0, sd = 0.5)
    ind_treat_effects = map_dbl(group, function(x) treat_effects[x])
    a <- rbinom(n, 1, expit(0.5 + x +z + ind_treat_effects))
  } else{
    a <- rbinom(n, 1, expit(0.5 + x +z ))
  }
  
  #Generate cause indicator epsilon
  epsilon = rbinom(n,1,expit(x+0.5))+1
  
  
  #Generate groups and random intercepts for each group
  if(re_dist == "norm"){
    rand_effects = rnorm(nclust, mean = 0, sd = 0.5)
  } else if(re_dist == "gamma"){
    rand_effects = rgamma(nclust, shape = 0.25) - 0.25
  }
  
  ind_effects = map_dbl(group, function(x) rand_effects[x])
  
  #Parameters for generating outcomes
  beta1 <- c(1, 0.5, -0.3)
  h1beta <- model.matrix(~x + z)
  h1psi <- model.matrix(~x)
  err1 = rnorm(sum(epsilon == 1),sd = 0.5)
  
  beta2 <- c(2, -0.1,0.2)
  h2beta <- model.matrix(~x + z)
  h2psi <- model.matrix(~x)
  err2 = rnorm(sum(epsilon == 2),sd = 0.5)
  
  #Generate counterfactuals for evaluating regimes
  T1_0 = exp(h1beta[epsilon== 1, ] %*% beta1 + ind_effects[epsilon ==1]+ err1)
  T1_1 = exp(h1beta[epsilon== 1, ] %*% beta1+ h1psi[epsilon ==1, ] %*% psi1 + ind_effects[epsilon ==1]+ err1)
  
  T2_0 = exp(h2beta[epsilon == 2, ] %*% beta2+ ind_effects[epsilon ==2] +err2)
  T2_1 = exp(h2beta[epsilon == 2, ] %*% beta2+ h2psi[epsilon ==2, ] %*% psi2 + ind_effects[epsilon ==2] +err2)
  
  T_0 = T_1 = rep(0,n)
  T_0[epsilon == 1] = T1_0
  T_0[epsilon == 2] = T2_0
  T_1[epsilon == 1] = T1_1
  T_1[epsilon == 2] = T2_1
  
  #Construct dataframe
  group = as.factor(group)
  data = tibble(x,z,a,group,epsilon, T_0, T_1) %>% mutate(Y = ifelse(a, T_1, T_0))
  
  #Create true blip associated to actual cause of failure
  data = data %>% mutate(blip = ifelse(epsilon ==1, cbind(1,x) %*% psi1, cbind(1,x) %*% psi2),opt = as.numeric(blip > 0))
  
  #Replace times by censoring time for those who were censored
  if(cens){  
    delta <- rbinom(n, 1, expit(-x +1.9-0.3*z))
    C <- rexp(n, 1/300)
    data = data %>% mutate(Y = ifelse(delta,Y,C), delta = delta)
  }
  
  data
}

bias_SE = function(psi_rep, psi){
  tab = matrix(0, nrow =2, ncol = length(psi))
  tab[1,] = apply(psi_rep, 2, mean) - psi
  tab[2,] = apply(psi_rep, 2, sd)
  
  tab
}

eval_DTR = function(train, test, psi1, psi2, psi1_hat, psi2_hat){
  
  #Model for probability of failing from cause 1
  cause_data = train %>% filter(delta == 1) %>% mutate(ind = as.numeric(epsilon ==1))
  cause_mod = glm(ind~ x+z, data = cause_data, family = binomial)
  
  #Get the optimal regime depending on true underlying cause
  #Note : have to make sure this works the way i think it does
  n_test = nrow(test)
  
  #Allocate treatment according to weighted rule
  test  = test %>% mutate(cause_prob = predict(cause_mod,newdata = test, type = "response"), blip_weighted = (cbind(1,x) %*% psi1_hat)*cause_prob + (cbind(1,x) %*% psi2_hat)*(1-cause_prob),
                          opt_weighted = as.numeric(blip_weighted > 0))
  
  #Allocate treatment according to greedy rule
  
  test = test %>% mutate(blip_greedy = ifelse(cause_prob > 0.5,cbind(1,x) %*% psi1_hat,cbind(1,x) %*% psi2_hat), opt_greedy = as.numeric(blip_greedy > 0))
  
  #Identify proportion of optimal treatment (POT)
  pot_w = sum(test$opt == test$opt_weighted)/n_test
  pot_g = sum(test$opt == test$opt_greedy)/n_test
  
  #Compute mean log survival time under both regimes given true data generating values of beta/psi
  test = test %>% mutate(T_w =opt_weighted*T_1 + (1-opt_weighted)*T_0, T_g = opt_greedy*T_1 + (1-opt_greedy)*T_0)
  mean_value_w = mean(log(test$T_w))
  mean_value_g = mean(log(test$T_g))
  
  list(weighted = list(pot = pot_w, value = mean_value_w, blip = test$blip_weighted), greedy = list(pot = pot_g, value = mean_value_g, blip = test$blip_greedy))
}

#Helper function to update metrics for DTRs

update_metrics = function(old_metrics, new_metrics, regimen){
  old = old_metrics[[regimen]]
  new = new_metrics[[regimen]]
  
  list(pot = rbind(old$pot, new$pot), value = rbind(old$value, new$value), blip = rbind(old$blip, t(new$blip)))
}

#Helper function to plot graph for each cause

plot_cause = function(test, blips, name){
  
  #Get 2.5%, 50% and 97.5% percentiles for blips of each individual in the test set
  blip_quantiles = t(apply(blips, 2, function(x) quantile(x, probs = c(0.025,0.5, 0.975)))) %>% as_tibble() 
  colnames(blip_quantiles) <- c("lower", "mid", "upper")
  
  dat = bind_cols(test, blip_quantiles) %>% arrange(blip)
  
  dat1 = dat %>% filter(epsilon == 1)
  dat2 = dat %>% filter(epsilon == 2)
  # Plotting
  
  colors <- c("Chosen Strategy" = "black", "True Blip" = "red")
  
  p1 =ggplot(dat1, aes(x = 1:nrow(dat1))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid, color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+ labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time",title ="Cause 1", color = "Legend",caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  p2 =ggplot(dat2, aes(x = 1:nrow(dat2))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid,color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+ labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time",title ="Cause 2",color = "Legend", caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  list(p1 = p1, p2 = p2)
}


plot_results = function(test, weighted_blips, greedy_blips){
  
  #Plot results for weighted and greedy strategies
  
  w_plots = plot_cause(test, weighted_blips, "(a) Weighted strategy")
  g_plots = plot_cause(test, greedy_blips, "(b) Greedy strategy")
  
  #Arrange plots
  blip_plot=ggarrange(w_plots$p1, g_plots$p1 ,w_plots$p2, g_plots$p2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
  
  blip_plot
}