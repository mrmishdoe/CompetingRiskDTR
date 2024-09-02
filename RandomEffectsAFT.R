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
aft = function(data, treat.mod, cens.mod, outcome.mod, cause, corstr = "exchangeable"){
  
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
  model = geem(formula = outcome.mod, data = data, weights = w, id = group, corstr = corstr)
  
  #Get main effect and interactions with treatment (all named elements of vector containing a)
  coefs = coef(model)
  
  coefs[grepl("a", names(coefs), fixed = T)]
}

sim_data = function(n, psi1, psi2, cens = T, treat_clust = F, re_dist = "norm", binX = F, total_var = 0.5, ICC = 0.5, nclust = 10, cens_low = T, treat_re= "medium"){
  #Covariates
  x <- rnorm(n, 0, 1)
  
  #Option for binary effect modifier
  if(binX){
    x<- rbinom(n,1, 0.5)
    }
  
  z = rnorm(n, 0, 2)
  
  #Clustering in treatment and outcome (we assume that clustering variable is the same)
  group = sample(nclust, n, replace = T)
  
  #Decomposition of total variance = tau^2 + sigma^2 
  tau = sqrt(total_var*ICC)
  sigma = sqrt(total_var - tau^2)
    
  #Clustering in the treatment allocation
  if(treat_clust){
    #standard deviation of treatment random effect
    sd_treat= case_when(treat_re == "low" ~ 0.1, treat_re == "medium" ~ 0.5, treat_re == "high"~1)
    
    treat_effects = rnorm(nclust, mean = 0, sd = sd_treat)
    ind_treat_effects = map_dbl(group, function(x) treat_effects[x])
    a <- rbinom(n, 1, expit(0.5 + x +z + ind_treat_effects))
  } else{
    a <- rbinom(n, 1, expit(0.5 + x +z ))
  }
  
  #Generate cause indicator epsilon
  epsilon = rbinom(n,1,expit(x+0.5))+1
  
  
  #Generate groups and random intercepts for each group
  if(re_dist == "norm"){
    rand_effects = rnorm(nclust, mean = 0, sd = tau)
  } else if(re_dist == "gamma"){
    rand_effects = rgamma(nclust, shape = tau^2) - tau^2
  }
  
  ind_effects = map_dbl(group, function(x) rand_effects[x])
  
  #Parameters for generating outcomes
  beta1 <- c(1, 0.1, -0.3) #1, 0.5
  h1beta <- model.matrix(~x + z)
  h1psi <- model.matrix(~x)
  err1 = rnorm(sum(epsilon == 1),sd = sigma)
  
  beta2 <- c(2, -0.1,0.2) #2, -0.1
  h2beta <- model.matrix(~x + z)
  h2psi <- model.matrix(~x)
  err2 = rnorm(sum(epsilon == 2),sd = sigma)
  
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
    intercept = ifelse(cens_low, 1.73, 0)
    delta <- rbinom(n, 1, expit(intercept-x -0.3*z)) #1.73 for 20% censoring, 0 for 50% censoring
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
  cause_mod = glm(ind~ x, data = cause_data, family = binomial)
  
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
  
  #For now, remove the log (so estimating the value based on survival time instead of log survival time)
  mean_value_w = mean(test$T_w)
  mean_value_g = mean(test$T_g)
  
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



simAFT = function(n_rep, n_train, n_test, psi1, psi2, models, corstr = "exchangeable", ...){
  
  #Generate test set
  test = sim_data(n_test, psi1, psi2, cens = F, ...)
  
  #Initialize raw result list
  raw = list()
  for(name in names(models)){
    raw[[name]] = list(psi1_hat = NULL, psi2_hat = NULL,metrics = list(weighted = list(pot = NULL, value = NULL, blip = NULL), greedy = list(pot =NULL, value = NULL, blip =NULL)))
  }
  
  #Main loop
  for(i in 1:n_rep){
    
    #Generate training set for estimating blip parameters
    train = sim_data(n_train, psi1, psi2, ...)
    
    #Estimate parameters for every model in list 
    for (name in names(models)){
      model = models[[name]]
      psi1_hat = aft(train, model$treat, model$cens, model$out,1, corstr)
      psi2_hat = aft(train, model$treat, model$cens, model$out,2, corstr)
      res_DTR = eval_DTR(train, test, psi1, psi2, psi1_hat, psi2_hat)
      
      cur = raw[[name]]
      updated_metrics = list(weighted = update_metrics(cur$metrics, res_DTR, "weighted"), greedy = update_metrics(cur$metrics, res_DTR, "greedy"))
      raw[[name]] = list(psi1_hat = rbind(cur$psi1_hat, psi1_hat), psi2_hat = rbind(cur$psi2_hat, psi2_hat), metrics = updated_metrics)
    }
    
    
  }
  
  #Final results
  res = list()
  
  #Initialize boxplots
  p0 = ggplot() + xlab("Model Specification")
  p1 = p0 + geom_hline(yintercept = psi1[1], linetype = 2) + labs(y = expression(hat(psi)[11]))
  p2 = p0 + geom_hline(yintercept = psi1[2], linetype = 2) + labs(y = expression(hat(psi)[12]))
  p3 = p0 + geom_hline(yintercept = psi2[1], linetype = 2) + labs(y = expression(hat(psi)[21]))
  p4 = p0 + geom_hline(yintercept = psi2[2], linetype = 2) + labs(y = expression(hat(psi)[22]))
  
  
  #Combine results across all iterations
  for(name in names(models)){
    
    #Blip plots
    metrics = raw[[name]]$metrics
    res[[name]]$blip_plot = plot_results(test, metrics$weighted$blip, metrics$greedy$blip)
    
    #POT and value
    #Optimal value with known cause of failure (for comparison)
    test = test %>% mutate(T_opt = opt*T_1 + (1-opt)*T_0)
    
    #Value with random assignment of treatment
    test = test %>% mutate(opt_rand = rbinom(n_test, 1, 0.5), T_rand = opt_rand*T_1 + (1-opt_rand)*T_0)
      
    res[[name]]$measures = list(weighted = c(POT = mean(metrics$weighted$pot), Value = mean(metrics$weighted$value) ), 
                                greedy = c(POT = mean(metrics$greedy$pot), Value = mean(metrics$greedy$value) ), 
                                Opt_Value = mean(test$T_opt), Rand_Value= mean(test$T_rand))
    
    #Bias and SE for parameter estimates
    res[[name]]$inference = bias_SE(cbind(raw[[name]]$psi1_hat, raw[[name]]$psi2_hat), c(psi1, psi2))
      
    #Box plots for parameters
      
    #generate
    p1 = p1 + geom_boxplot(data = raw[[name]]$psi1_hat[,1] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p2 = p2 + geom_boxplot(data = raw[[name]]$psi1_hat[,2] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p3 = p3 + geom_boxplot(data = raw[[name]]$psi2_hat[,1] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p4 = p4 + geom_boxplot(data = raw[[name]]$psi2_hat[,2] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
      
      
  }
  
  #Add boxplots to final res
  res$boxplots = list(p1, p2, p3,p4)
  
  #Return result list
  res
}

#Function to print output of AFT simulations
print_res = function(res, name, save= F, boxplot = F){
  
  #Save raw results
  if(save){
    saveRDS(res, file = paste("raw_results/",name,".rds", sep = ""))
  }
  
  #For each model, print out all info
  for(i in 1:(length(res) - 1) ){
    mod_name = names(res)[i]
    mod_list = res[[mod_name]]
    
    print(mod_list$blip_plot)
    print(mod_list$measures)
    print(mod_list$inference)
    
    #Save blip plot
    if(save){
      ggsave(paste("plots/", name, "_" ,mod_name, ".png",  sep = ""), mod_list$blip_plot, width=12, height = 12)
    }
  }
  
  if(boxplot){
    #Print boxplots
    for(i in 1:length(res$boxplots)){
      print(res$boxplots[[i]])
      
      if(save){
        ggsave(paste("plots/", name, "_box", i, ".png",  sep = ""), res$boxplots[[i]], width=12, height = 12)
      }
    }
    
  }
  
}