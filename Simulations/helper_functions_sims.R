#One stage competing risks AFT model
library(lme4)
library(tidyverse)
library(ggpubr)
library(geeM)

expit = function(x) exp(x)/(1 + exp(x))

#aft function with GEE's
aft = function(data, treat.mod, cens.mod, outcome.mod, cause, corstr = "exchangeable", treat = "a"){
  
  #Treatment model and censoring model 
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")
  
  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data[[treat]] - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  w = weights[index]
  
  data = data %>% mutate(w = w)
  
  #Fitting GEE with exchangeable correlation structure, geeM package requires data.frame
  
  data = data %>% arrange(group) %>% as.data.frame()
  model = geem(formula = outcome.mod, data = data, weights = w, id = group, corstr = corstr, sandwich = F)
  
  #Get main effect and interactions with treatment (all named elements of vector containing name of treatment)
  
  coefs = coef(model)
  
  blip = coefs[grepl(treat, names(coefs), fixed = T)]
  
  #list(mod = model, blip = blip)
  list(blip = blip, mod = model)
}


bias_SE = function(psi_rep, psi){
  tab = matrix(0, nrow =2, ncol = length(psi))
  tab[1,] = apply(psi_rep, 2, mean) - psi
  tab[2,] = apply(psi_rep, 2, sd)
  
  tab
}

eval_DTR = function(test, psi1_hat, psi2_hat, cause_prob, regime){
  
  #For censored and composite outcomes, psi1_hat refers to main cause of interest
  #Number of test observations
  n_test = nrow(test)
  
  #4 different choices, weighted, greedy, censor, composite (censor and composite both require same evaluation)
  if(regime == "weighted"){
    test  = test %>% mutate(blip_hat = (cbind(1,x) %*% psi1_hat)*cause_prob + (cbind(1,x) %*% psi2_hat)*(1-cause_prob),
                            opt_hat = as.numeric(blip_hat > 0), T_hat =opt_hat*T_1 + (1-opt_hat)*T_0)
    
  }else if(regime == "greedy"){
    test = test %>% mutate(blip_hat = ifelse(cause_prob > 0.5,cbind(1,x) %*% psi1_hat,cbind(1,x) %*% psi2_hat),
                           opt_hat = as.numeric(blip_hat > 0),T_hat =opt_hat*T_1 + (1-opt_hat)*T_0)
    
  }else if(regime == "other"){
    test = test %>% mutate(blip_hat = cbind(1,x) %*% psi1_hat,
                           opt_hat = as.numeric(blip_hat > 0),T_hat =opt_hat*T_1 + (1-opt_hat)*T_0)
  }
  
  pot = sum(test$opt == test$opt_hat)/n_test
  mean_value = mean(log(test$T_hat))
  
  
  list(pot = pot, value = mean_value, blip = test$blip_hat)
}
  

#Helper function to plot graph for each cause
get_plot_data = function(test, blips){
  
  #Get 2.5%, 50% and 97.5% percentiles for blips of each individual in the test set
  blip_quantiles = t(apply(blips, 2, function(x) quantile(x, probs = c(0.025,0.50, 0.975)))) %>% as_tibble() 
  colnames(blip_quantiles) <- c("lower", "mid", "upper")
  
  dat = bind_cols(test, blip_quantiles) %>% arrange(blip) #oracle
  
  dat
  
}

sim_data = function(n, psi1, psi2, cens = T, treat_clust = F, re_dist = "norm", binX = F, total_var = 0.5, 
                    ICC = 0.5, nclust = 50, cens_low = T, treat_re= "medium", epsilon_0= 0, 
                    K_delta = 0){
  #Covariates
  x <- rnorm(n, 0, 1)
  
  #Option for binary effect modifier
  if(binX){
    x<- sample(c(-2,2), n, replace = T)
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
  epsilon = rbinom(n,1,expit(x+0.5 + epsilon_0 + K_delta*a))+1
  
  #Generate groups and random intercepts for each group
  if(re_dist == "norm"){
    rand_effects = rnorm(nclust, mean = 0, sd = tau)
  } else if(re_dist == "gamma"){
    rand_effects = rgamma(nclust, shape = tau^2) - tau^2
  }
  
  ind_effects = map_dbl(group, function(x) rand_effects[x])
  
  #Parameters for generating outcomes
  beta1 <- c(1, 0.5, -0.3) #1, 0.5
  h1beta <- model.matrix(~x + z)
  h1psi <- model.matrix(~x)
  err1 = rnorm(sum(epsilon == 1),sd = sigma)
  
  beta2 <- c(2, -0.1,0.2) #2, -0.1
  h2beta <- model.matrix(~x + z)
  h2psi <- model.matrix(~x)
  err2 = rnorm(sum(epsilon == 2),sd = sigma)
  
  #Generate counterfactuals for evaluating regimes
  T1_0 = exp(h1beta[epsilon== 1, ] %*% beta1 + ind_effects[epsilon ==1]+ err1)
  T1_1 = exp(h1beta[epsilon== 1, ] %*% beta1 + h1psi[epsilon ==1, ] %*% psi1 + ind_effects[epsilon ==1]+ err1)
  
  T2_0 = exp(h2beta[epsilon == 2, ] %*% beta2 + ind_effects[epsilon ==2] +err2)
  T2_1 = exp(h2beta[epsilon == 2, ] %*% beta2 + h2psi[epsilon ==2, ] %*% psi2 + ind_effects[epsilon ==2] +err2)
  
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
  
  #Optimal value with known cause of failure (for comparison)
  data = data %>% mutate(T_opt = opt*T_1 + (1-opt)*T_0)
  
  #Value with random assignment of treatment
  data = data %>% mutate(opt_rand = rbinom(nrow(data), 1, 0.5), T_rand = opt_rand*T_1 + (1-opt_rand)*T_0)
  
  #Replace times by censoring time for those who were censored
  if(cens){  
    intercept = ifelse(cens_low, 1.73, 0)
    delta <- rbinom(n, 1, expit(intercept-x -0.3*z)) #1.73 for 20% censoring, 0 for 50% censoring
    C <- rexp(n, 1/300)
    data = data %>% mutate(Y = ifelse(delta,Y,C), delta = delta)
  }
  
  data
}

#Should turn this into two different functions to decouple estimation and plotting
simAFT = function(n_rep, n_train, n_test, psi1, psi2, models, corstr = "exchangeable",save = F, file = "default", 
                  cause_mod = T,...){
  
  #Generate test dataset
  test = sim_data(n_test, psi1, psi2, cens = F, ...)
  
  #### Replicates ####
  results = list()
  
  for(name in names(models)){
    
    #Parallel processing
    cl = makeCluster(5)
    clusterEvalQ(cl,{
      library(tidyverse)
      source("helper_functions_sims.R")
    })
    
    rep.est = function(i, models, corstr, n_train, ...){
      
      #Generate training set for estimating blip parameters
      train = sim_data(n_train, psi1, psi2, ...)%>% mutate(ind = as.numeric(epsilon ==1))
      
      psi1 = aft(train, models$treat, models$cens, models$out,1, corstr)$blip
      print("here")
      psi2 = aft(train, models$treat, models$cens, models$out,2, corstr)$blip
      
      if(cause_mod){
        cause_mod = glm(ind~x, data = train %>% filter(delta == 1) , family = binomial)
      }else{
        cause_mod = glm(ind~1, data = train %>% filter(delta == 1) , family = binomial)
      }
      
      cause_prob = predict(cause_mod,newdata = test, type = "response")
      
      #Estimates for censored and composite outcomes
      
      psi_censor = aft(train %>% mutate(delta = ifelse(delta == 1 & epsilon == 1, 1, 0)), 
                       models$treat, models$cens, models$out,1, corstr)$blip
      
      psi_composite = aft(train %>% mutate(epsilon = 1), 
                          models$treat, models$cens, models$out,1, corstr)$blip
      
      list(psi1_hat = psi1, psi2_hat = psi2, psi_ce =psi_censor, psi_co = psi_composite, cause_prob = cause_prob)  

    }
    
    
    rep_res = parLapply(cl,1:n_rep, rep.est, models = models[[name]], corstr = corstr, n_train = n_train, ...)
    stopCluster(cl) 
    
    results[[name]]= rep_res
  }
  
  #Save raw results
  if(save){
    saveRDS(results, file = paste("Results/",file,".rds", sep=""))
  }
  
  list(results = results, test = test)
  
}


summaryAFT= function(data, results){
  
  #Final results
  res = list()
  
  #Initialize boxplots
  p0 = ggplot() + xlab("Model Specification")
  p1 = p0 + geom_hline(yintercept = psi1[1], linetype = 2) + labs(y = expression(hat(psi)[11]))
  p2 = p0 + geom_hline(yintercept = psi1[2], linetype = 2) + labs(y = expression(hat(psi)[12]))
  p3 = p0 + geom_hline(yintercept = psi2[1], linetype = 2) + labs(y = expression(hat(psi)[21]))
  p4 = p0 + geom_hline(yintercept = psi2[2], linetype = 2) + labs(y = expression(hat(psi)[22]))
  
  #Compile raw results
  for(name in names(models)){
    
    raw = results[[name]]
    n_boot = length(raw)
    n = nrow(data)
    
    psi1_res = matrix(0, nrow = n_boot, ncol = 2)
    psi2_res = matrix(0, nrow = n_boot, ncol = 2)
    psice_res = matrix(0, nrow = n_boot, ncol = 2)
    psico_res = matrix(0, nrow = n_boot, ncol = 2)
    
    weighted_blip= matrix(0, nrow = n_boot, ncol = n)
    greedy_blip = matrix(0, nrow = n_boot, ncol = n)
    censor_blip = matrix(0, nrow = n_boot, ncol = n)
    comp_blip = matrix(0, nrow = n_boot, ncol = n)
    
    pot_w = rep(0, n_boot)
    pot_g = rep(0, n_boot)
    pot_ce = rep(0, n_boot)
    pot_co = rep(0, n_boot)
    
    value_w = rep(0, n_boot)
    value_g = rep(0, n_boot)
    value_ce = rep(0, n_boot)
    value_co = rep(0, n_boot)
    
    
    for(i in 1:n_boot){
      psi1_res[i,] = raw[[i]]$psi1_hat
      psi2_res[i,] = raw[[i]]$psi2_hat
      psice_res[i,] = raw[[i]]$psi_ce
      psico_res[i,] = raw[[i]]$psi_co
      
      eval_w = eval_DTR(data, raw[[i]]$psi1_hat, raw[[i]]$psi2_hat, raw[[i]]$cause_prob, "weighted")
      eval_g = eval_DTR(data, raw[[i]]$psi1_hat, raw[[i]]$psi2_hat, raw[[i]]$cause_prob, "greedy")
      
      eval_ce = eval_DTR(data, raw[[i]]$psi_ce, NULL, NULL, "other")
      eval_co = eval_DTR(data, raw[[i]]$psi_co, NULL, NULL, "other")
      
      weighted_blip[i,]= eval_w$blip
      greedy_blip[i,] = eval_g$blip
      censor_blip[i,] = eval_ce$blip
      comp_blip[i,] = eval_co$blip
      
      pot_w[i] = eval_w$pot
      pot_g[i] = eval_g$pot
      pot_ce[i] = eval_ce$pot
      pot_co[i] = eval_co$pot
      
      value_w = eval_w$value
      value_g = eval_g$value
      value_ce = eval_ce$value
      value_co = eval_co$value
    }
    
    #Combine results
    res_temp = list()
    
    #Blip plots
    res_temp$data_w = get_plot_data(data, weighted_blip)
    res_temp$data_g = get_plot_data(data, greedy_blip)
    res_temp$data_ce = get_plot_data(data, censor_blip)
    res_temp$data_co = get_plot_data(data, comp_blip)
    
    
    #POT and value
    res_temp$measures =  list(weighted = c(POT = mean(pot_w), Value = mean(value_w)), 
                         greedy = c(POT = mean(pot_g), Value = mean(value_g) ), 
                         censor = c(POT = mean(pot_ce), Value = mean(value_ce) ),
                         comp = c(POT = mean(pot_co), Value = mean(value_co) ),
                         Opt_Value = mean(log(data$T_opt)), Rand_Value= mean(log(data$T_rand)))
    
    
    res_temp$inference = bias_SE(cbind(psi1_res, psi2_res), c(psi1, psi2))
    
    #Fill boxplots
    p1 = p1 + geom_boxplot(data = psi1_res[,1] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p2 = p2 + geom_boxplot(data = psi1_res[,2] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p3 = p3 + geom_boxplot(data = psi2_res[,1] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    p4 = p4 + geom_boxplot(data = psi2_res[,2] %>% as_tibble() %>% mutate(scenario = which(names(models) == name)), aes(x = scenario,y = value))
    
    
    res[[name]] = res_temp
    
  }
  
  res$boxplots = list(p1=p1, p2=p2, p3=p3, p4=p4)
  #Return result list
  res
}

plotAFT = function(dat, name){
  
  dat1 = dat %>% filter(epsilon == 1)
  dat2 = dat %>% filter(epsilon == 2)
  
  colors <- c("Chosen Regime" = "black", "Oracle Regime" = "red")  
  
  p1 =ggplot(dat1, aes(x = 1:nrow(dat1))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid, color = "Chosen Regime"), linewidth= 1.2) + geom_line(aes(y = blip, color = "Oracle Regime"), linewidth = 1.2)+labs(x="Observation number (Ordered)", y = "Median estimated benefit",title ="Cause 1", color = "Legend",caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  p2 =ggplot(dat2, aes(x = 1:nrow(dat2))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid,color = "Chosen Regime"), linewidth= 1.2) + geom_line(aes(y = blip, color = "Oracle Regime"), linewidth = 1.2)+labs(x="Observation number (Ordered)", y = "Median estimated benefit",title ="Cause 2",color = "Legend", caption = name) + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  list(p1 = p1, p2 = p2)
}