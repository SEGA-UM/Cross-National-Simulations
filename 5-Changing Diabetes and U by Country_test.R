library(tidyverse)
library(parallel)
library(dplyr)
library(simsurv)
library(flexsurv)
library(survival)
library(data.table)
library(boot)
library(ggplot2)

#### 1. Import mortality statistics ####
# df_Mort_period contains the probability of death by sex and age in each country
df_Mort_period = read.csv("/Users/kelvin/Library/CloudStorage/GoogleDrive-kelvinzh@umich.edu/My Drive/SEGA/Projects/Cross-national Simulation/Data/periodlifetables_cleaned_20250605.csv")

##### 2. Input of incidence rate of diabetes and prevalence of dementia by country ####
# Age-standardized diabetes incidence rates are from Global Burden of Disease (GBD) in 2016
# Dementia prevalence are from harmonized estimates in HCAP
Info = function(C){
  DiabetesI = switch(C,
                     "US" = 0.00399,
                     "India" = 0.00236,
                     "England" = 0.00255,
                     "China" = 0.00230,
                     "Mexico" = 0.00425,
                     NA)
  
# Baseline hazard of diabetes lambda0 is calibrated so that the incidence rate over the follow-up period is close to age-standardized incidence rate in real world studies
  baseline = switch(C,
                    "US" = 0.000285,
                    "India" = 0.000159,
                    "England" = 0.000173,
                    "China" = 0.000155,
                    "Mexico" = 0.000307,
                    NA)
  
  df_Mort_subset = df_Mort_period %>% 
    filter(Country == C & Age >= 60)
  
  df_Mort_subset = df_Mort_subset %>% 
    dplyr::rename(
      Age.2016 = Age
    )
  
  df_Mort_subset$nq20_death = 1 - df_Mort_subset$np20_period
  
  DementiaP = switch(C,
                     "US" = 0.115, 
                     "India" = 0.088,
                     "England" = 0.139,
                     "China" = 0.051,
                     "Mexico" = 0.087,
                     NA)
  
  # We assume that countries with generally lower life expectancy will have larger effect sizes of diabetes and U on mortality
  g_Diabetes.Duration = switch(C,
               "England" = log(1.01), 
               "US" = log(1.02),
               "Mexico" = log(1.03),
               "China" = log(1.04),
               "India" = log(1.05),
               NA)
  
  g_U = switch(C,
               "England" = log(1.2), 
               "US" = log(1.4),
               "Mexico" = log(1.6),
               "China" = log(1.8),
               "India" = log(2),
               NA)
  
  return(list(DiabetesI = DiabetesI,
              baseline = baseline,
              df_Mort_subset = df_Mort_subset,
              DementiaP = DementiaP,
              g_Diabetes.Duration = g_Diabetes.Duration,
              g_U = g_U))
}

#### 3. Set parameters and causal structures ####
# Scenario 1: no collider stratification; Scenario 2: collider stratification; Scenario 3: collider stratification with interaction
paramGen = function(C,I){
  simInputs = data.frame(Scen = c(1,2,3), 
                         # Cumulative Mortality Parameters
                         # Effect of diabetes duration on mortality
                         g_Diabetes.Duration = c(Info(C)$g_Diabetes.Duration,Info(C)$g_Diabetes.Duration,Info(C)$g_Diabetes.Duration),       
                         
                         # Effect of U on mortality
                         g_U = c(0,Info(C)$g_U,Info(C)$g_U),   
                         
                         # Effect of U and diabetes duration interaction on mortality
                         g_Diab.durationU = c(0,0,log(1.02)),  
             
                         # Parameters of odds of dementia
                         # Effect of prevalent diabetes on dementia
                         b_Diabetes = log(2),  
                         
                         # Effect of U on dementia
                         b_U = log(1.2),    
                         
                         # Effect of age group on dementia
                         b_Age = log(1.2)      
  )
  
  return(list(Scenario = simInputs$Scen[I],
              gamma.param = simInputs[I,c("g_Diabetes.Duration", "g_U","g_Diab.durationU")],
              beta.param = simInputs[I,c("b_Diabetes", "b_U","b_Age")]))
}

#### 4. Generate output for one run ####
output = function(C,I){
  Info_by_country = Info(C)
  Inputs = paramGen(C,I)
  
  # For each birth year, we sample 30,000 participants. The birth years are ranged from 1927 to 1951
  df = data.frame(ID = 1:750000)
  df$Country = C
  df$Birthyear = sample(rep(1927:1951, each = 30000))

  # Age in 2016 is the duration between birthyear and the year of 2016
  df$Age.2016 = 2016 - df$Birthyear
  
  # The sample is divided into five age groups: 65-69, 70-74, 75-79, 80-84, 85-89
  df$Age.2016_cat = cut(df$Age.2016,
                      breaks = c(65, 70, 75, 80, 85, 90),
                      labels = c("1", "2", "3","4","5"),
                      right = FALSE)
  df$Age.2016_cat = as.numeric(df$Age.2016_cat)

  # Half of the sample is female
  df$Female = rbinom(n = 750000, size = 1, prob = 0.5)

  # U is a normally distributed variable
  df$U = rnorm(n = 750000)

  # Set the hazard function of incident diabetes
  simsurv.param = data.frame(
    baseline = rep(Info_by_country$baseline, nrow(df)),
    growth = rep(0.08, nrow(df)),
    peak_age = rep(65, nrow(df)),
    decline = rep(0.03, nrow(df))
    )

  cov = data.frame(
    id = df$ID)

  haz_fun = function(t,x,betas){
    current_age <- 20 + t
  
  hazard = ifelse(
    current_age < betas[["peak_age"]],
    betas[["baseline"]] * exp(betas[["growth"]] * (current_age - 20)),
    betas[["baseline"]] * exp(betas[["growth"]] * (betas[["peak_age"]] - 20) - betas[["decline"]] * (current_age - betas[["peak_age"]]))
    )
  
  return(hazard)
  }

  # Use simsurv() to simulate the time-to-event during the follow-up period
  simdat = simsurv(
    hazard = haz_fun,
    x = cov,
    betas = simsurv.param,
    maxt = 69
    )

  simdat = simdat %>% 
    rename(ID = id)

  # We right censored in 2016, so we only need time to event here.
  simdat_subset = subset(simdat, select = c("ID", "eventtime"))
  df_new = merge(simdat_subset, df, by = "ID")

  # Folllow-up time is the duration between age 20 to the year of 2016.
  df_new$FollowUptime = 2016 - (df_new$Birthyear + 20)

  # If the time to event is equal or less than your censor time, then Diabetes = 1
  df_new$Diabetes = ifelse(df_new$eventtime<=df_new$FollowUptime,1,0)

  # Time to event (diabetes/censor) is either the generated time to event or follow-up time, whichever is shorter
  df_new$Time.to.Event = pmin(df_new$eventtime, df_new$FollowUptime)

  # Diabetes duration is the length of time between censoring and time to event
  df_new$Diabetes.duration = df_new$FollowUptime - df_new$Time.to.Event

  df_new = df_new %>% 
   left_join(Info_by_country$df_Mort_subset %>% select(Age.2016, Female, nq20_death),
             by = c("Age.2016","Female"))

  # For each age and sex combination, we calibrate the intercept (g0) in the formula of death probability so that the overall probability of death matches the mortality statistics
  Age.sex_groups = unique(df_new[,c("Female", "Age.2016","nq20_death")])

  g0_results = data.frame(
    Female = numeric(),
    Age.2016 = numeric(),
    g0 = numeric()
    )

  for (i in seq_len(nrow(Age.sex_groups))){
   sub_df = subset(df_new,
                   Female == Age.sex_groups$Female[i] & Age.2016 == Age.sex_groups$Age.2016[i])
  
   find_intercept_surv <- function(g0) {
     p_surv_formula <- g0 + Inputs$gamma.param$g_Diabetes.Duration*sub_df$Diabetes.duration + Inputs$gamma.param$g_U*sub_df$U + Inputs$gamma.param$g_Diab.durationU*sub_df$Diabetes.duration*sub_df$U
    
     p_surv_est_individual <- exp(p_surv_formula) / (1 + exp(p_surv_formula))
     p_surv_known  <- Age.sex_groups$nq20_death[i]
     p_surv_est_pop <- mean(p_surv_est_individual)
    
     (p_surv_known - p_surv_est_pop)^2
     }
  
  opt_result <- optim(
    par = 0, 
    fn = find_intercept_surv,
    method = "BFGS"
    )
  
  g0_subgroup <- opt_result$par
  
  g0_results = rbind(
    g0_results,
    data.frame(
      Female = Age.sex_groups$Female[i],
      Age.2016 = Age.sex_groups$Age.2016[i],
      g0 = g0_subgroup
      )
  )
  }

  # Generate probability of survival as a function of diabetes duration and U
  df_new = df_new %>% 
   left_join(g0_results, by = c("Female", "Age.2016"))

  lin.pred = with(df_new, exp(df_new$g0 + Inputs$gamma.param$g_Diabetes.Duration*df_new$Diabetes.duration + Inputs$gamma.param$g_U*df_new$U + Inputs$gamma.param$g_Diab.durationU*df_new$Diabetes.duration*df_new$U))

  df_new$prob_death_est = with(df_new, lin.pred/ (1+lin.pred))

  df_new$Survival = rbinom(nrow(df_new), 1L, 1 - df_new$prob_death_est)

  # Test if the generated marginal probability of death in each subgroup is close to the mortality statistics
  df_surv_aggregated = df_new %>% 
    group_by(Female, Age.2016) %>% 
    summarise(
      obs_death = 1 - mean(Survival),
      pred_death = mean(nq20_death),
      perc_surv.diff = (obs_death - pred_death)*100,
      .groups = "drop"
      ) %>% 
    select(Female, Age.2016, perc_surv.diff)

  df_surv_aggregated = df_surv_aggregated %>% 
   distinct(Female, Age.2016, perc_surv.diff) 

  # Optimize the intercept to balance the marginal prevalence of dementia
  find_intercept_dementia = function(b0) {
   p_dementia_formula = b0 + Inputs$beta.param$b_Diabetes*df_new$Diabetes + Inputs$beta.param$b_U*df_new$U + Inputs$beta.param$b_Age*df_new$Age.2016_cat
  
   p_dementia_est_individual = exp(p_dementia_formula)/(1+exp(p_dementia_formula))
   p_dementia_est_overall = mean(p_dementia_est_individual)
   p_dementia_known = Info_by_country$DementiaP
  
   p_dementia_known - p_dementia_est_overall
   }

  b0 = uniroot(find_intercept_dementia, interval = c(-20,20))$root

  lin.pred = with(df_new,exp(b0 + Inputs$beta.param$b_Diabetes*Diabetes + Inputs$beta.param$b_U*U + Inputs$beta.param$b_Age*Age.2016_cat)) 

  df_new$p_dementia = with(df_new, lin.pred/ (1+lin.pred))

  df_new$Dementia  = rbinom(nrow(df_new), 1L, df_new$p_dementia)

  # Filter out those who survive in 2016
  df_survived = subset(df_new, Survival == 1)

  # Stratified sample 75000 individuals proportionate to the age distribution in the survived sample
  frac_need = 75000/nrow(df_survived)
  df_sample <- df_survived %>% 
    group_by(Age.2016_cat) %>% 
    sample_frac(size = frac_need) %>% 
    ungroup()

  # Estimated effect of diabetes duration on dementia in the survived population
  Model_pop = glm(Dementia ~ Diabetes + as.factor(Age.2016_cat), data = df_survived, family = binomial(link = "logit"))
  coef = coef(Model_pop)
  ci = confint.default(Model_pop)
  Model.sum_pop = data.frame((cbind(coef, ci)[2:6,]))
  colnames(Model.sum_pop) = c("Beta", "Lower_CI", "Upper_CI")
  
  # Estimated effect of diabetes duration on dementia in the sample
  Model_sample = glm(Dementia ~ Diabetes + as.factor(Age.2016_cat), data = df_sample, family = binomial(link = "logit"))
  coef = coef(Model_sample)
  ci = confint.default(Model_sample)
  Model.sum_sample = data.frame((cbind(coef, ci)[2:6,]))
  colnames(Model.sum_sample) = c("Beta", "Lower_CI", "Upper_CI")

  # Standardize estimate coefficient
  df_sample$Age.2016_cat = factor(df_sample$Age.2016_cat)

  # Standardized effect of diabetes on dementia among survived samples
  # Weights are from participants aged 65 to 89 in HRS-HCAP+ELSA-HCAP+LASI-DAD+Mex-Cog+CHARLS-HCAP
  age_levels = levels(df_sample$Age.2016_cat)
  std_wts = setNames(c(0.3784, 0.2559, 0.1897, 0.1217, 0.0543),
                     age_levels)

  grid = expand.grid(Diabetes = c(0, 1),
                     Age.2016_cat = age_levels,
                     KEEP.OUT.ATTRS = FALSE)

  std_or_fun = function(data, idx) {
    d = data[idx, ]                                   
    fit = glm(Dementia ~ Diabetes + factor(Age.2016_cat),
              data = d, family = binomial)
  
    grid$pred = predict(fit, newdata = grid, type = "response")
    wide = pivot_wider(grid, names_from = Diabetes, values_from = pred)
  
    p0 = sum(wide$`0` * std_wts)                      # std risk, Diab 0
    p1 = sum(wide$`1` * std_wts)                      # std risk, Diab 1
  
    log((p1/(1-p1)) / (p0/(1-p0)))                    # log(OR)
    }


  boot_out = boot(df_sample, statistic = std_or_fun,
                  R = 500, parallel = "no")
  ci = boot.ci(boot_out, type = "perc") 

  Model.std <- tibble(
    beta_Diabetes_std = boot_out$t0,
    beta_Diabetes_std_LCI = ci$percent[4],
    beta_Diabetes_std_UCI = ci$percent[5]
    )

  df.check = data.frame(
    #Summary statistics in the overall population
    # Mean of age in the overall population
    mean_age_all = mean(df_new$Age.2016), 
    
    # SD of age in the overall population
    sd_age_all = sd(df_new$Age.2016),     
    
    # Proportion of females in the overall population
    prop_female_all = mean(df_new$Female),        
    
    # Proportion of diabetes in the overall population
    prop_diabetes_all = mean(df_new$Diabetes), 
    
    # Mean of U in the overall population
    mean_U_all = mean(df_new$U),       
    
    # SD of U in the overall population
    sd_U_all = sd(df_new$U),     
    
    # Mean of diabetes duration among overall diabetic population
    mean_diabetes.duration_all = mean(df_new$Diabetes.duration[df_new$Diabetes == 1]), 
    
    # Incidence rate of diabetes in the overall population
    incidence.rate_diabetes = sum(df_new$Diabetes)/sum(df_new$FollowUptime), 
    
    # Proportion of survivals in the overall population
    prop_survival_all = mean(df_new$Survival),  
    
    # Mean of differences between generated and input marginal probabilities of death in each subgroup
    mean_perc_surv.diff = mean(df_surv_aggregated$perc_surv.diff),
    
    # Minimal probability of death
    min_p_death = min(df_new$prob_death_est),
    
    # Maximum probability of death
    max_p_death = max(df_new$prob_death_est),
    
    # Proportion of dementia in the overall population
    prop_dementia_all = mean(df_new$Dementia),                     
  
    # Minimal probability of dementia
    min_p_dementia = min(df_new$p_dementia),
    
    # Maximum probability of dementia
    max_p_dementia = max(df_new$p_dementia),
    
    #Summary statistics in the analytic sample
    # Mean of age in the analytic sample
    mean_age_sample = mean(df_sample$Age.2016),     
    
    # SD of age in the analytic sample
    sd_age_sample = sd(df_sample$Age.2016),   
    
    # Oldest age in the analytic sample
    max_age_sample = max(df_sample$Age.2016),  
    
    # Proportion of females in the analytic sample
    prop_female_sample = mean(df_sample$Female),        
    
    # Proportion of U in the analytic sample
    mean_U_sample = mean(df_sample$U),       
    
    # SD of U in the analytic sample
    sd_U_sample = sd(df_sample$U),          
    
    # Proportion of diabetes in the analytic sample
    prop_diabetes_sample = mean(df_sample$Diabetes),    
    
    # Mean of diabetes duration among diabetic sample
    mean_diabetes.duration_sample = mean(df_sample$Diabetes.duration[df_sample$Diabetes == 1]), 
    
    # Proportion of dementia in the analytic sample
    prop_dementia_sample = mean(df_sample$Dementia)           
    )

  return(fullresults = list(
    Model.sum_pop = Model.sum_pop,
    Model.sum_sample = Model.sum_sample,
    Model.std = Model.std,
    df.check = df.check))
  }

#### 5. Iterate the process for multiple times ####
multiple_output = function(N,C,I,seeds){
  stopifnot(N <= length(seeds))
  
  results <- parallel::mclapply(
    X = seq_len(N),
    FUN = function(i) {
      set.seed(seeds[i])
      ans = output(C, I)
      ans$df.check$seed = seeds[i]
      ans
    },
    mc.cores = 3,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE # respect the seed list
  )
  
  # Keep a record of which seed was used for each result
  names(results) <- paste0("seed_", seeds[seq_len(N)])
  attr(results, "seeds") <- seeds[seq_len(N)]
  
  saveRDS(results, file = sprintf("Results_Scene%s_%s_N%s_ChangingDiab.RDA", I, C, N))
  invisible(results)
  }


#### 6. Perform Analysis ####
setwd("/Users/kelvin/Library/CloudStorage/GoogleDrive-kelvinzh@umich.edu/My Drive/SEGA/Projects/Cross-national Simulation/Results/Results_20250711/5-Changing Diabetes and U by Country")
my_random_seeds = scan("/Users/kelvin/Library/CloudStorage/GoogleDrive-kelvinzh@umich.edu/My Drive/SEGA/Projects/Cross-national Simulation/Code/Code_20250711/seeds.txt", what = numeric(), sep = ",",quiet = TRUE)

# multiple_output(N = 250, "US", I=1)
# multiple_output(N = 250, "India", I=1)
# multiple_output(N = 250, "England", I=1)
# multiple_output(N = 250, "China", I=1)
# multiple_output(N = 250, "Mexico", I=1)


multiple_output(N = 250, "US", I=2, seeds = my_random_seeds)
multiple_output(N = 250, "India", I=2, seeds = my_random_seeds)
multiple_output(N = 250, "England", I=2, seeds = my_random_seeds)
multiple_output(N = 250, "China", I=2, seeds = my_random_seeds)
multiple_output(N = 250, "Mexico", I=2, seeds = my_random_seeds)

multiple_output(N = 250, "US", I=3, seeds = my_random_seeds)
multiple_output(N = 250, "India", I=3, seeds = my_random_seeds)
multiple_output(N = 250, "England", I=3, seeds = my_random_seeds)
multiple_output(N = 250, "China", I=3, seeds = my_random_seeds)
multiple_output(N = 250, "Mexico", I=3, seeds = my_random_seeds)


Outputlist <- list(
  England_2 = readRDS("Results_Scene2_England_N50_ChangingDiab.RDA"),
  US_2      = readRDS("Results_Scene2_US_N50_ChangingDiab.RDA"),
  Mexico_2  = readRDS("Results_Scene2_Mexico_N50_ChangingDiab.RDA"),
  China_2   = readRDS("Results_Scene2_China_N50_ChangingDiab.RDA"),
  India_2   = readRDS("Results_Scene2_India_N50_ChangingDiab.RDA"),
  
  England_3 = readRDS("Results_Scene3_England_N50_ChangingDiab.RDA"),
  US_3      = readRDS("Results_Scene3_US_N50_ChangingDiab.RDA"),
  Mexico_3  = readRDS("Results_Scene3_Mexico_N50_ChangingDiab.RDA"),
  China_3   = readRDS("Results_Scene3_China_N50_ChangingDiab.RDA"),
  India_3   = readRDS("Results_Scene3_India_N50_ChangingDiab.RDA")
)

countries <- c("England","US","Mexico","China","India")
scenarios <- c(2,3)
combinations <- expand.grid(C = countries, I = scenarios)

# Accumulate results
Final_results <- rbindlist(lapply(seq_len(nrow(combinations)), function(i) {
  C <- combinations$C[i]
  I <- combinations$I[i]
  key <- paste0(C, "_", I)
  Output <- Outputlist[[key]]
  
  Results <- rbindlist(lapply(Output, `[[`, "df.check"), idcol = "source")
  Results$Country <- C
  Results$Scenario <- I
  
  Results$beta_Diabetes_pop <- sapply(Output, function(x) x$Model.sum_pop[1, 1])
  Results$beta_Diabetes_pop_LCI <- sapply(Output, function(x) x$Model.sum_pop[1, 2])
  Results$beta_Diabetes_pop_UCI <- sapply(Output, function(x) x$Model.sum_pop[1, 3])
  Results$Coverage_Diabetes_pop <- ifelse(
    Results$beta_Diabetes_pop_LCI < log(2) & Results$beta_Diabetes_pop_UCI > log(2), 1, 0
  )
  
  Results$beta_Diabetes <- sapply(Output, function(x) x$Model.sum_sample[1, 1])
  Results$beta_Diabetes_LCI <- sapply(Output, function(x) x$Model.sum_sample[1, 2])
  Results$beta_Diabetes_UCI <- sapply(Output, function(x) x$Model.sum_sample[1, 3])
  Results$Coverage_Diabetes <- ifelse(
    Results$beta_Diabetes_LCI < log(2) & Results$beta_Diabetes_UCI > log(2), 1, 0
  )
  
  Results$beta_Diabetes_std <- sapply(Output, function(x) x$Model.std$beta_Diabetes_std)
  Results$beta_Diabetes_std_LCI <- sapply(Output, function(x) x$Model.std$beta_Diabetes_std_LCI)
  Results$beta_Diabetes_std_UCI <- sapply(Output, function(x) x$Model.std$beta_Diabetes_std_UCI)
  Results$Coverage_Diabetes_std <- ifelse(
    Results$beta_Diabetes_std_LCI < log(2) & Results$beta_Diabetes_std_UCI > log(2), 1, 0
  )
  return(Results)
}))

summary_stats <- Final_results[, .(
  mean_beta_Diabetes_pop = mean(beta_Diabetes_pop),
  sd_beta_Diabetes_pop = sd(beta_Diabetes_pop),
  Coverage_beta_Diabetes_pop = mean(Coverage_Diabetes_pop),
  
  mean_beta_Diabetes = mean(beta_Diabetes),
  sd_beta_Diabetes = sd(beta_Diabetes),
  Coverage_beta_Diabetes = mean(Coverage_Diabetes),
  
  mean_beta_Diabetes_std = mean(beta_Diabetes_std),
  sd_beta_Diabetes_std = sd(beta_Diabetes_std),
  Coverage_beta_Diabetes_std = mean(Coverage_Diabetes_std)
  ), 
  by = .(Country, Scenario)]

summary_stats$Bias_pop = (summary_stats$mean_beta_Diabetes_pop - (log(2)))/log(2)
summary_stats$Bias = (summary_stats$mean_beta_Diabetes - (log(2)))/log(2)
summary_stats$Bias_std = (summary_stats$mean_beta_Diabetes_std - (log(2)))/log(2)

writexl::write_xlsx(summary_stats, "summary_stats_ChangingDiab and U.xlsx")


