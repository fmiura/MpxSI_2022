############################################################
#Aim: to fit a Bayesian hierarchical model (normal dist) to MPX serial interval data in NL
#Final edit: 22 November 2022
#Editor: Fumi Miura & Jantien Backer
############################################################
###Procedure
#0. Package 
#1. NL data 
#2. Stan code
#3. Model fit
#4. Visualization
############################################################

###0. Package -----
library(rstan)
library(tidyverse)
library(glue)
library(loo)

###1. NL data -----
#raw data
GGD_raw_data <- read_csv("GGD_raw_data.csv")
GGD_summary <- GGD_raw_data %>% 
  mutate(GGD_index = as.character(GGD_index)) %>% 
  group_by(GGD_index) %>% 
  count() %>% 
  ungroup() %>% 
  add_row(GGD_index = "0", 
          n = nrow(GGD_raw_data),
          .before = TRUE)

#input data for Stan
data <- list(
  N = nrow(GGD_raw_data), #no. of observed SI values
  K = max(GGD_raw_data$GGD_index), #no. of GGD
  X = GGD_raw_data$AbsSI, #reported serial interval values
  GID = GGD_raw_data$GGD_index #GGD id (=each GGD has different ID)
)

###2. Stan code -----
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

stan_normal <- '
  data{
    int N;
    int K;
    real X[N];
    int<lower=1, upper=K> GID[N];
  }
  
  parameters {
    real mu;
    real d_k[K];
    real<lower=0> s_r;
    real<lower=0> s_X;
  }
  
  transformed parameters{
    real mu_k[K];
    for (k in 1:K){
    mu_k[k] = mu + d_k[k];
    }
  }
  
  model{
    for(k in 1:K){
    d_k[k] ~ normal(0, s_r);
    }
    for(n in 1:N){
    X[n] ~ normal(mu_k[GID[n]], s_X);
    }
    s_r ~ cauchy(0, 10);
    s_X ~ cauchy(0, 10);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = normal_lpdf(X[n] | mu_k[GID[n]], s_X);
    }
  }
'

###3. Model fit -----
fit_GGD_normal <- stan(model_code=stan_normal, 
                       data=data, 
                       seed=123, 
                       warmup = 1000, 
                       iter = 20000,
                       thin = 5,
                       control = list(adapt_delta = 0.975, max_treedepth = 15))

#summary(fit_GGD_normal);traceplot(fit_GGD_normal, pars= paste("mu_k[", 1:data$K,"]",sep = "")); pairs(fit_GGD_normal, pars = c("mu", "s_r","s_X", "lp__"))
#write_rds(fit_GGD_normal,"fit_GGD_normal.rds")
loo(extract_log_lik(fit_GGD_normal))$loo
waic(extract_log_lik(fit_GGD_normal))$waic

summary_normal_df <- as_tibble(summary(fit_GGD_normal)$summary, rownames = "par") %>% 
  # filter all par that start with 'mu'
  filter(grepl(par, pattern = "^mu")) %>% 
  select(par, mean, `2.5%`, `97.5%`) %>% 
  rename(lower = `2.5%`,
         upper = `97.5%`) %>% 
  mutate(GGD_index = as.character(0:data$K)) %>% #str_extract_all(par, pattern = "[0-9]", simplify = TRUE)[,1]) %>% 
  left_join(GGD_summary)

summary_normal_df$GGD_index[summary_normal_df$GGD_index=="0"] <- "All pooled"
write_rds(summary_normal_df, "summary_normal_df2.rds")

###4. posterior samples -----
extract_fit_GGD_normal <- rstan::extract(fit_GGD_normal)

pos_normal_df <- tibble(mu.0 = extract_fit_GGD_normal$mu,
                        mu = extract_fit_GGD_normal$mu_k) %>% 
  do.call(data.frame, .) %>%
  as_tibble %>% 
  pivot_longer(cols = starts_with("mu"), names_to = "GGD_index", values_to = "pos", names_prefix = "mu.") %>% 
  left_join(GGD_summary %>% select(-n)) 

pos_normal_df$GGD_index[pos_normal_df$GGD_index=="0"] <- "All pooled"
write_rds(pos_normal_df,"pos_normal_df2.rds")