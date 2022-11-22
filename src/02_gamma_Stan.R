############################################################
#Aim: to fit a Bayesian hierarchical model (gamma dist) to MPX serial interval data in NL
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

stan_code_gamma <- '
  data{
    int N;
    int K;
    real X[N];
    int<lower=1, upper=K> GID[N];
  }
  
  parameters {
    real<lower=0> mu;
    real d_k[K];
    real<lower=0> s_r;
    real<lower=0> s_X;
  }
  
  transformed parameters{
    real<lower=0> mu_k[K];
    real<lower=0> shape[K];
    real<lower=0> rate[K];

    for (k in 1:K){
    mu_k[k] = mu + d_k[k];
    shape[k] = square(mu_k[k] / s_X);
    rate[k] = mu_k[k] / square(s_X);
    }
  }
  
  model{
    for(k in 1:K){
    d_k[k] ~ normal(0, s_r);
    }
    for(n in 1:N){
    X[n] ~ gamma(shape[GID[n]], rate[GID[n]]);
    }
    s_r ~ cauchy(0, 10);
    s_X ~ cauchy(0, 10);
  }
  
  generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] = gamma_lpdf(X[n] | shape[GID[n]], rate[GID[n]]);
    }
  }
'

###3. Model fit -----
fit_GGD_gamma <- stan(model_code=stan_code_gamma, 
                      data=data, 
                      seed=123, 
                      warmup = 1000, 
                      iter = 20000,
                      thin = 5,
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
#summary(fit_GGD_gamma); traceplot(fit_GGD_gamma, pars= paste("mu_k[", 1:data$K,"]",sep = "")); pairs(fit_GGD_gamma, pars = c("mu", "s_r","s_X", "lp__"))
#write_rds(fit_GGD_gamma,"fit_GGD_gamma.rds")
loo(extract_log_lik(fit_GGD_gamma))$loo
waic(extract_log_lik(fit_GGD_gamma))$waic

summary_gamma_df <- as_tibble(summary(fit_GGD_gamma)$summary, rownames = "par") %>% 
  # filter all par that start with 'mu'
  filter(grepl(par, pattern = "^mu")) %>% 
  select(par, mean, `2.5%`, `97.5%`) %>% 
  rename(lower = `2.5%`,
         upper = `97.5%`) %>% 
  mutate(GGD_index = as.character(0:data$K)) %>% #str_extract_all(par, pattern = "[0-9]", simplify = TRUE)[,1]) %>% 
  left_join(GGD_summary)

summary_gamma_df$GGD_index[summary_gamma_df$GGD_index=="0"] <- "All pooled"
write_rds(summary_gamma_df, "summary_gamma_df2.rds")

###4. posterior samples -----
extract_fit_GGD_gamma <- rstan::extract(fit_GGD_gamma)

pos_gamma_df <- tibble(mu.0 = extract_fit_GGD_gamma$mu,
                        mu = extract_fit_GGD_gamma$mu_k) %>% 
  do.call(data.frame, .) %>%
  as_tibble %>% 
  pivot_longer(cols = starts_with("mu"), names_to = "GGD_index", values_to = "pos", names_prefix = "mu.") %>% 
  left_join(GGD_summary %>% select(-n)) 

pos_gamma_df$GGD_index[pos_gamma_df$GGD_index=="0"] <- "All pooled"
write_rds(pos_gamma_df,"pos_gamma_df2.rds")