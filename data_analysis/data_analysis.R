#load libraries
library(tidyverse)
library(here)
library(glmmTMB)
library(DHARMa)
library(nlme)
library(emmeans)
source("data_visualization/data_visualization.R")
set.seed(0904)

## Metabolics ####
### gas type & respiratory gas recovery time GLM ####
resp_mod_df <- all_resp %>% 
  pivot_longer(cols = O2_tobase:CO2_tobase,
               names_to = "gas_type",
               values_to = "recover_time"
  ) %>% 
  subset(select = -c(Bow_rate))

resp_mod_df_dl <- resp_mod_df %>% filter(sp == "Beluga")
resp_mod_df_tt <- resp_mod_df %>% filter(sp == "Dolphin")

resp_mod_dl <- glm(recover_time ~ gas_type + Trial + Animal + Trial*gas_type + Trial*Animal, 
                   data = resp_mod_df_dl, 
                   family = Gamma(link ="inverse"))
resp_mod_tt <- glm(recover_time ~ gas_type + Trial + Animal + Trial*gas_type + Trial*Animal, 
                   data = resp_mod_df_tt, 
                   family = Gamma(link = "inverse"))

### LRT: gas type ####
resp_mod_Ngas_dl <- glm(recover_time ~ Trial + Animal + Trial*Animal,
                        data = resp_mod_df_dl, 
                        family = Gamma(link = "inverse"))

likelihood_ratio_test_dl <- lmtest::lrtest(resp_mod_Ngas_dl, resp_mod_dl)

if (likelihood_ratio_test_dl$"Pr(>Chisq)"[2] < 0.05) {
  cat("Reject the null hypothesis. The full model is significantly better than 
        the null model.\n")
} else {
  cat("Fail to reject the null hypothesis. The null model is sufficient.\n")
}

resp_mod_Ngas_tt <- glm(recover_time ~ Trial + Animal + Trial*Animal,
                        data = resp_mod_df_tt, 
                        family = Gamma(link = "inverse"))

likelihood_ratio_test_tt <- lmtest::lrtest(resp_mod_Ngas_tt, resp_mod_tt)

if (likelihood_ratio_test_tt$"Pr(>Chisq)"[2] < 0.05) {
  cat("Reject the null hypothesis. The full model is significantly better than 
        the null model.\n")
} else {
  cat("Fail to reject the null hypothesis. The null model is sufficient.\n")
}

### species recovery times GLMM ####
resp_mod_df_all <- resp_mod_df %>% filter(Trial == "x3 Swim" | Trial == "SAB")
resp_mod_df_all$sp <- relevel(resp_mod_df_all$sp, ref = "Dolphin")

resp_mod_all <- glmmTMB(recover_time ~ Trial*gas_type + Trial*sp + gas_type*sp + (1|Animal), 
                        data = resp_mod_df_all, 
                        family = Gamma(link ="log"))

recov_coefs <- MASS::mvrnorm(n = 1000, 
                             mu = fixef(resp_mod_all)$cond,
                             Sigma = vcov(resp_mod_all)$cond) %>%
  as.data.frame()

head(recov_coefs)
exp(mean(recov_coefs$spBeluga))
exp(sd(recov_coefs$spBeluga))

## Blood gases ####
blood_at_draw0_dl <- function(blood_coefs) { #create grid of pH values to predict for when draw time is held at 0
  grid <- expand_grid(
    Draw_Time = 0,
    Trial = fct_drop(unique(all_blood_dl$Trial)),
    Animal = fct_drop(unique(all_blood_dl$Animal))
  )
  X <- model.matrix(~ Draw_Time + Trial + Animal + Draw_Time:Trial, 
                    data = grid)
  
  #grid$pH <- exp(mean(X %*% ph_coefs)) #calculates the response variable value on the response scale
  grid$log_blood <- (X %*% blood_coefs)[, 1] 
  
  result <- grid %>% 
    group_by(Trial) %>%
    summarize(mean_blood_log = mean(log_blood), 
              blood_intercept = exp(mean_blood_log), 
              .groups = "drop_last") %>%
    pull(blood_intercept)
  
  names(result) <- unique(grid$Trial)
  result
}
blood_recov_dl <- function(blood_coefs) {
  grid <- expand_grid(
    Draw_Time = c(0, 1000), #gets start and end times (first and last time)
    Trial = fct_drop(unique(all_blood_dl$Trial)),
    Animal = fct_drop(unique(all_blood_dl$Animal))
  )
  X <- model.matrix(~ Draw_Time + Trial + Animal + Draw_Time:Trial, 
                    data = grid)
  grid$log_blood <- (X %*% blood_coefs)[, 1] #just taking the y-value (pH)
  grid$blood <- exp(grid$log_blood)
  
  result <- grid %>% 
    group_by(Trial, Animal) %>% 
    summarize(dblood_dt = (blood[2] - blood[1]) / (1000 / 60),
              .groups = "drop_last") %>% 
    summarize(dblood_dt = mean(dblood_dt)) %>% 
    pull(dblood_dt)
  names(result) <- unique(grid$Trial)
  result
}
blood_at_draw0_tt <- function(blood_coefs) { #create grid of pH values to predict for when draw time is held at 0
  grid <- expand_grid(
    Draw_Time = 0,
    Trial = fct_drop(unique(all_blood_tt$Trial)),
    Animal = fct_drop(unique(all_blood_tt$Animal))
  )
  X <- model.matrix(~ Draw_Time + Trial + Animal + Draw_Time:Trial, 
                    data = grid)
  
  grid$log_blood <- (X %*% blood_coefs)[, 1] 
  
  result <- grid %>% 
    group_by(Trial) %>%
    summarize(mean_blood_log = mean(log_blood), 
              blood_intercept = exp(mean_blood_log), 
              .groups = "drop_last") %>%
    pull(blood_intercept)
  
  names(result) <- levels(grid$Trial)
  result
}
blood_recov_tt <- function(blood_coefs) {
  grid <- expand_grid(
    Draw_Time = c(0, 600), #gets start and end times (first and last time) -- different than belugas
    Trial = fct_drop(unique(all_blood_tt$Trial)),
    Animal = fct_drop(unique(all_blood_tt$Animal))
  )
  X <- model.matrix(~ Draw_Time + Trial + Animal + Draw_Time:Trial, 
                    data = grid)
  grid$log_blood <- (X %*% blood_coefs)[, 1] #just taking the y-value (pH)
  grid$blood <- exp(grid$log_blood)
  
  result <- grid %>% 
    group_by(Trial, Animal) %>% 
    summarize(dblood_dt = (blood[2] - blood[1]) / (600 / 60),
              .groups = "drop_last") %>% 
    summarize(dblood_dt = mean(dblood_dt)) %>% 
    pull(dblood_dt)
  
  names(result) <- levels(grid$Trial)
  result
}

### pO2/pCO2 and time at surface GLM ####
all_blood_dl <- all_blood_dl %>% droplevels()
all_blood_tt <- all_blood_tt %>% droplevels()

all_blood_dl$Animal_Date <- interaction(all_blood_dl$Animal, all_blood_dl$Date, drop = TRUE)

#belugas O2
po2_rec_dl <- glm(PO2 ~ Draw_Time * Trial + Animal, 
                  data = all_blood_dl, 
                  family = Gamma(link = "log"))

#predict to get average rate of change and y-int
o2_dl_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estimate coefs of beluga ph recovery model
                             mu = coefficients(po2_rec_dl),
                             Sigma = vcov(po2_rec_dl))

simulated_o2_yint_dl <- apply(o2_dl_coefs, MARGIN = 1, FUN = blood_at_draw0_dl) #runs above function for each of the 1000 sims
o2_yint_mean_dl <- rowMeans(simulated_o2_yint_dl)
o2_yint_sd_dl <- apply(simulated_o2_yint_dl, 1, sd) #gets sd by trial type
o2_yint_mean_dl - 1.96 * o2_yint_sd_dl #CI low
o2_yint_mean_dl + 1.96 * o2_yint_sd_dl #CI high

simulated_o2_recov_dl <- apply(o2_dl_coefs, MARGIN = 1, FUN = blood_recov_dl)
o2_recov_mean_dl <- rowMeans(simulated_o2_recov_dl) #gets mean by trial type
o2_recov_sd_dl<- apply(simulated_o2_recov_dl, 1, sd) #gets sd by trial type
o2_recov_mean_dl - 1.96 * o2_recov_sd_dl #CI low
o2_recov_mean_dl + 1.96 * o2_recov_sd_dl #CI high


#belugas CO2
pco2_rec_dl <- glm(PCO2 ~ Draw_Time * Trial + Animal, 
                   data = all_blood_dl, 
                   family = Gamma(link = "log"))

#predict to get average rate of change and y-int
co2_dl_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estimate coefs of beluga ph recovery model
                              mu = coefficients(pco2_rec_dl),
                              Sigma = vcov(pco2_rec_dl))
#yint
simulated_co2_yint_dl <- apply(co2_dl_coefs, MARGIN = 1, FUN = blood_at_draw0_dl) #runs above function for each of the 1000 sims
co2_yint_mean_dl <- rowMeans(simulated_co2_yint_dl)
co2_yint_sd_dl <- apply(simulated_co2_yint_dl, 1, sd) #gets sd by trial type
co2_yint_mean_dl - 1.96 * co2_yint_sd_dl #CI low
co2_yint_mean_dl + 1.96 * co2_yint_sd_dl #CI high

#rate of change
simulated_co2_recov_dl <- apply(co2_dl_coefs, MARGIN = 1, FUN = blood_recov_dl)
co2_recov_mean_dl <- rowMeans(simulated_co2_recov_dl) #gets mean by trial type
co2_recov_sd_dl<- apply(simulated_co2_recov_dl, 1, sd) #gets sd by trial type
co2_recov_mean_dl - 1.96 * co2_recov_sd_dl #CI low
co2_recov_mean_dl + 1.96 * co2_recov_sd_dl #CI high

#dolphins O2
po2_rec_tt <- glm(PO2 ~ Draw_Time * Trial + Animal, 
                  data = all_blood_tt, 
                  family = Gamma(link = "log"))

o2_tt_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estiate coefs of beluga ph recovery model
                             mu = coefficients(po2_rec_tt),
                             Sigma = vcov(po2_rec_tt))

#y-int
simulated_o2_yint_tt <- apply(o2_tt_coefs, MARGIN = 1, FUN = blood_at_draw0_tt) #runs above function for each of the 1000 sims
o2_yint_mean_tt <- rowMeans(simulated_o2_yint_tt)
o2_yint_sd_tt <- apply(simulated_o2_yint_tt, 1, sd) #gets sd by trial type
o2_yint_mean_tt - 1.96 * o2_yint_sd_tt #CI low
o2_yint_mean_tt + 1.96 * o2_yint_sd_tt #CI high

#rate of change
simulated_o2_recov_tt <- apply(o2_tt_coefs, MARGIN = 1, FUN = blood_recov_tt)
o2_recov_mean_tt <- rowMeans(simulated_o2_recov_tt) #gets mean by trial type
o2_recov_sd_tt <- apply(simulated_o2_recov_tt, 1, sd) #gets sd by trial type
o2_recov_mean_tt - 1.96 * o2_recov_sd_tt #CI low
o2_recov_mean_tt + 1.96 * o2_recov_sd_tt #CI high

#dolphins CO2
pco2_rec_tt <- glm(PCO2 ~ Draw_Time * Trial + Animal, 
                   data = all_blood_tt, 
                   family = Gamma(link = "log"))

co2_tt_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estiate coefs of beluga ph recovery model
                              mu = coefficients(pco2_rec_tt),
                              Sigma = vcov(pco2_rec_tt))

#y-int
simulated_co2_yint_tt <- apply(co2_tt_coefs, MARGIN = 1, FUN = blood_at_draw0_tt) #runs above function for each of the 1000 sims
co2_yint_mean_tt <- rowMeans(simulated_co2_yint_tt)
co2_yint_sd_tt <- apply(simulated_co2_yint_tt, 1, sd) #gets sd by trial type
co2_yint_mean_tt - 1.96 * co2_yint_sd_tt #CI low
co2_yint_mean_tt + 1.96 * co2_yint_sd_tt #CI high

#rate of change
simulated_co2_recov_tt <- apply(co2_tt_coefs, MARGIN = 1, FUN = blood_recov_tt)
co2_recov_mean_tt <- rowMeans(simulated_co2_recov_tt) #gets mean by trial type
co2_recov_sd_tt <- apply(simulated_co2_recov_tt, 1, sd) #gets sd by trial type
co2_recov_mean_tt - 1.96 * co2_recov_sd_tt #CI low
co2_recov_mean_tt + 1.96 * co2_recov_sd_tt #CI high

## Blood pH ####
### pCO2 vs. lactate and blood pH GLM ####
blood_dl <- all_blood %>% filter(sp == "Beluga")
blood_tt <- all_blood %>% filter(sp == "Dolphin")

co2_mod_dl <- glm(pH ~ PCO2 + Trial + Animal + Trial*PCO2, 
                  data = blood_dl, 
                  family = gaussian(link = "log"))

co2_mod_tt <- glm(pH ~ PCO2 + Trial + Animal + Trial*PCO2, 
                  data = blood_tt, 
                  family = gaussian(link = "log"))

lact_mod_dl <- glm(pH ~ Lactate + Trial + Animal + Trial*Lactate, 
                   data = blood_dl, 
                   family = gaussian(link = "log"))

lact_mod_tt <- glm(pH ~ Lactate + Trial + Animal + Trial*Lactate, 
                   data = blood_tt, 
                   family = gaussian(link = "log"))

combo_mod_dl <- glm(pH ~ Lactate + PCO2 + Trial + Animal + Trial*Lactate + Trial*PCO2, 
                    data = blood_dl, 
                    family = gaussian(link = "log"))

combo_mod_tt <- glm(pH ~ Lactate + PCO2 + Trial + Animal + Trial*Lactate + Trial*PCO2, 
                    data = blood_tt, 
                    family = gaussian(link = "log"))

#compare performance across models
performance::compare_performance(co2_mod_dl, lact_mod_dl, combo_mod_dl, co2_mod_tt, lact_mod_tt, combo_mod_tt)

#test if inclusion of CO2 in beluga model made a difference in combo model (compared to just lact mod). The combined model with CO2 performs sig better. 
likelihood_ratio_test1 <- lmtest::lrtest(combo_mod_dl, lact_mod_dl)
if (likelihood_ratio_test1$"Pr(>Chisq)"[2] < 0.05) {
  cat("Reject the null hypothesis. The full model is significantly better than 
        the null model.\n")
} else {
  cat("Fail to reject the null hypothesis. The null model is sufficient.\n")
}

likelihood_ratio_test2 <- lmtest::lrtest(combo_mod_dl, co2_mod_dl) #combo model with lact performs sig better than just co2 model
if (likelihood_ratio_test2$"Pr(>Chisq)"[2] < 0.05) {
  cat("Reject the null hypothesis. The full model is significantly better than 
        the null model.\n")
} else {
  cat("Fail to reject the null hypothesis. The null model is sufficient.\n")
}

### blood pH and time at surface GLM ####
#belugas pH
ph_rec_dl <- glm(pH ~ Draw_Time * Trial + Animal, #can add interaction term between trial and draw time if want to get rate of change for each trial type
                 data = all_blood_dl, 
                 family = Gamma(link = "log"))

#predict to get average rate of change and y-int
ph_dl_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estimate coefs of beluga ph recovery model
                             mu = coefficients(ph_rec_dl),
                             Sigma = vcov(ph_rec_dl))

#yint
simulated_ph0_dl <- apply(ph_dl_coefs, MARGIN = 1, FUN = blood_at_draw0_dl) #runs above function for each of the 1000 sims
ph_yint_mean_dl <- rowMeans(simulated_ph0_dl)
ph_yint_sd_dl <- apply(simulated_ph0_dl, 1, sd) #gets sd by trial type
ph_yint_mean_dl - 1.96 * ph_yint_sd_dl #CI low
ph_yint_mean_dl + 1.96 * ph_yint_sd_dl #CI high

#rate of change
simulated_ph_recov_dl <- apply(ph_dl_coefs, MARGIN = 1, FUN = blood_recov_dl)
ph_recov_mean_dl <- rowMeans(simulated_ph_recov_dl) #gets mean by trial type
ph_recov_sd_dl <- apply(simulated_ph_recov_dl, 1, sd) #gets sd by trial type
ph_recov_mean_dl - 1.96 * ph_recov_sd_dl #CI low
ph_recov_mean_dl + 1.96 * ph_recov_sd_dl #CI high

#dolphins pH
ph_rec_tt <- glm(pH ~ Draw_Time * Trial + Animal, 
                 data = all_blood_tt, 
                 family = Gamma(link = "log"))

#predict to get average rate of change and y-int
ph_tt_coefs <- MASS::mvrnorm(n = 1000, #number of simulations to estiate coefs of beluga ph recovery model
                             mu = coefficients(ph_rec_tt),
                             Sigma = vcov(ph_rec_tt))

#y-int
simulated_ph0_tt <- apply(ph_tt_coefs, MARGIN = 1, FUN = blood_at_draw0_tt) #runs above function for each of the 1000 sims
ph_yint_mean_tt <- rowMeans(simulated_ph0_tt)
ph_yint_sd_tt <- apply(simulated_ph0_tt, 1, sd) #gets sd by trial type
ph_yint_mean_tt - 1.96 * ph_yint_sd_tt #CI low
ph_yint_mean_tt + 1.96 * ph_yint_sd_tt #CI high

#rate of change
simulated_ph_recov_tt <- apply(ph_tt_coefs, MARGIN = 1, FUN = blood_recov_tt)
ph_recov_mean_tt <- rowMeans(simulated_ph_recov_tt) #gets mean by trial type
ph_recov_sd_tt <- apply(simulated_ph_recov_tt, 1, sd) #gets sd by trial type
ph_recov_mean_tt - 1.96 * ph_recov_sd_tt #CI low
ph_recov_mean_tt + 1.96 * ph_recov_sd_tt #CI high

## Ventilation ####
### ventilation metrics and time at surface asymptotic models ####
vent_dl <- dat_freq_rec %>% filter(sp == "Beluga" & elapsed_time < 10) %>% droplevels() %>% drop_na()
vent_dl$mean_ibi[vent_dl$mean_ibi == "NaN"] <- NA

vent_tt <- dat_freq_rec %>% filter(sp == "Dolphin" & elapsed_time < 10) %>% droplevels() %>% drop_na()

#beluga models
length_k <- 1 + (length(unique(vent_dl$Animal))-1) + (length(unique(vent_dl$Trial))-1)

#breath frequency
vent_rate_dl <- gnls(
  model = breath_freq_min ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_dl, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = bel_rest_hline_f, 
    R0 = 6, 
    k = rep(0.1, length_k)), 
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
)

#breath duration 
vent_dur_dl <- gnls(
  model = mean_breath_dur ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_dl, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = bel_rest_hline_d, 
    R0 = 1.5, 
    k = rep(0.1, length_k)), 
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
)

#IBI 
vent_ibi_dl <- gnls(
  model = mean_ibi ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_dl, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = bel_rest_hline_i, 
    R0 = 8, 
    k = rep(0.1, length_k)), 
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
)

# dolphin models
length_k <- 1+ (length(unique(vent_tt$Animal))-1) + (length(unique(vent_tt$Trial))-1)

#breath frequency
vent_rate_tt <- gnls(
  model = breath_freq_min ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_tt, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = dol_rest_hline_f, 
    R0 = 6, 
    k = rep(0.1, length_k)), 
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
)

#breath duration 
vent_dur_tt <- gnls(
  model = mean_breath_dur ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_tt, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = dol_rest_hline_d, 
    R0 = 1, 
    k = rep(0.1, length_k)), 
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
)

#IBI
vent_ibi_tt <- gnls(
  model = mean_ibi ~ Asym + (R0 - Asym) * exp(-k * elapsed_time), 
  data = vent_tt, 
  
  params = list(
    Asym ~ 1, 
    R0 ~ 1,
    k ~ Animal + Trial
  ),
  
  start = list(
    Asym = dol_rest_hline_i, 
    R0 = 12, 
    k = rep(0.001, length_k)),
  
  control = gnlsControl(
    maxIter = 100, 
    tolerance = 1e-5, 
    nlsTol = 1e-5, 
    minScale = 1e-10),
  
  correlation = corAR1(form = ~ elapsed_time | Date/Animal)
  
)

## Vasculature ####
### vasculature metrics and time at surface GLMM ####
vasc_dl <- dat_f %>% 
  filter(sp == "Beluga" & Animal != "Whisper" & Trial != "Rest") %>% 
  droplevels() %>%
  group_by(Animal, Trial, Date) %>%
  mutate(img_num = as.numeric(img_num)) %>%
  arrange(img_num, .by_group = TRUE)

vasc_tt <- dat_f %>% 
  mutate(img_num = as.numeric(img_num)) %>%
  filter(sp == "Dolphin" & Trial != "Rest") %>% 
  filter(Date != 60723 & Date != 71223 & Date != 62723 & Date != 71123 & Date != 62823) %>% #filter out dates to get regular sampling
  filter(img_num < 10) %>% #filter out extended samples to support regularly spaced autocorrelation
  droplevels() %>%
  group_by(Animal, Trial, Date) %>%
  arrange(img_num, .by_group = TRUE)

# beluga
#mean temp ~ trial type 
vasc_mean_dl <- glmmTMB(surf_mean ~ Trial*img_num + Animal + ar1(factor(img_num) + 0 | Date/Animal), 
                        data = vasc_dl, 
                        family = Gamma(link = "log"))

#max temp ~ trial type
vasc_max_dl <- glmmTMB(max_max ~ Trial*img_num + Animal + ar1(factor(img_num) + 0 | Date/Animal), 
                       data = vasc_dl, 
                       family = Gamma(link = "log"))

#percent perfusion ~ trial type 
vasc_dl$prop_vasc[vasc_dl$prop_vasc == 0] <- 0 # change to just above 0 for beta distribution in model
vasc_perf_dl <- glmmTMB(prop_vasc ~ Trial*img_num + ar1(factor(img_num) + 0 | Date/Animal), 
                        ziformula = ~Trial, #models the probability of 0 by trial type
                        data = vasc_dl, 
                        family = beta_family(link = "logit"))
mean_tukey <- emmeans::emmeans(vasc_mean_tt, ~Trial)
pairs(mean_tukey, adjust = "bonferroni")

# dolphin
#mean temp ~ trial type 
vasc_mean_tt <- glmmTMB(surf_mean ~ Trial*img_num + Animal + ar1(factor(img_num) + 0 | Date/Animal), 
                        data = vasc_tt, 
                        family = Gamma(link = "log"))

mean_tukey <- emmeans::emmeans(vasc_mean_tt, ~Trial)
pairs(mean_tukey, adjust = "bonferroni")

#max temp ~ trial type 
vasc_max_tt <- glmmTMB(max_max ~ Trial*img_num + Animal + ar1(factor(img_num) + 0 | Date/Animal), 
                       data = vasc_tt, 
                       family = Gamma(link = "log"))

max_tukey <- emmeans::emmeans(vasc_max_tt, ~Trial)
pairs(max_tukey, adjust = "bonferroni")

#percent perfusion ~ trial type 
vasc_tt$prop_vasc[vasc_tt$prop_vasc == 0] <- 0 # change to just above 0 for beta distribution in model
vasc_perf_tt <- glmmTMB(prop_vasc ~ Trial*img_num + Animal + ar1(factor(img_num) + 0 | Date/Animal), 
                        ziformula = ~Trial, #models the probability of 0 by trial type
                        data = vasc_tt, 
                        family = beta_family(link = "logit"))

perf_tukey <- emmeans::emmeans(vasc_perf_tt, ~Trial)
pairs(perf_tukey, adjust = "bonferroni")

