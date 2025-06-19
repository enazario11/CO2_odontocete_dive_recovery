#load libraries
library(tidyverse)
library(here)

#Dolphin data ####
dolph_resp <- read.csv(here("data/dolphins/VO2 Data Antis Sorted.csv")); dolph_resp$id <- paste0(dolph_resp$Animal, dolph_resp$Trial, dolph_resp$Date)
dolph_resp <- dolph_resp %>% filter(Test != 1)

dolph_vent <- read.csv(here("data/dolphins/Ventilation Duration Data Updated.csv")); dolph_vent$id <- paste0(dolph_vent$Animal, dolph_vent$Trial, dolph_vent$Date)
#filter out outliers (from antibiotics and general)
dolph_vent <- dolph_vent %>% 
  filter(id != "Rainx1 Swim70621" & id != "Rainx1 Swim120721" & id != "Rainx1 Swim32822" & id != "Rainx2 Swim112321" & id != "Rainx2 Swim40522" & id != "Donleyx1 Swim72821" & id != "Donleyx2 Swim90721" & id != "Donleyx3 Swim11322" & id != "Donleyx5 Swim20722")

dolph_blood <- read.csv(here("data/dolphins/i-STAT Data.csv")); dolph_blood$id <- paste0(dolph_blood$Animal, dolph_blood$Trial, dolph_blood$Date)
dolph_blood <- dolph_blood %>% filter(id != "RainRest52622" & id != "RainRest33122")

dolph_flir <- read.csv(here("data/dolphins/FLIR_dolphins_data.csv")); dolph_flir$id <- paste0(dolph_flir$dolphin, dolph_flir$session_type, dolph_flir$date); dolph_flir <- dolph_flir %>% rename("Animal" = "dolphin", "Trial" = "session_type", "Date" = "date")

#Beluga data ####
bel_resp <- read.csv(here("data/belugas/DL_Metabolics_Data.csv")); bel_resp$id <- paste0(bel_resp$Animal, bel_resp$Trial, bel_resp$Date)

bel_vent <- read.csv(here("data/belugas/DL_Vent_Data.csv")); bel_vent$id <- paste0(bel_vent$Animal, bel_vent$Trial, bel_vent$Date)

bel_blood <- read.csv(here("data/belugas/DL_Blood_Data.csv")); bel_blood$id <- paste0(bel_blood$Animal, bel_blood$Trial, bel_blood$Date)

bel_flir <- read.csv(here("data/belugas/DL_Vasc_Data.csv")); bel_flir$id <- paste0(bel_flir$Animal, bel_flir$Trial, bel_flir$Date)

#Combined DFs ####
## metabolics ####
all_resp <- dolph_resp %>% subset(select = -c(Speed_S1:Speed_S5, Test)) %>% mutate(sp = "Dolphin")
bel_temp <- bel_resp %>% subset(select = -c(Speed_S1:Speed_S5)) %>% mutate(sp = "Beluga")
all_resp <- rbind(all_resp, bel_temp)
all_resp$Animal <- factor(all_resp$Animal, levels = c("Maple", "Qinu", "Nunavik", "Shila", "Donley", "Rain"))
all_resp$Trial <- factor(all_resp$Trial, levels = c("Rest", "x1 Swim", "x2 Swim", "x3 Swim", "x5 Swim", "SAB"))
all_resp$sp <- factor(all_resp$sp, levels = c("Beluga", "Dolphin"), labels = c("Beluga", "Dolphin"))

## ventilation ####
dolph_temp <- dolph_vent %>% mutate(sp = "Dolphin")
bel_temp <- bel_vent %>% mutate(sp = "Beluga")
all_vent <- rbind(dolph_temp, bel_temp)
all_vent$Animal <- factor(all_vent$Animal, levels = c("Maple", "Qinu","Nunavik", "Shila", "Donley", "Rain"))
all_vent$Trial <- factor(all_vent$Trial, levels = c("Rest", "x1 Swim", "x2 Swim", "x3 Swim", "x5 Swim", "SAB"))
all_vent$sp <- factor(all_vent$sp, levels = c("Beluga", "Dolphin"), labels = c("Beluga", "Dolphin"))

  #calculate elapsed time into surface recovery and breath frequency per 60 sec bins
#elapsed time
dat_freq <- all_vent %>%
  group_by(id) %>%
  arrange(id) %>% 
  ungroup()
dat_freq$sp <- factor(dat_freq$sp, levels = c("Dolphin", "Beluga"))
dat_freq <- dat_freq %>%
  group_by(id) %>%
  filter(Surface_Number == max(Surface_Number))

dat_freq$elapsed <- NA
dat_freq_rec <- NULL
for(i in 1:length(unique(dat_freq$id))){
  ids <- unique(dat_freq$id)
  curr_id = ids[i]
  df_temp <- dat_freq %>% filter(id == curr_id)
  
  df_temp$elapsed = df_temp$Start_Time - min(df_temp$Start_Time, na.rm = TRUE)

  dat_freq_rec <- rbind(dat_freq_rec, df_temp)
}

#breath frequency
binwidth_s <- 30
vent_binned <- dat_freq_rec %>% 
  filter(Breath_Duration > 0) %>% #filters out typo where NA should have been
  ungroup() %>% 
  group_by(sp, id) %>% 
  mutate(elapsed_bin = elapsed %/% binwidth_s) %>% 
  group_by(sp, id, Animal, Date, Trial, elapsed_bin) %>% 
  summarize(breath_freq_min = n() * 60 / binwidth_s,
            mean_ibi = mean(IBI, na.rm = TRUE),
            mean_breath_dur = mean(Breath_Duration, na.rm = TRUE),
            .groups = "drop")

vent_binned$breath_freq_min[is.na(vent_binned$breath_freq_min)] <- 0
dat_freq_rec <- vent_binned %>% filter(Trial != "Rest" & Trial != "x2 Swim") %>% mutate(elapsed_time = (elapsed_bin*30)/60)


## blood ####
dolph_temp <- dolph_blood %>% mutate(sp = "Dolphin") %>% subset(select = -c(Dur_S1:Dur_S5, Fasted, SAB_num, X))
bel_temp <- bel_blood %>% mutate(sp = "Beluga") %>% filter(Breathhold == "N") %>% subset(select = -c(Dur_S1:Bow_rate, Breathhold))
all_blood <- rbind(dolph_temp, bel_temp)
all_blood$Animal <- factor(all_blood$Animal, levels = c("Maple", "Qinu","Nunavik", "Shila", "Whisper", "Donley", "Rain"))
all_blood$Trial <- factor(all_blood$Trial, levels = c("Rest", "x1 Swim", "x3 Swim", "x5 Swim", "SAB"))
all_blood$sp <- factor(all_blood$sp, levels = c("Beluga", "Dolphin"), labels = c("Beluga", "Dolphin"))

## vasculature ####
dolph_temp <- dolph_flir %>% mutate(sp = "Dolphin") %>% subset(select= -c(X, X.1, n_SAB))
bel_temp <- bel_flir %>% mutate(sp = "Beluga") %>% subset(select = -c(bow_dur))
all_vasc <- rbind(dolph_temp, bel_temp)
all_vasc$Animal <- factor(all_vasc$Animal, levels = c("Maple", "Qinu","Nunavik", "Shila", "Whisper", "Donley", "Rain"))
all_vasc$Trial <- factor(all_vasc$Trial, levels = c("Rest", "x1 Swim", "x3 Swim", "x5 Swim", "SAB"))
all_vasc$sp <- factor(all_vasc$sp, levels = c("Beluga", "Dolphin"), labels = c("Beluga", "Dolphin"))

dolph_temp <- all_vasc %>% filter(sp == "Dolphin")
bel_temp <- all_vasc %>% filter(sp == "Beluga")

#establish all beluga resting points
bel_rest <- bel_temp %>% 
  filter(Trial == "Rest" | img_num == "rest") %>%
  mutate(Trial = "Rest")

bel_temp2 <- bel_temp %>%
  filter(Trial != "Rest" & img_num != "rest")

all_vasc <- rbind(bel_rest, bel_temp2, dolph_temp)

dat_f <- all_vasc %>%
  mutate(point_mean = as.numeric(point_mean), 
         surf_mean = as.numeric(surf_mean), 
         prop_vasc = vasc_area/fluke_area, 
         max3 = as.numeric(max3)) %>%
  rowwise() %>%
  mutate(max_max = max(c_across(max1:max5))) %>%
  ungroup()

dat_f$Trial <- factor(dat_f$Trial, levels = c("Rest", "x1 Swim", "x3 Swim", "x5 Swim", "SAB"))

#standardize bel prop area by resting prop area for the session
dat_f$std_area = NA
dolph_temp <- dat_f %>% filter(sp == "Dolphin") %>% mutate(std_area = prop_vasc)

bel_temp <- dat_f %>% filter(sp == "Beluga")
bel_mod <- NULL
for(i in 1:length(unique(bel_temp$id))){
  #group images by id (animal_trial_date)
  id_curr = unique(bel_temp$id)[i]
  df_curr = bel_temp %>% filter(id == id_curr)
  
  #identify resting area
  if(any(str_detect(df_curr$img_num, "rest")) == FALSE){
    df_curr <- df_curr %>% mutate(std_area = NA)
  }
  
  if(any(str_detect(df_curr$img_num, "rest")) == TRUE){
    rest_area <- df_curr %>% filter(img_num == "rest") 
    not_rest_area <- df_curr %>% filter(img_num != "rest")
    
    #standardize by the resting area value
    rest_val <- rest_area$prop_vasc
    
    if(rest_val > 0){
      not_rest_area <- not_rest_area %>% mutate(std_area = prop_vasc / rest_val)
      rest_area <- rest_area %>% mutate(std_area = prop_vasc / rest_val)
      
      df_curr <- rbind(rest_area, not_rest_area)
    }
    
    if(rest_val == 0){
      not_rest_area <- not_rest_area %>% mutate(std_area = NA)
      rest_area <- rest_area %>% mutate(std_area = NA)
      
      df_curr <- rbind(rest_area, not_rest_area)
    }
    
  }
  
  bel_mod <- rbind(bel_mod, df_curr)
  
}

dat_f <- rbind(dolph_temp, bel_mod)

## Lap speed ####
dolph_speed <- dolph_resp %>% subset(select = c(Animal, Date, Trial, Speed_S1:Bow_rate)) %>% mutate(sp = "Dolphin")

blood_speed <- bel_blood %>% subset(select =  c(Animal, Date, Trial, Speed_S1, Speed_S2, Speed_S3, Bow_rate)) %>% mutate(sp = "Beluga")

bel_speed <- bel_resp %>% subset(select = c(Animal:Trial, Speed_S1:Bow_rate)) %>% mutate(sp = "Beluga")
blood_speed$Speed_S4 = "NA"
blood_speed$Speed_S5 = "NA"
bel_speed <- rbind(bel_speed, blood_speed)

all_speed <- rbind(dolph_speed, bel_speed)

