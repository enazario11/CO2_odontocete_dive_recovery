#load libraries
library(tidyverse)
library(here)
library(tidyquant)
library(ggpubr)
library(kableExtra)
library(tinytex)
library(patchwork)
library(ggpubr)
source(here("data_process/data_process.R"))

#create new theme
ms_theme <- function(){
  theme_tq() %+replace%
    theme(axis.title = element_text(size = 16), 
          axis.text = element_text(size = 14, color = "black"), 
          strip.text = element_text(size = 16, color = "white"), 
          legend.title = element_text(size = 16), 
          legend.text = element_text(size = 14), 
          title = element_text(size = 16, color = "#072025", face = "bold"))
}

col_pal <- c("#BFD1D9","#95B8BF", "#729CA6", "#4C6C73", "#2A3A40")
col_pal2 <- c("#ABBCBD", "#BFD1D9","#95B8BF", "#729CA6", "#4C6C73", "#2A3A40")


## Metabolics figure (Nazario et al., 2025 Figure 1) ####
#NOTE: silhouettes and arrow at 0min added in Adobe Illustrator 
#get resting values
resp_rest <- all_resp %>% filter(Trial == "Rest") %>% group_by(sp) %>% 
  summarise(mean_rest_o = mean(VO2_r, na.rm = TRUE), 
            sd_rest_o = sd(VO2_r, na.rm = TRUE),
            mean_rest_c = mean(VCO2_r), 
            sd_rest_c = sd(VCO2_r, na.rm = TRUE))

resp_rest_dl <- resp_rest %>% filter(sp == "Beluga")
resp_rest_tt <- resp_rest %>% filter(sp == "Dolphin")

sum_resp <- all_resp %>%
  group_by(sp, Trial) %>%
  summarise(perc_c_c = mean(VCO2, na.rm = T), 
            sd_perc_c_c = sd(VCO2, na.rm = T), 
            perc_o_c = mean(VO2, na.rm = T), 
            sd_perc_o_c = sd(VO2, na.rm = T), 
            perc_c_r = mean(VCO2_r, na.rm = T), 
            sd_perc_c_r = sd(VCO2_r, na.rm = T), 
            perc_o_r = mean(VO2_r, na.rm = T), 
            sd_perc_o_r = sd(VO2_r, na.rm = T), 
            mean_time_cost = 0, 
            mean_recover_c = mean(CO2_tobase), 
            sd_recover_c = sd(CO2_tobase),
            mean_recover_o = mean(O2_tobase), 
            sd_recover_o = sd(O2_tobase), 
            n = n())

#plot -- O2 
o2_rect <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf), 
                      ymin = c((resp_rest_dl$mean_rest_o-resp_rest_dl$sd_rest_o),
                               (resp_rest_tt$mean_rest_o-resp_rest_tt$sd_rest_o)), 
                      ymax = c((resp_rest_dl$mean_rest_o+resp_rest_dl$sd_rest_o),
                               (resp_rest_tt$mean_rest_o+resp_rest_tt$sd_rest_o)))
o2_line <- data.frame(yintercept = c(resp_rest_dl$mean_rest_o, resp_rest_tt$mean_rest_o), sp = c("Beluga", "Dolphin"))

o2_resp <- sum_resp %>% filter(Trial != "Rest") %>%
  ggplot(aes(color = Trial)) + 
  geom_rect(data = transform(o2_rect, sp = c("Beluga", "Dolphin")), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+  
  geom_hline(data = o2_line, aes(yintercept = yintercept), linewidth = 1, linetype = "dashed", color = "grey60")+
  geom_errorbar(aes(x = mean_time_cost, ymin = perc_o_c - sd_perc_o_c, ymax = perc_o_c + sd_perc_o_c, group = Trial), position = position_dodge(width = 2), size = 1, width = 0) +
  geom_point(aes(x = mean_time_cost, y = perc_o_c, group = factor(Trial)), position = position_dodge(width = 2), size = 4) + #plot cost points
  geom_errorbar(aes(x = mean_recover_o/60, ymin = perc_o_r - sd_perc_o_r, ymax = perc_o_r + sd_perc_o_r, group = Trial), position = position_dodge(width = 70), size = 1, width = 0) +
  geom_point(aes(x = mean_recover_o/60, y = perc_o_r), size = 4) + #plot rest points
  facet_wrap(~sp, scales = "free_y")+
  labs(x = " ", 
       y = expression(bold(paste(VO[2]["  "]('ml O'[2]%.%'kg'^'-1'%.%min^'-1')))), parse = TRUE)+
  xlim(-1, 10)+
  scale_color_manual(values = col_pal)+
  ms_theme()+
  theme(legend.position = "top",
        legend.justification = "left",
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())

#plot -- CO2
co2_rect <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf), 
                       ymin = c((resp_rest_dl$mean_rest_c-resp_rest_dl$sd_rest_c),
                                (resp_rest_tt$mean_rest_c-resp_rest_tt$sd_rest_c)),
                       ymax = c((resp_rest_dl$mean_rest_c+resp_rest_dl$sd_rest_c),
                                (resp_rest_tt$mean_rest_c+resp_rest_tt$sd_rest_c)))
co2_line <- data.frame(yintercept = c(resp_rest_dl$mean_rest_c, resp_rest_tt$mean_rest_c), sp = c("Beluga", "Dolphin"))

co2_resp <- sum_resp %>% filter(Trial != "Rest") %>%
  ggplot(aes(color = Trial)) + 
  geom_rect(data = transform(co2_rect, sp = c("Beluga", "Dolphin")), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(data = co2_line, aes(yintercept = yintercept), linewidth = 1, linetype = "dashed", color = "grey60")+
  geom_errorbar(aes(x = mean_time_cost, ymin = perc_c_c - sd_perc_c_c, ymax = perc_c_c + sd_perc_c_c, group = Trial), position = position_dodge(width = 2), size = 1, width = 0) +
  geom_point(aes(x = mean_time_cost, y = perc_c_c, group = factor(Trial)), position = position_dodge(width = 2), size = 4) + #plot cost points
  geom_errorbar(aes(x = mean_recover_c/60, ymin = perc_c_r - sd_perc_c_r, ymax = perc_c_r + sd_perc_c_r, group = Trial), position = position_dodge(width = 70), size = 1, width = 0) +
  geom_point(aes(x = mean_recover_c/60, y = perc_c_r), size = 4) + #plot rest points
  facet_wrap(~sp, scales = "free_y")+
  labs(x = expression(bold("Surface recovery time (min)")), 
       y = expression(bold(paste(VCO[2]["  "]('ml CO'[2]%.%'kg'^'-1'%.%min^'-1')))), parse = TRUE)+
  xlim(-1, 10)+
  scale_color_manual(values = col_pal)+
  ms_theme()+
  theme(legend.position = "none")

o2_resp/co2_resp

## Blood gas figure (Nazario et al., 2025 Figure 2) ####
#NOTE: plots combined and model y-intercept and rate of change text added in Adobe Illustrator
# get resting values 
rest_b_tt <- all_blood %>%
  filter(Trial == "Rest" & sp == "Dolphin") %>% 
  summarise(mean_lactate = mean(Lactate), 
            sd_lact = sd(Lactate),
            mean_O2 = mean(PO2), 
            sd_o2 = sd(PO2),
            mean_CO2 = mean(PCO2), 
            sd_CO2 = sd(PCO2),
            mean_pH = mean(pH), 
            sd_ph = sd(pH),
            mean_HCO3 = mean(HCO3), 
            sd_hco3 = sd(HCO3),
            mean_TCO2 = mean(TCO2), 
            sd_TCO2 = sd(TCO2),
            mean_sO2 = mean(sO2), 
            sd_so2 = sd(sO2))

rest_b_dl <- all_blood %>%
  filter(Trial == "Rest" & sp == "Beluga") %>%
  summarise(mean_lactate = mean(Lactate), 
            sd_lact = sd(Lactate),
            mean_O2 = mean(PO2), 
            sd_o2 = sd(PO2),
            mean_CO2 = mean(PCO2), 
            sd_CO2 = sd(PCO2),
            mean_pH = mean(pH), 
            sd_ph = sd(pH),
            mean_HCO3 = mean(HCO3), 
            sd_hco3 = sd(HCO3),
            mean_TCO2 = mean(TCO2), 
            sd_TCO2 = sd(TCO2),
            mean_sO2 = mean(sO2), 
            sd_so2 = sd(sO2))

#get 1 min bins means
all_blood <- all_blood %>%
  mutate(draw_time_min = Draw_Time/60,
         draw_time_min = round(draw_time_min, digits = 0))
all_blood_dl <- all_blood %>% filter(sp == "Beluga" & Trial != "Rest") %>%
  mutate(draw_time_min = Draw_Time/60, 
         draw_time_min = round(draw_time_min, digits = 0))
all_blood_tt <- all_blood %>% filter(sp == "Dolphin" & Trial != "Rest") %>%
  mutate(draw_time_min = Draw_Time/60, 
         draw_time_min = round(draw_time_min, digits = 0))

blood_sum <- all_blood %>%
  group_by(sp, Trial, draw_time_min) %>%
  summarise(mean_lact = mean(Lactate, na.rm = T), 
            mean_o2 = mean(PO2, na.rm = T), 
            sd_o2 = sd(PO2, na.rm = T),
            mean_co2 = mean(PCO2, na.rm = T), 
            mean_hco3 = mean(HCO3, na.rm =T), 
            mean_tco2 = mean(TCO2, na.rm = T),
            sd_tco2 = sd(TCO2, na.rm = T),
            mean_ph = mean(pH, na.rm = T), 
            sd_ph = sd(pH, na.rm = T),
            mean_so2 = mean(sO2, na.rm = T), 
            n = n())

#plot pO2
po2_rect_dl <- data.frame(xmin = c(-Inf), xmax = c(Inf), 
                          ymin = c((rest_b_dl$mean_O2-rest_b_dl$sd_o2)), 
                          ymax = c((rest_b_dl$mean_O2+rest_b_dl$sd_o2)))

po2_rect_tt <- data.frame(xmin = c(-Inf), xmax = c(Inf), 
                          ymin = c((rest_b_tt$mean_O2-rest_b_tt$sd_o2)), 
                          ymax = c((rest_b_tt$mean_O2+rest_b_tt$sd_o2)))

po2_r_dl <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Beluga") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_o2, fill = sp)) +
  geom_rect(data = po2_rect_dl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_dl$mean_O2, linewidth = 1, linetype = "dashed", color = "grey20")+
  geom_point(data = all_blood_dl, aes(x = draw_time_min, y = PO2, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_dl, aes(x = draw_time_min, y = PO2, group = Animal), alpha = 0.5, color = "grey")+
  geom_line(aes(color = sp), size = 1.5, alpha = 0.8) +
  #geom_errorbar(aes(ymin = mean_o2 - sd_o2, ymax = mean_o2+sd_o2))+
  geom_point(size = 5, shape = 21, color = "black", stroke= 1) +
  facet_grid(sp~Trial)+ 
  labs(x = "Surface recovery time (min)", 
       y = expression(bold(paste("pO"[2], " (mmHg)"))), parse = TRUE, 
       fill = "Species: ")+
  guides(color = FALSE)+
  scale_fill_manual(values = c("#D9D8D0"))+
  scale_color_manual(values = c("#D9D8D0"))+
  ms_theme()+
  ylim(20, 60)+
  theme(legend.position = "none")

po2_r_tt <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Dolphin") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_o2, fill = sp)) +
  geom_rect(data = po2_rect_tt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_tt$mean_O2, linewidth = 1, linetype = "dashed", color = "grey20")+
  geom_point(data = all_blood_tt, aes(x = draw_time_min, y = PO2, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_tt, aes(x = draw_time_min, y = PO2, group = Animal), alpha = 0.5, color = "grey")+
  geom_line(aes(color = sp), size = 1.5, alpha = 0.8) +
  #geom_errorbar(aes(ymin = mean_o2 - sd_o2, ymax = mean_o2+sd_o2))+
  geom_point(size = 5, shape = 21, color = "black", stroke= 1) +
  facet_grid(sp~Trial)+ 
  labs(x = "Surface recovery time (min)", 
       y = expression(bold(paste("pO"[2], " (mmHg)"))), parse = TRUE, 
       fill = "Species: ")+
  guides(color = FALSE)+
  scale_fill_manual(values = c("grey40"))+
  scale_color_manual(values = c("grey40"))+
  ms_theme()+
  #ylim(40, 85)+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 9, by = 3))+
  theme(legend.position = "none")

#plot pCO2
co2_rect_dl <- data.frame(xmin = c(-Inf), xmax = c(Inf),
                          ymin = c(rest_b_dl$mean_CO2-rest_b_dl$sd_CO2),
                          ymax = c((rest_b_dl$mean_CO2+rest_b_dl$sd_CO2)))

co2_rect_tt <- data.frame(xmin = c(-Inf), xmax = c(Inf),
                          ymin = c((rest_b_tt$mean_CO2-rest_b_tt$sd_CO2)),
                          ymax = c((rest_b_tt$mean_CO2+rest_b_tt$sd_CO2)))


co2_r_dl <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Beluga") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_co2, fill = sp)) +
  geom_rect(data = co2_rect_dl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_dl$mean_CO2, linewidth = 1, linetype = "dashed", color = "grey20")+
  geom_point(data = all_blood_dl, aes(x = draw_time_min, y = PCO2, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_dl, aes(x = draw_time_min, y = PCO2, group = Animal), alpha = 0.5, color = "grey")+
  geom_line(aes(color = sp), size = 1.5, alpha = 0.8) +
  #geom_errorbar(aes(ymin = mean_tco2 - sd_tco2, ymax = mean_tco2+sd_tco2))+
  geom_point(size = 5, shape = 21, color = "black", stroke = 1) +
  facet_grid(sp~Trial)+ 
  labs(x = "Surface recovery time (min)", 
       y = expression(bold(paste("pCO"[2], " (mmHg)"))), parse = TRUE)+
  scale_fill_manual(values = c("#D9D8D0"))+
  scale_color_manual(values = c("#D9D8D0"))+
  ms_theme()+
  ylim(50, 90)+
  theme(legend.position = "none")

co2_r_tt <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Dolphin") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_co2, fill = sp)) +
  geom_rect(data = co2_rect_tt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_tt$mean_CO2, linewidth = 1, linetype = "dashed", color = "grey20")+
  geom_point(data = all_blood_tt, aes(x = draw_time_min, y = PCO2, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_tt, aes(x = draw_time_min, y = PCO2, group = Animal), alpha = 0.5, color = "grey")+
  geom_line(aes(color = sp), size = 1.5, alpha = 0.8) +
  #geom_errorbar(aes(ymin = mean_tco2 - sd_tco2, ymax = mean_tco2+sd_tco2))+
  geom_point(size = 5, shape = 21, color = "black", stroke = 1) +
  facet_grid(sp~Trial)+ 
  labs(x = "Surface recovery time (min)", 
       y = expression(bold(paste("pCO"[2], " (mmHg)"))), parse = TRUE)+
  scale_fill_manual(values = c("grey40"))+
  scale_color_manual(values = c("grey40"))+
  ms_theme()+
  ylim(40, 80)+
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 9, by = 3))+
  theme(legend.position = "none")

all_o2 <- po2_r_dl/po2_r_tt
all_o2

all_co2 <- co2_r_dl/co2_r_tt
all_co2

## Blood ph figure (Nazario et al., 2025 Figure 3) ####
#NOTE: modeled y-intercept and rate of change text added in Adobe Illustrator

#get resting values
ph_rect_dl <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf),
                         ymin = c((rest_b_dl$mean_pH-rest_b_dl$sd_ph)),
                         ymax = c((rest_b_dl$mean_pH+rest_b_dl$sd_ph)))

ph_rect_tt <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf),
                         ymin = c((rest_b_tt$mean_pH-rest_b_tt$sd_ph)),
                         ymax = c((rest_b_tt$mean_pH+rest_b_tt$sd_ph)))

#plot 
ph_dl <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Beluga") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_ph, color = Trial)) +
  geom_rect(data = ph_rect_dl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_dl$mean_pH, linewidth = 1, linetype = "dashed", color = "grey60")+
  #geom_errorbar(aes(ymin = mean_ph - sd_ph, ymax = mean_ph+sd_ph), color = "grey40")+
  geom_point(data = all_blood_dl, aes(x = draw_time_min, y = pH, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_dl, aes(x = draw_time_min, y = pH, group = Animal), alpha = 0.5, color = "grey")+
  geom_point(aes(color = Trial), size = 4) +
  geom_line(aes(color = Trial), size = 1) +
  facet_wrap(~Trial, ncol = 3)+
  xlab("Surface recovery time (min)") +
  ylab("pH")+
  ylim(7.16, 7.3)+
  labs(title = "Beluga")+
  scale_color_manual(values = c(col_pal[1], col_pal[3], col_pal[5]))+
  ms_theme()

ph_tt <- blood_sum %>%
  filter(Trial != "Rest" & sp == "Dolphin") %>%
  ggplot(mapping = aes(x = draw_time_min, y = mean_ph, color = Trial)) +
  geom_rect(data = ph_rect_tt, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(yintercept = rest_b_tt$mean_pH, linewidth = 1, linetype = "dashed", color = "grey60")+
  #geom_errorbar(aes(ymin = mean_ph - sd_ph, ymax = mean_ph+sd_ph), color = "grey40")+
  geom_point(data = all_blood_tt, aes(x = draw_time_min, y = pH, shape = Animal), alpha = 0.8, color = "grey", size = 3)+
  geom_line(data = all_blood_tt, aes(x = draw_time_min, y = pH, group = Animal), alpha = 0.5, color = "grey")+
  geom_point(aes(color = Trial), size = 4) +
  geom_line(aes(color = Trial), size = 1) +
  facet_wrap(~Trial, ncol = 5)+
  xlab("Surface recovery time (min)") +
  ylab("pH")+
  ylim(7.22, 7.4)+
  labs(title = "Dolphin")+
  scale_x_continuous(breaks=seq(0,9,by=3))+
  scale_color_manual(values = c(col_pal[1], col_pal[3:5]))+
  ms_theme()

all_ph <- ggarrange(ph_dl, ph_tt, nrow = 2, ncol = 1, legend = "none")
all_ph

## Ventilation figure (Nazario et al., 2025 Figure 4) ####
dat_freq_rec$Trial <- factor(dat_freq_rec$Trial, levels = c("Rest","x1 Swim", "x3 Swim", "x5 Swim", "SAB"))

#beluga plots
#rest
bel_rest <- vent_binned %>%
  filter(Trial == "Rest" & sp == "Beluga") 

#rest freq 
bel_rest_hline_f = mean(bel_rest$breath_freq_min, na.rm = TRUE)
bel_rest_sd_f = sd(bel_rest$breath_freq_min, na.rm = TRUE)

#freq
bel_freq <- dat_freq_rec %>% filter(sp == "Beluga" & 
                                      elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = breath_freq_min, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = bel_rest_hline_f - bel_rest_sd_f, ymax = bel_rest_hline_f + bel_rest_sd_f), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = bel_rest_hline_f, color = "grey60", linetype = "dashed", size = 1) +
  geom_point(alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE, linewidth = 2)+
  #facet_wrap(~Trial_Type)+
  #xlim(1, 11)+
  ylim(1,14)+
  xlab("Surface recovery time (min)")+
  ylab("Breath frequency per 60 sec")+
  labs(title = "Beluga ventilation")+
  ms_theme()+
  scale_color_manual(values = c(col_pal[3], col_pal[5]))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none", 
        legend.justification = "left")

#dur
bel_rest_hline_d = mean(bel_rest$mean_breath_dur, na.rm = TRUE)
bel_rest_sd_d = sd(bel_rest$mean_breath_dur, na.rm = TRUE)

bel_dur <- dat_freq_rec %>% filter(sp == "Beluga" & 
                                     elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = mean_breath_dur, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = bel_rest_hline_d - bel_rest_sd_d, ymax = bel_rest_hline_d + bel_rest_sd_d), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = bel_rest_hline_d, color = "grey60", linetype = "dashed", size = 1) +
  geom_point(alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE, linewidth = 2)+
  #facet_wrap(~Trial_Type)+
  #xlim(1, 11)+
  ylim(0,4)+
  xlab("Surface recovery time (min)")+
  ylab("Breath duration (s)")+
  ms_theme()+
  scale_color_manual(values = c(col_pal[3], col_pal[5]))+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none")

#IBI
bel_rest_hline_i = mean(bel_rest$mean_ibi, na.rm = TRUE)
bel_rest_sd_i = sd(bel_rest$mean_ibi, na.rm = TRUE)

bel_ibi <- dat_freq_rec %>% filter(sp == "Beluga" & 
                                     elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = mean_ibi, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = bel_rest_hline_i - bel_rest_sd_i, ymax = bel_rest_hline_i + bel_rest_sd_i), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = bel_rest_hline_i, color = "grey60", linetype = "dashed", size = 1) +
  geom_point( alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE,  linewidth = 2)+
  # facet_wrap(~Trial_Type)+
  #xlim(1, 11)+
  ylim(0,100)+
  xlab("Surface recovery time (min)")+
  ylab("Inter-breath Interval (s)")+
  ms_theme()+
  scale_color_manual(values = c(col_pal[3], col_pal[5]))+
  theme(legend.position = "none")

#dolphin plots
#rest
dol_rest <- vent_binned %>%
  filter(Trial == "Rest" & sp == "Dolphin") 

#freq
dol_rest_hline_f = mean(dol_rest$breath_freq_min, na.rm = TRUE)
dol_rest_sd_f = sd(dol_rest$breath_freq_min, na.rm = TRUE)

dol_freq <- dat_freq_rec %>% filter(sp == "Dolphin" & 
                                      elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = breath_freq_min, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = dol_rest_hline_f - dol_rest_sd_f, ymax = dol_rest_hline_f + dol_rest_sd_f), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = dol_rest_hline_f, color = "grey60", linetype = "dashed", size = 1) +
  geom_point(alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE, linewidth = 2)+
  #facet_wrap(~Trial_Type, ncol = 4)+
  #xlim(1, 11)+
  ylim(1,14)+
  xlab("Surface recovery time (min)")+
  ylab("Breath frequency per 60 sec")+
  labs(title = "Dolphin ventilation", color = "Trial type")+
  ms_theme()+
  scale_color_manual(values = col_pal[2:5])+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none")
dol_freq2 <- dol_freq + theme(legend.position = "top", legend.justification = "left")

#dur
dol_rest_hline_d = mean(dol_rest$mean_breath_dur, na.rm = TRUE)
dol_rest_sd_d = sd(dol_rest$mean_breath_dur, na.rm = TRUE)

dol_dur <- dat_freq_rec %>% filter(sp == "Dolphin" & 
                                     elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = mean_breath_dur, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = dol_rest_hline_d - dol_rest_sd_d, ymax = dol_rest_hline_d + dol_rest_sd_d), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = dol_rest_hline_d, color = "grey60", linetype = "dashed", size = 1) +
  geom_point(alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE, linewidth = 2)+
  #facet_wrap(~Trial_Type, ncol = 4)+
  #xlim(1, 11)+
  ylim(0,4)+
  xlab("Surface recovery time (min)")+
  ylab("Breath duration (s)")+
  ms_theme()+
  scale_color_manual(values = col_pal[2:5])+
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none")

#ibi
dol_rest_hline_i = mean(dol_rest$mean_ibi, na.rm = TRUE)
dol_rest_sd_i = sd(dol_rest$mean_ibi, na.rm = TRUE)

dol_ibi <- dat_freq_rec %>% filter(sp == "Dolphin" & 
                                     elapsed_time < 10) %>%
  ggplot(aes(x = elapsed_time, y = mean_ibi, color = Trial))+ 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = dol_rest_hline_i - dol_rest_sd_i, ymax = dol_rest_hline_i + dol_rest_sd_i), fill = "grey90", color = "grey90") +
  geom_hline(yintercept = dol_rest_hline_i, color = "grey60", linetype = "dashed", size = 1) +
  geom_point(alpha = 0.2, size = 1.5)+
  geom_smooth(se = FALSE, linewidth = 2)+
  #facet_wrap(~Trial_Type, ncol = 4)+
  #xlim(1, 11)+
  ylim(0,100)+
  xlab("Surface recovery time (min)")+
  ylab("Inter-breath Interval (s)")+
  ms_theme()+
  scale_color_manual(values = col_pal[2:5])+
  theme(legend.position = "none")

#combine
bel_vent_plots <- bel_freq/bel_dur/bel_ibi
dol_vent_plots <- dol_freq/dol_dur/dol_ibi

all_vent_plot <- ggarrange(bel_vent_plots, dol_vent_plots, ncol = 2, common.legend = TRUE, legend.grob = get_legend(dol_freq2))
all_vent_plot

## Vasculature figure (Nazario et al., 2025 Figure 5) ####
#format data
dat_f <- dat_f %>% mutate(sp = as.factor(sp))
levels(dat_f$sp) <- c("Beluga", "Dolphin")

#rest hlines by species
#beluga
bel_rest_vasc <- dat_f %>% 
  filter(Trial == "Rest" & sp == "Beluga") %>% 
  summarise(mean_mean = mean(surf_mean, na.rm = TRUE),
            sd_mean = sd(surf_mean, na.rm = TRUE), 
            mean_max = mean(max_max, na.rm = TRUE), 
            sd_max = sd(max_max, na.rm = TRUE), 
            mean_vaso = mean(prop_vasc, na.rm = TRUE), 
            sd_vaso = sd(prop_vasc, na.rm = TRUE))

#dolphin
dol_rest_vasc <- dat_f %>% 
  filter(Trial == "Rest" & sp == "Dolphin") %>% 
  summarise(mean_mean = mean(surf_mean, na.rm = TRUE),
            sd_mean = sd(surf_mean, na.rm = TRUE), 
            mean_max = mean(max_max, na.rm = TRUE), 
            sd_max = sd(max_max, na.rm = TRUE), 
            mean_vaso = mean(prop_vasc, na.rm = TRUE), 
            sd_vaso = sd(prop_vasc, na.rm = TRUE))

#rest rectangles
mean_rect <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf), 
                        ymin = c((bel_rest_vasc$mean_mean-bel_rest_vasc$sd_mean),
                                 (dol_rest_vasc$mean_mean-dol_rest_vasc$sd_mean)), 
                        ymax = c((bel_rest_vasc$mean_mean+bel_rest_vasc$sd_mean),
                                 (dol_rest_vasc$mean_mean+dol_rest_vasc$sd_mean)))
mean_line <- data.frame(yintercept = c(bel_rest_vasc$mean_mean, dol_rest_vasc$mean_mean),
                        sp = c("Beluga", "Dolphin"))

max_rect <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf), 
                       ymin = c((bel_rest_vasc$mean_max-bel_rest_vasc$sd_max),
                                (dol_rest_vasc$mean_max-dol_rest_vasc$sd_max)), 
                       ymax = c((bel_rest_vasc$mean_max+bel_rest_vasc$sd_max),
                                (dol_rest_vasc$mean_max+dol_rest_vasc$sd_max)))
max_line <- data.frame(yintercept = c(bel_rest_vasc$mean_max, dol_rest_vasc$mean_max), 
                       sp = c("Beluga", "Dolphin"))

vaso_rect <- data.frame(xmin = c(-Inf, -Inf), xmax = c(Inf, Inf), 
                        ymin = c((bel_rest_vasc$mean_vaso-bel_rest_vasc$sd_vaso),
                                 (dol_rest_vasc$mean_vaso-dol_rest_vasc$sd_vaso)), 
                        ymax = c((bel_rest_vasc$mean_vaso+bel_rest_vasc$sd_vaso),
                                 (dol_rest_vasc$mean_vaso+dol_rest_vasc$sd_vaso)))
vaso_line <- data.frame(yintercept = c(bel_rest_vasc$mean_vaso, dol_rest_vasc$mean_vaso),
                        sp = c("Beluga", "Dolphin"))

#values summarized into 3min bins
dat_f <- dat_f %>%
  mutate(time_min = img_sec %/% 180)

sum_dat_f <- dat_f %>%
  group_by(sp, Trial, time_min) %>%
  summarise(mean_mean = mean(surf_mean, na.rm = T), 
            mean_max = mean(max_max, na.rm = T), 
            mean_vaso = mean(prop_vasc, na.rm = T))

#mean surface plot
mean_surf_p <- sum_dat_f %>% 
  filter(Trial != "Rest") %>%
  ggplot(aes(time_min*3, mean_mean)) + #multiply by three to translate time interval grouping done above
  geom_rect(data = transform(mean_rect, sp = c("Beluga", "Dolphin")), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(data = mean_line, aes(yintercept = yintercept), color = "grey60", linetype = "dashed", size = 1) + 
  geom_point(size = 4, aes(fill = Trial, color = Trial), pch = 21) + 
  geom_line(linewidth = 1, aes(color = Trial)) + 
  facet_wrap(~sp, scales = "free_x") + 
  scale_fill_manual(values = c(col_pal[1], col_pal[3:5]))+
  scale_color_manual(values = c(col_pal[1], col_pal[3:5]))+
  xlab("Surface recovery time (min)")+
  ylab("Mean temperature (°C)")+
  scale_x_continuous(breaks=seq(0,20,by=3))+
  ms_theme()+
  theme(legend.position = "top", 
        legend.justification = "left", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())

#max surface plot
max_surf_p <- sum_dat_f %>% 
  filter(Trial != "Rest") %>%
  ggplot(aes(time_min*3, mean_max)) + #multiply by three to translate time interval grouping done above
  geom_rect(data = transform(max_rect, sp = c("Beluga", "Dolphin")), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(data = max_line, aes(yintercept = yintercept), color = "grey60", linetype = "dashed", size = 1) + 
  geom_point(size = 4, aes(fill = Trial, color = Trial), pch = 21) + 
  geom_line(linewidth = 1, aes(color = Trial)) + 
  facet_wrap(~sp, scales = "free_x") + 
  scale_fill_manual(values = c(col_pal[1], col_pal[3:5]))+
  scale_color_manual(values = c(col_pal[1], col_pal[3:5]))+
  xlab("Surface recovery time (min)")+
  ylab("Max temperature (°C)")+
  scale_x_continuous(breaks=seq(0,20,by=3))+
  ms_theme()+
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())

#vaso plot
vaso_p <- sum_dat_f %>% 
  filter(Trial != "Rest") %>%
  ggplot(aes(time_min*3, mean_vaso*100)) + #multiply by three to translate time interval grouping done above
  geom_rect(data = transform(vaso_rect, sp = c("Beluga", "Dolphin")), aes(xmin = xmin, xmax = xmax, ymin = ymin*100, ymax = ymax*100), fill = "grey90", color = "grey90", inherit.aes = FALSE)+
  geom_hline(data = vaso_line, aes(yintercept = yintercept*100), color = "grey60", linetype = "dashed", size = 1) + 
  geom_point(size = 4, aes(fill = Trial, color = Trial), pch = 21) + 
  geom_line(linewidth = 1, aes(color = Trial)) + 
  facet_wrap(~sp, scales = "free_x") + 
  scale_fill_manual(values = c(col_pal[1], col_pal[3:5]))+
  scale_color_manual(values = c(col_pal[1], col_pal[3:5]))+
  xlab("Surface recovery time (min)")+
  ylab("Percent perfusion (%)")+
  scale_x_continuous(breaks=seq(0,20,by=3))+
  ms_theme()+
  theme(legend.position = "none")

#combined
all_vasc_plot <- mean_surf_p/max_surf_p/vaso_p
all_vasc_plot

# Supplemental Tables ####
### Metabolics table #####
#NOTE: x2 swims are included here and not in the analysis.
resp_tab <- all_resp %>%
  group_by(sp, Trial) %>%
  summarise(`Mean rest VO2` = mean(VO2_r, na.rm = TRUE), 
            `SD rest VO2` = sd(VO2_r, na.rm = TRUE), 
            `Mean rest VCO2` = mean(VCO2_r, na.rm = TRUE), 
            `SD rest VCO2` = sd(VCO2_r, na.rm = TRUE),
            `Mean swim VO2` = mean(VO2, na.rm = TRUE), 
            `SD swim VO2` = sd(VO2, na.rm = TRUE), 
            `Mean swim VCO2` = mean(VCO2, na.rm = TRUE), 
            `SD swim VCO2` = sd(VCO2, na.rm = TRUE), 
            `Mean O2 recovery time` = mean(O2_tobase/60, na.rm = TRUE), 
            `SD O2 recovery time` = sd(O2_tobase/60, na.rm = TRUE), 
            `Mean CO2 recovery time` = mean(CO2_tobase/60, na.rm = TRUE), 
            `SD CO2 recovery time` = sd(CO2_tobase/60, na.rm = TRUE), 
            n = n())

kableExtra::kbl(resp_tab) %>%
  kable_paper("striped") %>%
  add_header_above(c(" ", " " = 1, "Resting rates (ml/kg*min)" = 4, "Exercise cost rates (ml/kg*min)" = 4, "Return to Rest (min)" = 4, " " = 1)) 

### Blood table ####
blood_tab <- all_blood %>%
  group_by(sp, Trial) %>%
  summarise(`Mean lactate` = mean(Lactate), 
            `SD lactate` = sd(Lactate),
            `Mean pO2` = mean(PO2), 
            `SD pO2` = sd(PO2),
            `Mean pCO2` = mean(PCO2), 
            `SD pCO2` = sd(PCO2),
            `Mean pH` = mean(pH), 
            `SD ph` = sd(pH),
            `Mean HCO3` = mean(HCO3), 
            `SD Hco3` = sd(HCO3),
            `Mean TCO2` = mean(TCO2), 
            `SD TCO2` = sd(TCO2),
            `Mean sO2` = mean(sO2), 
            `SD sO2` = sd(sO2), 
            n = n()) %>%
  mutate_if(is.numeric, ~ round(.,3))

kableExtra::kbl(blood_tab, digits = 3) %>%
  kable_paper("striped")

### Ventilation table ####
#NOTE: x2 swims are included here and not in the analysis.
dolph_vent2 <- dolph_vent %>%
  mutate(Elapsed = Start_Time - head(Start_Time, n = 1, na.rm = TRUE),
         Recover_per = ifelse(Trial_Type == "Rest", "Resting", 
                              ifelse(Trial_Type == "x1 Swim" & Elapsed > 150, "Resting", 
                                     ifelse(Trial_Type == "x2 Swim" & Elapsed > 150, "Resting", 
                                            ifelse(Trial_Type == "x3 Swim" & Elapsed > 150, "Resting",
                                                   ifelse(Trial_Type == "x5 Swim" & Elapsed > 150, "Resting",
                                                          ifelse(Trial_Type == "SAB" & Elapsed > 150, "Resting", "Dive Cost")))))), 
         sp = "Dolphin")
bel_vent2 <- bel_vent %>%
  mutate(Elapsed = Start_Time - head(Start_Time, n = 1, na.rm = TRUE), 
         Recover_per= ifelse(Trial_Type == "Rest", "Resting", 
                             ifelse(Trial_Type == "x3 Swim" & Elapsed > 150, "Resting", 
                                    ifelse(Trial_Type == "SAB" & Elapsed > 150, "Resting", "Dive cost"))), 
         sp = "Beluga")

swim_vent <- rbind(dolph_vent2, bel_vent2)
swim_vent$Animal <- factor(swim_vent$Animal, levels = c("Maple", "Qinu", "Donley", "Rain"))

swim_vent <- swim_vent %>% group_by(id, Recover_per) %>% mutate(freq = (max(Breath_Number)/max(Elapsed, na.rm= TRUE))*60) %>% ungroup()
swim_vent_tt <- swim_vent %>% filter(sp == "Dolphin")
swim_vent_tt$Trial_Type <- factor(swim_vent_tt$Trial_Type, levels = c("Rest", "x1 Swim", "x2 Swim", "x3 Swim", "x5 Swim", "SAB"))

swim_vent_dl <- swim_vent %>% filter(sp == "Beluga")
swim_vent_dl$Trial_Type <- factor(swim_vent_dl$Trial_Type, levels = c("Rest", "x3 Swim", "SAB"))

vent_tab <- rbind(swim_vent_tt, swim_vent_dl)

vent_tab <- swim_vent %>%
  group_by(sp, Trial_Type) %>%
  summarise(`Mean breath frequency (breaths/min)` = mean(freq, na.rm = TRUE), 
            `SD breath freqeuncy (breaths/min` = sd(freq, na.rm = TRUE), 
            `Mean breath duration (s)` = mean(Breath_Duration, na.rm = TRUE), 
            `SD  breath duration (s)` = sd(Breath_Duration, na.rm = TRUE), 
            `Mean IBI (s)` = mean(IBI, na.rm = TRUE), 
            `SD IBI (s)` = sd(IBI, na.rm = TRUE),
            n = length(unique(id))) %>%
  mutate_if(is.numeric, ~ round(.,3))

kableExtra::kbl(vent_tab) %>%
  kable_paper("striped")

### Vasculature table ####
vaso_tab <- dat_f %>% 
  group_by(sp, Trial) %>% 
  summarise(`Mean surface temperature (°C)` = mean(surf_mean, na.rm = TRUE),
            `SD surface temperature (°C)` = sd(surf_mean, na.rm = TRUE), 
            `Mean maximum surface temperature (°C)` = mean(max_max, na.rm = TRUE),
            `SD maximum surface temperature (°C)` = sd(max_max, na.rm = TRUE),
            `Mean proportion vasodilation` = mean(prop_vasc, na.rm = TRUE)*100, 
            `SD proportion vasodilation` = sd(prop_vasc, na.rm = TRUE)*100, 
            n = length(unique(id))) %>%
  mutate_if(is.numeric, ~ round(.,3))

kableExtra::kbl(vaso_tab) %>%
  kable_paper("striped")


### Speed figure ####
speed_l <- gather(all_speed, speed_num, speed, Speed_S1:Speed_S5, factor_key = TRUE) %>%
  filter(Trial != "Rest") %>%
  subset(select = -c(Bow_rate)) %>%
  mutate(speed_num = gsub("\\D", "", speed_num))

sum <- speed_l %>%
  filter(Trial != "SAB") %>%
  group_by(sp, Trial, speed_num) %>%
  summarise(mean = mean(as.numeric(speed), na.rm = TRUE), 
            sd = sd(as.numeric(speed), na.rm = TRUE)) 

ggplot(sum, mapping = aes(x = speed_num, y = mean)) +
  geom_bar(stat = "identity", position = position_dodge(1), aes(fill = Trial), width = 1.1, alpha = 0.9) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, group = Trial), width = 0.2, position = position_dodge(1))+
  facet_wrap(~sp, nrow = 1, scales = "free_x") +
  xlab("Swim Number") +
  labs(fill = "Trial", 
       y = expression(bold(paste("Mean speed "('m'%.%'s'^'-1'))))) +
  scale_fill_manual(values = col_pal)+
  ms_theme()+
  theme(axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        strip.text = element_text(size = 12), 
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12))




