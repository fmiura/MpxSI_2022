############################################################
#Aim: to visualize MPX serial interval data
#Final edit: 22 November 2022
#Editor: Fumi Miura & Jantien Backer
############################################################
###Procedure
#0. Package 
#1. Import data 
#2. ggplot
############################################################

###0. Package -----
#library(rstan)
#library(ggforce)
#library(magrittr)
library(patchwork)
library(tidyverse)
library(glue)
library(ggridges)
library(cowplot)

###1. Import data -----
#raw data
Anonym_All_data_score <- read_csv("Anonym_All_data_score.csv")
Anonym_gantt_data <- read_rds("Anonym_gantt_data.rds")
pos_gamma_df <- read_rds("pos_gamma_df.rds")
summary_gamma_df <- read_rds("summary_gamma_df.rds")
pos_normal_df <- read_rds("pos_normal_df.rds")
summary_normal_df <- read_rds("summary_normal_df.rds")
#reorder factor
pos_normal_df <- pos_normal_df %>% mutate(GGD_index = fct_relevel(GGD_index, "All pooled","1", "2", "3", "4", "5", "6", "7", "8", "9"))
pos_gamma_df <- pos_gamma_df %>% mutate(GGD_index = fct_relevel(GGD_index, "All pooled","1", "2", "3", "4", "5", "6", "7", "8", "9"))
#Fig-1A: histogram of reported serial interval values -----
# two pretty color schemes
bivariate_color_scale_1 <- tibble(
  "9" = "#3F2949", 
  "8" = "#77324C",
  "7" = "#AE3A4E", 
  "6" = "#435786",
  "5" = "#806A8A", 
  "4" = "#BC7C8F",
  "3" = "#4885C1", 
  "2" = "#89A1C8",
  "1" = "#CABED0" 
)

bivariate_color_scale_2 <- tibble(
  "9" = "#3b4994", 
  "8" = "#8c62aa",
  "7" = "#be64ac", 
  "6" = "#5698b9",
  "5" = "#a5add3", 
  "4" = "#dfb0d6",
  "3" = "#5ac8c8", 
  "2" = "#ace4e4",
  "1" = "#e8e8e8" 
)

# main histogram 
main_plot <- ggplot(data = Anonym_All_data_score,
                    aes(x = SI, fill = as.character(score)))+
  geom_bar(alpha = 1, 
           width = 1) +
  scale_fill_manual(values = bivariate_color_scale_2) +
  labs(x="Time between symptom onsets (days)", y ="Count", fill = "Reliability") +
  theme_bw() +
  guides(fill = "none") +
  scale_x_continuous(breaks = seq(-20, 30 , by=2)) +
  scale_y_continuous(breaks = seq(0, 15, by=1), 
                     minor_breaks = 0:15,
                     limits = c(0, 14),
                     expand = c(0, 0))
# legend panel
legend <- expand_grid(x = 1:3,
                      y = 1:3) %>% 
  mutate(score = 3*(y-1) + x) %>% 
  ggplot(aes(x, y, fill = as.character(score))) + 
  geom_tile() +
  coord_equal() +
  scale_fill_manual(values = bivariate_color_scale_2) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("unreliable\nsymptom onset", "plausible\nsymptom onset", "reliable\nsymptom onset")) +
  scale_y_continuous(breaks = 1:3,
                     labels = c("unlikely case pair", "likely case pair,\nmultiple contacts", "likely case pair,\nsingle contact")) +
  labs(#title = "Reliability",
    # x = "symptom onset",
    # y = "case pair") +
    x = NULL,
    y = NULL) +
  # geom_text(x = -0.2, y = 3.7, label = "case pair", hjust = 0.5) +
  # geom_text(x = 3.7, y = -0.2, label = "symptom onset", vjust = 0.5, angle = -90) +
  guides(fill = "none") +
  theme_minimal() +
  theme(axis.title.y.left = element_text(angle = 0,
                                         size = 9,
                                         hjust = 1,
                                         #color = "#ffc000",
                                         margin = margin(r = -60, l = 10)),
        axis.title.x.bottom = element_text(angle = -90,
                                           size = 9,
                                           hjust = 1,
                                           vjust = 1,
                                           margin = margin(t = -30, b = 10)),
        axis.text.x.bottom = element_text(angle = -90,
                                          #size = 12,
                                          hjust = 0,
                                          vjust = 0.5),
        #axis.text.y.left = element_text(size = 12),
        panel.grid = element_blank())

Fig_1A <- ggdraw(main_plot) +
  draw_plot(legend, x = 0.5, y = 0.4, width = 0.6, height = 0.6)

ggsave("Fig_1A.tiff", plot=Fig_1A, units="cm", width=15, height=9, dpi=300)

#Fig-1B: transmission pairs in Amsterdam -----
##using time from the last moment of exposure
Fig_1B <- ggplot(data=Anonym_gantt_data) +
  geom_linerange(mapping = aes(x=pair_ID, ymin=infector_t, ymax=infectee_t, color=Presymp_status), size=0.5)+
  geom_point(mapping = aes(x=pair_ID, y=infector_t, color=Presymp_status), size=2)+
  geom_point(mapping = aes(x=pair_ID, y=infectee_t, color=Presymp_status), size=2, shape=17)+
  geom_point(mapping = aes(x=pair_ID, y=exactEX_t, color=Presymp_status), size=3, shape=4)+
  geom_linerange(mapping = aes(x=pair_ID, ymin=startEX_t, ymax=endEX_t, color=Presymp_status), size=2, alpha=0.5)+
  scale_colour_manual(values=c(RColorBrewer::brewer.pal(4, "RdBu")[1],RColorBrewer::brewer.pal(4, "RdBu")[4],RColorBrewer::brewer.pal(4, "RdBu")[11]))+
  coord_flip()+
  theme_bw()+
  labs(y="Time since the end of exposure (days)", x ="Pair ID", col = "Pre-symptomatic")

ggsave("Fig_1B.tiff", plot=Fig_1B, units="cm", width=15, height=9, dpi=300)

#Fig-1(A)(B)
Fig_1AB <- (Fig_1A/Fig_1B) + plot_annotation(tag_levels = "a")
ggsave("Fig_1AB.tiff", plot=Fig_1AB, units="cm", width=15, height=14, dpi=300)

#Fig-2: pooled estimates of serial intervals by GGD -----
##normal 
Fig2_poolSI_norm <- ggplot() +
  geom_density_ridges(data=pos_normal_df, aes(y = GGD_index, x = pos), col = NA, scale = 1,alpha = 0.4, fill = "blue")+
  geom_errorbarh(data=summary_normal_df,aes(y=GGD_index, x=mean, xmin=lower, xmax=upper, height = 0), size=0.7)+
  geom_point(data=summary_normal_df,aes(y=GGD_index, x=mean), size=2)+
  geom_vline(xintercept=summary_normal_df$mean[1], lty=2, size=0.3) +  # add a dotted line at x = (the pooled mean) after flip
  theme_bw()  +
  xlim(c(0,35))+
  ylab("Study region") + xlab("Mean serial interval (95% CrI)") +
  geom_text(data = mutate_if(summary_normal_df, 
                             is.numeric, round, 1),
            aes(label = glue("(N={n}) {mean} [{lower}, {upper}]"), x = 33, y = (1:10)+0.15)
  )
ggsave("Fig2_poolSI_norm.tiff", plot=Fig2_poolSI_norm, units="cm", width=5*5, height=3*5, dpi=300)
##gamma
FigS2_poolSI_gamma <- ggplot() +
  geom_density_ridges(data=pos_gamma_df, aes(y = GGD_index, x = pos), col = NA, scale = 1,alpha = 0.4, fill = "blue")+
  geom_errorbarh(data=summary_gamma_df,aes(y=GGD_index, x=mean, xmin=lower, xmax=upper, height = 0), size=0.7)+
  geom_point(data=summary_gamma_df,aes(y=GGD_index, x=mean), size=2)+
  geom_vline(xintercept=summary_gamma_df$mean[1], lty=2, size=0.3) +  # add a dotted line at x = (the pooled mean) after flip
  theme_bw()  +
  xlim(c(0,35))+
  ylab("Study region") + xlab("Mean serial interval (95% CrI)") +
  geom_text(data = mutate_if(summary_gamma_df, 
                             is.numeric, round, 1),
            aes(label = glue("(N={n}) {mean} [{lower}, {upper}]"), x = 33, y = (1:10)+0.15)
  )