################################################################################
# Purpose: model used in analysis for
#          Submitted manuscript [Nature Food]
# NOTE   : this code is a guide for transparency and
#          reproducibility
################################################################################ 

################################################################################
# quantile regression on shrimp microplastic (MP)
################################################################################

# Key libraries used in this model
library(quantreg)  # quantile regressio analysis
library(olsrr)   # OLSR model
library(tidyverse)  # ggplot and data wrangling
library(dplyr)  # Data organizing
library(ggpubr) # organising plots
library(strucchange) #
library(changepoint) # change point analysis
library(readxl)
library(quantregForest)
library(qgam)
library(mgcViz)
library(ellipse) # correlation plot
library(RColorBrewer) # corelation plot
library(corrplot)   
library(ggplot2)
library(ggtext)
library(grid) # For proper margin units
library(gratia) # diagnostic plots
library(patchwork) # diagnostics
library(changepoint)

# loading data 
imp.dat<-read_excel("data/shrimp_data.xlsx",sheet = "shrimp")
summary(imp.dat)

################################################################################ 
# Length
################################################################################
# fit ols model-----------------------------------------------------------------
ols_l <- lm(MPperind ~ AvgL, data= imp.dat)
ols_l_pred <- predict(ols_l, newdata = imp.dat)
summary(ols_l)
AIC(ols_l)

# fit linear QR model-----------------------------------------------------------
fit_l_99 <- rq(MPperind ~ AvgL, tau = 0.99, data = imp.dat)
summary (fit_l_99,se="boot")  
AIC(fit_l_99)
# generate predicted values
pred_l_99 <- predict(fit_l_99, newdata = data.frame
                     (AvgL = seq(min(imp.dat$AvgL),
                                 max(imp.dat$AvgL), 
                                 length.out = 10000)), 
                     interval = "confidence")
# Merge the predicted values with the original data
pred_l_data99 <- data.frame(AvgL = seq(min(imp.dat$AvgL),
                                       max(imp.dat$AvgL),
                                       length.out = 10000), 
                            MPperind = pred_l_99[, 1], lower = pred_l_99[, 2],
                            upper = pred_l_99[, 3])

# setting parameters for deriving equation from the quantile summary and coefs.
coefs_l <- coef(fit_l_99)
intercept_l <- round(coefs_l[1], digit = 2) # intercept
a_l <- round(coefs_l[2], digit = 2) # Second coef.
b_l <- round(coefs_l[3], digit = 3) # Third coef.

# fit non-linear QR model-------------------------------------------------------
fit_nl_99 <- rq(MPperind ~ AvgL + I(AvgL^2), tau = 0.99, data = imp.dat)
summary (fit_nl_99,se="boot")  
AIC(fit_nl_99)
summary(fit_nl_99)
# generate predicted values
pred_nl_99 <- predict(fit_nl_99, newdata = data.frame
                      (AvgL = seq(min(imp.dat$AvgL),
                                  max(imp.dat$AvgL), 
                                  length.out = 10000)), 
                      interval = "confidence")
# Merge the predicted values with the original data
pred_nl_data99 <- data.frame(AvgL = seq(min(imp.dat$AvgL),
                                        max(imp.dat$AvgL),
                                        length.out = 10000), 
                             MPperind = pred_nl_99[, 1], 
                             lower = pred_nl_99[, 2], upper = pred_nl_99[, 3])

# setting parameters for deriving equation from the quantile summary and coefs.
coefs_nl <- coef(fit_nl_99)
intercept_nl <- round(coefs_nl[1], digit = 2) # intercept
a_nl <- round(coefs_nl[2], digit = 2) # Second coef.
b_nl <- round(coefs_nl[3], digit = 3) # Third coef.

# fit ols model-----------------------------------------------------------------
ols_nl <- lm(MPperind ~ AvgL + I(AvgL^2), data= imp.dat)
ols_nl_pred <- predict(ols_nl, newdata = imp.dat)
summary(ols_nl)
AIC(ols_nl)

# qgam model--------------------------------------------------------------------
# fitting model at 0.90 quantile for shrimp length -----------------------------
gam_model <- gam(
  MPperind ~
    s(AvgL),
  data = imp.dat
)


fitL1 <- qgam(list(form = MPperind ~ s(AvgL, k = 20, bs = "ts"), ~ s(AvgL)),
              data = imp.dat, qu = 0.99)
fitlv <- getViz(fitL1)
fitgam <- getViz(gam_model)
summary(fitgam)


################################################################################
# change point analysis ########################################################
################################################################################

change_point_1 <- summary(ols_l)$coefficients[1,2]
change_point_2 <- summary(fit_l_99)$coefficients[1,2]
change_point_3 <- summary(fit_nl_99)$coefficients[1,2]
change_point_4 <- summary(fitlv)$coefficients[1,2]
change_point_5 <- summary(ols_nl)$coefficients[1,2]

fitL1$coefficients[1]
fitlv$formula
ols_l$coefficients


################################################################################
AIC(ols_l, fit_l_99, fit_nl_99)
################################################################################
# R1 measure
# pseudo-R^2 measure suggested by Koenker and Machado (1999)
# R1(τ) should lie in [0,1]
# where 1 would correspond to a perfect fit since the numerator which consists 
# of the weighted sum of deviations would be zero.
################################################################################
fitimp1 <- rq(MPperind_log ~ AvgL_log, tau =0.99, data = imp.dat)
fitimp0 <- rq(MPperind_log~1, tau =0.99, data = imp.dat)
rho <- function(u,tau=.5)u*(tau - (u < 0))
R1imp <- 1 - fitimp1$rho/fitimp0$rho

fitW1 <- rq(MPperind ~ AvgL + I(AvgL^2), tau =0.99, data = imp.dat)
fitW0 <- rq(MPperind ~ 1, tau =0.99, data = imp.dat)
R1_W <- 1 - fitW1$rho/fitW0$rho


R1imp
R1_W

library(Qtools)

GOFTest(fit_l_99)
GOFTest(fit_nl_99)

################################################################################
# combined model using qgam approach at 0.99 quantile
################################################################################
Q_qgam01 <- qgam(MPperind ~ 
                   s(Longitude, Latitude, bs="ts") +
                   s(AvgL, bs = "ts", k = 8) +
                   s(AvgW, bs= "ts", k = 8) +
                   s(SpeciesNo, bs ="ad", k = 8) +
                   as.factor(Continent)+
                   te(Years,bs="ts", d =1),
                 argGam = list(select=F),
                 data = imp.dat, qu=0.99)

Q_qgam02 <- qgam(MPperind ~ 
                   s(Longitude, Latitude, bs="ts") +
                   s(AvgL, bs = "ts", k = 8) +
                   s(AvgW, bs= "ts", k = 8) +
                   s(SpeciesNo, bs ="ad", k = 8) +
                   as.factor(Continent)+
                   te(Years,bs="ts", d =1),
                 argGam = list(select=F),
                 data = imp.dat, qu=0.5)

Q_qgam02 <- gam(MPperind ~ 
                  s(Longitude, Latitude, bs="ts") +
                  s(AvgL, bs = "ts", k = 8) +
                  s(AvgW, bs= "ts", k = 8) +
                  s(SpeciesNo, bs ="ad", k = 8) +
                  as.factor(Continent)+
                  te(Years,bs="ts", d =1),
                argGam = list(select=F),
                data = imp.dat
)
AIC(Q_qgam02)
e <- getViz(Q_qgam01)
a <- getViz(Q_qgam02)
e1<-e #gamViz object from above, edit variable name in plots as appropriate
a1 <- a
summary_model1 <- summary(a1)
summary_model2 <- summary(e1)
AIC(a1, e1)
print(plot(e1, allTerms = TRUE), pages = 1)

p_table1 <- data.frame(summary_model1$p.table) 
p_table2 <- data.frame(summary_model2$p.table) 

################################################################################
# Quantile Residual Diagnostics#################################################
################################################################################
# Identifies patterns in residuals across quantiles:
appraise(Q_qgam02, method = "simulate") 
appraise(Q_qgam02, method = "simulate") 

# Generate the diagnostic plots for 0.99 quantile
diag_plot01<- appraise(Q_qgam01, method = "simulate",
                       point_col = "#1E88E5", 
                       point_alpha = 0.6,
                       line_col = "#D81B60") &
  labs(title = NULL) &  # Removes all plot titles
  theme_minimal(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    text = element_text(family = "Arial"),
    plot.title = element_blank()  # Explicitly removes title space
  )


# Add panel labels (a-d) using patchwork
labeled_plot_1 <- diag_plot01 + 
  plot_layout(ncol = 2) &  # Arranges 4 plots in 2x2 grid
  plot_annotation(tag_levels = 'a'  # Uses lowercase letters
  ) & # Optional styling
  theme(plot.tag = element_text(size = 10, face = "bold"),
        plot.tag.position = c(0, 1))  # Top-left corner


# Generate the diagnostic plots for 0.5 quantile
diag_plo02 <- appraise(Q_qgam02, method = "simulate",
                       point_col = "#1E88E5",
                       point_alpha = 0.6,
                       line_col = "#D81B60") &
  labs(title = NULL) &
  theme_minimal(base_size = 9) &
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    text = element_text(family = "Arial"),
    plot.title = element_blank()
  )

# Add panel labels (a-d) using patchwork
labeled_plot <- diag_plo02 + 
  plot_layout(ncol = 2) &  # Arranges 4 plots in 2x2 grid
  plot_annotation(tag_levels = 'a'  # Uses lowercase letters
  ) & # Optional styling
  theme(plot.tag = element_text(size = 10, face = "bold"),
        plot.tag.position = c(0, 1))  # Top-left corner

################################################################################
## Visualization ###############################################################
################################################################################
# Ordinary least square plot####################################################
p_a <- ggplot(data = imp.dat, mapping = aes(x = AvgL,
                                            y = MPperind)) +
  geom_point(size= 1.5, color= "grey80") +
  xlab("Shrimp length (cm)") +
  ylab(bquote("Micoplastic abundance (MP/ind)")) +
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,
              color="grey20",lineend = "butt")+ggtitle("a")+
  geom_vline(xintercept = change_point_1, linetype ="dashed", color = "black", 
             alpha = 0.3)+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "MP = 0.6 + 0.8*SL",
           hjust = 0, vjust = 1, color= "black", size = 3)+ggtitle("a")+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = 0.3",
           hjust = 0, vjust = 2.5, color= "black", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))

p_b <- ggplot(data = imp.dat, aes(x = AvgL, y = MPperind)) +
  geom_point(size= 1.5, color= "grey80") +  # draw raw data points
  xlab("Shrimp length (cm)") +
  ylab(bquote("Microplastic concentration (MP/ind)"))+
  geom_ribbon(data = pred_l_data99, aes(ymin = lower, ymax = upper),
              fill = "#FA5447", alpha = 0.1)+
  geom_line(data = pred_l_data99, aes(y = MPperind), color = "#FA5447") +
  geom_vline(xintercept = change_point_2, linetype ="dashed", color = "black",
             alpha = 0.3)+  
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "MP = 11.1 + 1.7*SL",
           hjust = 0, vjust = -8, color= "#FA5447", size = 3)+ggtitle("d")+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = 0.07",
           hjust = 0, vjust = -6.6, color= "black", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))


# plot non-linear least square##################################################
p_nlr <- ggplot(data = imp.dat, aes(x = AvgL, y = MPperind)) +
  geom_point(size= 1.5, color= "grey80") +  # draw raw data points
  xlab("Shrimp length (cm)") +
  ylab(bquote("Microplastic abundance (MP/ind)"))+
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,
              color="grey20",lineend = "butt", formula = y ~ x +I(x^2))+
  geom_vline(xintercept = change_point_5, linetype ="dashed", color = "black", 
             alpha = 0.3)+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "MP = 1.98 + 0.36*SL - 0.02*SL^2",
           hjust = 0, vjust = 1, color= "black", size = 3)+ggtitle("b")+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = 0.30",
           hjust = 0, vjust = 2.5, color= "black", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))+
  theme(axis.title.y = element_blank()
  )


# plot quantile#################################################################
p_c <- ggplot(data = imp.dat, aes(x = AvgL, y = MPperind)) +
  geom_point(size= 1.5, color= "grey80") +  # draw raw data points
  xlab("Shrimp length (cm)") +
  ylab(bquote("Microplastic abundance (MP/ind)"))+
  geom_ribbon(data = pred_nl_data99, aes(ymin = lower, ymax = upper), 
              fill = "#FA5447", alpha = 0.1)+
  geom_vline(xintercept = change_point_3, linetype ="dashed", color = "black",
             alpha = 0.3)+
  geom_line(data = pred_nl_data99, aes(y = MPperind), color = "#FA5447") +  
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "MP = 1.7 + 4.4*SL - 0.1*SL^2",
           hjust = 0, vjust = -2, color= "#FA5447", size = 3)+ggtitle("e")+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = 0.33",
           hjust = 0, vjust = -0.5, color= "#FA5447", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))+
  theme(axis.title.y = element_blank()
  )



# gam model#####################################################################
library(mgcv)
library(tidymv)
p_gam <- plot_smooths(
  model = gam_model,
  series = AvgL
)+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = 0.35",
           hjust = 0, vjust = 4, color= "black", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))+
  ggtitle("c")+
  xlab("Shrimp length (cm)")+
  theme(axis.title.y = element_blank())


# this sm function can allows the grobs to the glist to combine the graph
p_d <- plot(sm(fitlv, 1),trans=function(x) x+coef(fitlv)[1])+
  l_points(shape = 18, size= 1.9, color= "grey80")+l_ciPoly(fill= "#FA5447", 
                                                            alpha =0.1)+
  l_fitLine(color ="#FA5447")+ggtitle("f")+
  xlab("Shrimp length (cm)")+
  theme(axis.title.y = element_blank())+
  annotate("Text", x = min(imp.dat$AvgL),
           y = max(imp.dat$MPperind),
           label = "R^2 = -0.97",
           hjust = 0, vjust = -8, color= "black", size = 3)+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))
summary(fitlv)
AIC(fitlv)
AIC(gam_model)

figure_5 <- mgcViz::gridPrint(p_a,p_nlr, p_gam,
                              p_b,p_c,p_d, nrow=2)
ggsave("figure_comp_mod_spp.pdf",plot= figure_5, width = 310, device = "pdf",
       height = 160, units = "mm", dpi = 300)


# QGAM plot#####################################################################

plot_a<-plot(sm(e1,1),trans=function(x) x+coef(e1)[1]
)+l_fitRaster()+l_rug()+l_fitContour()+ggtitle("e")+
  guides(fill=guide_colorbar("MP/ind"))+
  scale_fill_distiller(palette = "RdBu",type = "div") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))

plot_b<-plot(sm(e1,2),trans=function(x) x+coef(e1)[1]
)+l_rug()+l_ciPoly(fill= "#FA5447", alpha =0.1)+
  l_fitLine(color ="#FA5447")+ggtitle("f")+
  xlab("Shrimp length (cm)")+
  ylab("Microplastic abundance (MP/ind)")+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))

plot_c<-plot(sm(e1,3),trans=function(x) x+coef(e1)[1]
)+l_rug()+l_ciPoly(fill= "#00AFBB", alpha =0.1)+
  l_fitLine(color ="#00AFBB")+ggtitle("g")+
  xlab("Shrimp weight (gm)")+
  ylab("Microplastic abundance (MP/ind)") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))

plot_d<-plot(sm(e1,4),trans=function(x) x+coef(e1)[1]
)+l_rug()+l_ciPoly(fill= "green", alpha =0.1)+
  l_fitLine(color ="green")+ggtitle("h")+
  xlab("Species number (no.)")+
  ylab("Microplastic abundance (MP/ind)") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black",
                                        fill = NA))

plot_e<-plot(sm(a1,1),trans=function(x) x+coef(a1)[1]
)+l_fitRaster()+l_rug()+l_fitContour()+ggtitle("a")+
  guides(fill=guide_colorbar("MP/ind"))+
  scale_fill_distiller(palette = "RdBu",type = "div") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black", 
                                        fill = NA))

plot_f<-plot(sm(a1,2),trans=function(x) x+coef(a1)[1]
)+l_rug()+l_ciPoly(fill= "#FA5447", alpha =0.1)+
  l_fitLine(color ="#FA5447")+ggtitle("b")+
  xlab("Shrimp length (cm)")+
  ylab("Microplastic abundance (MP/ind)")+
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black", 
                                        fill = NA))

plot_g<-plot(sm(a1,3),trans=function(x) x+coef(a1)[1]
)+l_rug()+l_ciPoly(fill= "#00AFBB", alpha =0.1)+
  l_fitLine(color ="#00AFBB")+ggtitle("c")+
  xlab("Shrimp weight (gm)")+
  ylab("Microplastic abundance (MP/ind)") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black", 
                                        fill = NA))

plot_h<-plot(sm(a1,4),trans=function(x) x+coef(a1)[1]
)+l_rug()+l_ciPoly(fill= "green", alpha =0.1)+
  l_fitLine(color ="green")+ggtitle("d")+
  xlab("Species number (no.)")+
  ylab("Microplastic abundance (MP/ind)") +
  theme(panel.background = element_rect(linewidth = 0.8, colour = "black", 
                                        fill = NA))

################################################################################
# Figure 4 plot
figure_4 <- mgcViz::gridPrint(plot_e,plot_f,plot_g, plot_h,
                              plot_a,plot_b,plot_c, plot_d,nrow=2)

################################################################################
# Changepoint and threshold detection###########################################
################################################################################

# Get fitted values and predictor data (0.99 QGAM)
fit_data <- data.frame(
  Length = Q_qgam01$model$AvgL,  # Replace AvgL with your actual length variable
  MP = fitted(Q_qgam01),
  Changepoint = FALSE
)


# 1. Run changepoint analysis
cpt <- changepoint::cpt.meanvar(fit_data$MP, method = "PELT", penalty = "BIC")
threshold_positions <- fit_data$Length[changepoint::cpts(cpt)]

# 2. Create SEGMENTATION dataframe (renamed from 'segments')
seg_df <- data.frame(
  start = c(min(fit_data$Length), threshold_positions),
  end = c(threshold_positions, max(fit_data$Length))
)

# 3. Calculate segment statistics
seg_df$mean_MP <- sapply(1:nrow(seg_df), function(i) {
  mean(fit_data$MP[fit_data$Length >= seg_df$start[i] & 
                     fit_data$Length <= seg_df$end[i]])
})

seg_df$n <- sapply(1:nrow(seg_df), function(i) {
  sum(fit_data$Length >= seg_df$start[i] & 
        fit_data$Length <= seg_df$end[i])
})

# 4. Filter meaningful thresholds (10% change & min 10 obs)
meaningful_thresholds <- seg_df %>% 
  mutate(MP_change = abs(mean_MP - lag(mean_MP))/mean_MP) %>%
  filter(MP_change > 0.10, n > 15) %>% 
  arrange(desc(MP_change))


cpt_plot <- ggplot(fit_data, aes(Length, MP)) +
  # Data points (Nature-style muted colors)
  geom_point(alpha = 0.8, size = 2, color = "#DEDEDE") +
  
  # Threshold lines (color-blind friendly palette)
  geom_vline(
    data = meaningful_thresholds,
    aes(xintercept = start),
    color = c("#009E73", "#56B4E9", "#F4A637"), # Nature Communications palette
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  #
  
  # Nature-style elements
  labs(
    x = "Shrimp length (cm)",
    y = "Microplastic particles per individual"
  ) +
  
  theme_minimal(base_size = 8) + # Nature typically uses 8pt base
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    plot.caption = element_text(
      hjust = 0, 
      size = 7 # Fixed margin error
    ),
    plot.margin = unit(c(4, 8, 4, 4), "mm") # Top, right, bottom, left
  ) +
  
  # Statistical annotation
  annotate(
    "text", 
    x = min(fit_data$Length), 
    y = max(fit_data$MP) * 0.95,
    label = paste(
      "PELT changepoint detection\n",
      "BIC penalty |",
      "Min. ΔMP > 10%"
    ),
    hjust = -1.8, vjust = 14, size = 2.5, color = "grey30"
  )