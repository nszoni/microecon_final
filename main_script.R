
# NFO ---------------------------------------------------------------------

# Final Project.
# Authors: Winkler Soma, Alexandrov Dániel, Nguyen Nam Son
# Date: 19-11-2020

# Setup -------------------------------------------------------------------

library(dplyr)
library(AER)
library(haven)
library(plm)
library(car)
library(summarytools)
library(stargazer)
library(ggplot2)
library(reshape2)
library(patchwork)

wd <- file.path('~', 'microecon_final')
setwd(wd)

df <- read_dta("Bonjour_data.dta")

glimpse(df)

#pooled

dfsum <- df[, c('bweight', 'earning', 'schyear', 'highqua', 'age', 'married', 'full', 'own_exp', 'sm16', 'sm18')] %>%
  dfSummary(., plain.ascii = FALSE, style = "grid", 
                                                    graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = "/tmp") 

dfsum$Missing <- NULL
view(dfsum, file = "~/microecon_final/dfsum.html")

# Methods -----------------------------------------------------------------

# simple OLS

ols1 <- lm(log(earning) ~ highqua + age + I(age**2/100), data=df)

iv1 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) | age + I(age**2/100) + twihigh,  data=df)

ols2 <- lm(log(earning) ~ highqua + age + I(age**2/100) + LNandSE + married + own_exp + part + bweight, data=df)

iv2 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) + LNandSE + married + own_exp + part + bweight| twihigh + age + I(age**2/100) + LNandSE + married + own_exp + part,  data=df)

#panel
pdf <- pdata.frame(df, index = c("family", "twinno"))
pdim(pdf)

fdr_ols1 <- plm(log(earning) ~ highqua + age + I(age**2/100), data=df, model = 'fd')

fdr_iv1 <- plm(log(earning) ~ highqua + age + I(age**2/100) | twihigh + age + I(age**2/100), data=df, model = 'fd')

fdr_ols2 <- plm(log(earning) ~ highqua + age + I(age**2/100) + LNandSE + married + own_exp + part + bweight, data=df, model = 'fd')

fdr_iv2 <- plm(log(earning) ~ highqua + age + I(age**2/100) + LNandSE + married + own_exp + part + bweight | twihigh + age + I(age**2/100) + LNandSE + married + own_exp + part + bweight, data=df, model = 'fd')

#comparing all the models together

stargazer(ols1, iv1, ols2, iv2, fdr_ols1, fdr_iv1, fdr_ols2, fdr_iv2, font.size = "small",
          align = TRUE, type = 'text',  out="models1.htm")

#smoking as IV
ivsm16 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) | age + I(age**2/100) + sm16,  data=df)
ivsm18 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) | age + I(age**2/100) + sm18,  data=df)

stargazer(ols1, ivsm16, ivsm18, type = 'text', out = "models2.htm")

# Correlations -----------------------------------------------------------------

pdf$dhighqua <- diff(pdf$highqua)
pdf$dmarried <- diff(pdf$married)
pdf$dbweight <- diff(pdf$bweight)
pdf$dexp_par <- diff(pdf$exp_par)
pdf$down_exp <- diff(pdf$own_exp)
pdf$dfull <- diff(pdf$full)
pdf$dpart <- diff(pdf$part)
pdf$dself <- diff(pdf$self)
pdf$dsm16 <- diff(pdf$sm16)
pdf$dsm18 <- diff(pdf$sm18)

cormatdf <- df[, c('highqua',
                   'bweight',
                   'full',
                   'part',
                   'self',
                   'married',
                   'own_exp',
                   'exp_par')]

cormat <- cormatdf %>% cor(use = "complete.obs") %>% round(., 2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

cormat <- reorder_cormat(cormat)
lower_tri <- get_lower_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1) )+
  coord_fixed()

p1 <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = "Correlation heatmap between non-differenced variables",
       subtitle = "for feature selection of control variables")

cormatdf <- pdf[, c('dhighqua',
                    'dbweight',
                    'dfull',
                    'dpart',
                    'dself',
                    'dmarried',
                    'down_exp',
                    'dexp_par')]

cormat <- cormatdf %>% cor(use = "complete.obs") %>% round(., 2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

cormat <- reorder_cormat(cormat)
lower_tri <- get_lower_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1) )+
  coord_fixed()

p2 <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = "Correlation heatmap of differenced variables",
       subtitle = "differencing reduced the correlation between regressors")

p1 + p2

#smoking as IV?

cormatdf <- pdf[, c('highqua',
                    'sm16',
                    'sm18')]

cormat <- cormatdf %>% cor(use = "complete.obs") %>% round(., 2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

cormat <- reorder_cormat(cormat)
lower_tri <- get_lower_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1) )+
  coord_fixed()

p3 <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = "Correlation heatmap of differenced variables",
       subtitle = "differencing reduced the correlation between education and smoking")

#differencing smoking
cormatdf <- pdf[, c('dhighqua',
                    'dsm16',
                    'dsm18')]

cormat <- cormatdf %>% cor(use = "complete.obs") %>% round(., 2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

cormat <- reorder_cormat(cormat)
lower_tri <- get_lower_tri(cormat)

# Melt the correlation matrix
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1) )+
  coord_fixed()

p4 <- ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(title = "Correlation heatmap of differenced variables",
       subtitle = "differencing reduced the correlation between education and smoking")

p3 + p4
