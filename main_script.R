
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

wd <- file.path('~', 'microecon_final')
setwd(wd)

df <- read_dta("Bonjour_data.dta")

glimpse(df)

aggregate(df[, c('bweight', 'earning', 'schyear', 'highqua', 'age', 'married', 'full', 'own_exp', 'sm16', 'sm18')],
                   list(df$twinno), function(x) c(round(mean(x, na.rm = TRUE), 2)))

#working group
df[, c('bweight', 'earning', 'schyear', 'highqua', 'age', 'married', 'full', 'own_exp', 'sm16', 'sm18')] %>%
  filter((df$full == 1) | (df$part == 1) | (df$self == 1)) %>%
  summarize_all(list(mean), na.rm = TRUE)

#pooled
df[, c('bweight', 'earning', 'schyear', 'highqua', 'age', 'married', 'full', 'own_exp', 'sm16', 'sm18')] %>%
  summarize_all(list(mean), na.rm = TRUE)


dfsum <- df %>% dfSummary(., plain.ascii = FALSE, style = "grid", 
                                                    graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = "/tmp") 

dfsum$Missing <- NULL
view(dfsum, file = "~/microecon_final/dfsum.html")

cormatdf <- df[, c('bweight',
                    'earning',
                    'schyear',
                    'highqua',
                    'age',
                    'own_exp')]

cormat <- cormatdf %>% cor(use = "complete.obs") %>% round(., 2)

# Reorder matrix according to coefficients and hierarchical clustering order
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  dd[is.na(dd)] <- 0
  dd[is.nan(dd)] <- 0
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

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

ggheatmap + 
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
  labs(title = "Correlation heatmap between the non-categorial variables",
       subtitle = "ordered according to hclust")
