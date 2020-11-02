
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
library(Hmisc)
library(corrplot)

wd <- file.path('~', 'microecon_final')
setwd(wd)

df <- read_dta("Bonjour_data.dta")

# Descriptive stats -------------------------------------------------------

dfsum <- df[, c('bweight', 'earning', 'schyear', 'highqua', 'age', 'married', 'full', 'own_exp', 'sm16', 'sm18')] %>%
  dfSummary(., plain.ascii = FALSE, style = "grid", 
                                                    graph.magnif = 0.75, valid.col = FALSE, tmp.img.dir = "/tmp") 

dfsum$Missing <- NULL
view(dfsum, file = "~/microecon_final/dfsum.html")

# Feature Selection -----------------------------------------------------------------

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

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#basic features

cormatdf <- df[, c('highqua',
                   'bweight',
                   'full',
                   'part',
                   'self',
                   'married',
                   'own_exp',
                   'exp_par')]

res1 <- rcorr(as.matrix(cormatdf))
res11 <- flattenCorrMatrix(res1$r, res1$P)

write.table(subset(res11, row == 'highqua'), file = "~/microecon_final/res1.txt", sep="\t")

# Insignificant correlation are crossed
corrplot(res1$r, type="upper", order="hclust", p.mat = res1$P, sig.level = 0.05, tl.col = "black", tl.srt = 45)

#differenced features

cormatdf <- pdf[, c('dhighqua',
                    'dbweight',
                    'dfull',
                    'dpart',
                    'dself',
                    'dmarried',
                    'down_exp',
                    'dexp_par')]

res2 <- rcorr(as.matrix(cormatdf))
res22 <- flattenCorrMatrix(res2$r, res2$P)

write.table(subset(res22, row == 'dhighqua'), file = "~/microecon_final/res2.txt", sep="\t")

# Insignificant correlation are crossed
corrplot(res2$r, type="upper", order="hclust", p.mat = res2$P, sig.level = 0.05, tl.col = "black", tl.srt = 45)

#smoking as IV?

cormatdf <- pdf[, c('highqua',
                    'sm16',
                    'sm18')]

res3 <- rcorr(as.matrix(cormatdf))
res33 <- flattenCorrMatrix(res3$r, res3$P)

write.table(subset(res33, row == 'highqua'), file = "~/microecon_final/res3.txt", sep="\t")

# Insignificant correlation are crossed
p3 <- corrplot(res3$r, type="upper", order="hclust", 
               p.mat = res3$P, sig.level = 0.05, tl.col = "black", tl.srt = 45)

#differencing smoking
cormatdf <- pdf[, c('dhighqua',
                    'dsm16',
                    'dsm18')]

res4 <- rcorr(as.matrix(cormatdf))
res44 <- flattenCorrMatrix(res4$r, res4$P)

write.table(subset(res44, row == 'dhighqua'), file = "~/microecon_final/res4.txt", sep="\t")

# Insignificant correlation are crossed
p4 <- corrplot(res4$r, type="upper", order="hclust", 
               p.mat = res4$P, sig.level = 0.05, tl.col = "black", tl.srt = 45)

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
          align = TRUE, type = 'text',  out = "models1.htm")

#smoking as IV
ivsm16 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) | age + I(age**2/100) + sm16,  data=df)
ivsm18 <- ivreg(log(earning) ~ highqua + age + I(age**2/100) | age + I(age**2/100) + sm18,  data=df)

stargazer(ols1, ivsm16, ivsm18, type = 'text', out = "models2.htm")

