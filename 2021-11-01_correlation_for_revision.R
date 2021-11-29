# correlation for revision
library(data.table)

setwd("~/Desktop/climate_niche_seed_dispersal/seed_dispersal")
trait.dir <- paste0(getwd(),"/trait_data")
niche_data <- list.files(paste0(trait.dir), "niche.csv")

all_niche_tables <- lapply(paste0(trait.dir, "/",niche_data), fread)
all_clades <- do.call( rbind, all_niche_tables)
#---------------
# Correlation ai ~ precipitation
x<-lm(all_clades$mean_aridity~all_clades$mean_prec)
summary(x)
#all:
#  lm(formula = all_clades$mean_aridity ~ all_clades$mean_prec)
#
#Residuals:
#  Min       1Q   Median       3Q      Max 
#-1.92063 -0.18795 -0.03186  0.15905  1.65119 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           1.74372    0.04808   36.27   <2e-16 ***
#  all_clades$mean_prec  1.02691    0.00680  151.01   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.3158 on 3173 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared:  0.8779,	Adjusted R-squared:  0.8778 
#F-statistic: 2.28e+04 on 1 and 3173 DF,  p-value: < 2.2e-16

#---------------
# Correlation ai ~ temperature
y<-lm(all_clades$mean_aridity~all_clades$mean_temp)
summary(y)
#Call:
#  lm(formula = all_clades$mean_aridity ~ all_clades$mean_temp)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-6.9324 -0.3251  0.2193  0.5316  2.8273 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           -23.015      4.243  -5.425 6.23e-08 ***
#  all_clades$mean_temp    5.637      0.748   7.536 6.31e-14 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.8955 on 3173 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared:  0.01758,	Adjusted R-squared:  0.01727 
#F-statistic: 56.78 on 1 and 3173 DF,  p-value: 6.306e-14

#---------------
# Correlation precipitation ~ temperature
z<-lm(all_clades$mean_prec~all_clades$mean_temp)
summary(z)
#Call:
#  lm(formula = all_clades$mean_prec ~ all_clades$mean_temp)
#
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-6.0596 -0.2565  0.1203  0.4509  1.8240 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -68.1773     3.6681  -18.59   <2e-16 ***
#  all_clades$mean_temp  13.2583     0.6467   20.50   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.7745 on 3174 degrees of freedom
#Multiple R-squared:  0.1169,	Adjusted R-squared:  0.1167 
#F-statistic: 420.3 on 1 and 3174 DF,  p-value: < 2.2e-16
