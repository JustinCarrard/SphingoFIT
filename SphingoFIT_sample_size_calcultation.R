#==========================================================================================
# Sample size calculation for SpingoFIT
# Author: Denis Infanger
#==========================================================================================
#------------------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

#------------------------------------------------------------------------------------------
# Helper functions
#------------------------------------------------------------------------------------------

# Function to calculate the noncentrality parameter of the t-distribution for the sample size calculation

delta <- function(n2, k = 1, d) {
  d*sqrt((k*n2)/(1 + k))
}

# Function to calculate sample size for independent t-test on the follow-up values with unequal allocation ratio k

calc_power <- function(k, d, alpha, power) { # Page 48 in Shein-Chung Chow et al. (2018), 3rd ed.
  # k = allocation ratio (n_intervention/n_controls)
  # d = Cohen's d (effect size)
  # alpha = significance level (e.g. 0.05)
  # power = statistical power (e.g. 0.8)
  
  optim_fun <- function(nC, k, d, alpha, power) {
    (1 - pt(qt(1 - alpha/2, df = nC*(k + 1) - 2), df = nC*(k + 1) - 2, ncp = delta(nC, k = k,  d = d)) +
       pt(-qt(1 - alpha/2, df = nC*(k + 1) - 2), df = nC*(k + 1) - 2, ncp = delta(nC, k = k,  d = d))) - power
  }
  
  uniroot(optim_fun, lower = 1.01, upper = 1e6, d = d, k = k, alpha = alpha, power = power)$root
}

#------------------------------------------------------------------------------------------
# Sample size calculation
#------------------------------------------------------------------------------------------
# Set drop out rate

dropout <- 0.1

# The assumed standard deviation: Average SD of sphingolipids on log2-scale of COmPLETE Health study participants

assumed_sd <- 0.475

# Convert published differences on original scale to log2-scale (Kasumov et al. 2015)
# This assumes that the published raw data follows a log-normal distribution and that we transform it using log2 (not log!)

lnmean_1 <- 5.692 # Mean of lognormal distribution
lnsd_1 <- 0.199*sqrt(13) # SD of lognormal distribution (SEM is reported)
lnmean_2 <- 4.892
lnsd_2 <- 0.246*sqrt(13)

mean_1 <- log(lnmean_1^2/sqrt(lnsd_1^2 + lnmean_1^2))/log(2) # Mean on log2-scale (normal distribution when lognormal is log2(x) transformed)
# sd_1 <- (sqrt(log(1 + lnsd_1^2/lnmean_1^2)))/log(2) # SD on log2-scale (normal distribution)

mean_2 <- log(lnmean_2^2/sqrt(lnsd_2^2 + lnmean_2^2))/log(2)
# sd_2 <- (sqrt(log(1 + lnsd_2^2/lnmean_2^2)))/log(2)

geomean_g1 <- 2^(mean_1) # Geometric mean 1
geomean_g2 <- 2^(mean_2) # Geometric mean 2

geomean_g1/geomean_g2 # Geometric mean ratio (does not depend on base of logarithm used)

diff_normalscale <- abs(mean_2 - mean_1) # Differences of means on log2-scale (normal distributions)

rho <- seq(0.05, 0.95, 0.05)
pow <- c(0.8, 0.9)

# Calculate sample size according to Borm et al. (2007): 10.1016/j.jclinepi.2007.02.006

required_n_ancova <- sapply(rho, FUN = function(x) {
  sapply(pow, function(y) {
    ceiling(power.t.test(
      delta = diff_normalscale
      , sd = sqrt(assumed_sd^2*(1 - x^2))
      , sig.level = 0.05
      , power = y
      , type = "two.sample"
      , alternative = "two.sided"
      , strict = TRUE
    )$n*(1/(1 - dropout)))
  })
})

# Set up data frame to store results for plotting

plot_dat <- data.frame(
  rho = rep(rho, times = length(pow))
  # , detectable_diff = 2^(as.vector(t(detectable_diff_ancova)))
  , reqN = as.vector(t(required_n_ancova))
  , power = factor(paste0(100*rep(pow, each = length(rho)), "%"))
)

# Create plot

theme_set(theme_bw())
p <- ggplot(data = subset(plot_dat, power %in% "80%"), aes(x = rho, y = reqN, group = power)) +
  # geom_hline(yintercept = exp(detectable_diff_ttest), col = "red", linetype = 2) +
  geom_point(aes(colour = power), size = 4) +
  geom_line(aes(colour = power)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    y = "Required N per group"
    , x = expression("Correlation coefficient "*rho*" between pre and post measurements")
    , title = paste0("Assumed standard deviation: ", signif(assumed_sd, 3), "; assumed geometric mean ratio: ", round(2^(diff_normalscale), 2), "; assumed dropout: ", dropout*100, "%")
  ) +
  scale_colour_manual(name = "Power", values = c("#008FD0", "#F07E00")) +
  theme(
    axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, margin=margin(0,12,0,0)),
    axis.title.x=element_text(colour = "black", size = 17, margin=margin(12,0,0,0)),
    # axis.title.y=element_text(size=15,hjust=0.5, vjust=1),
    axis.text.x=element_text(colour = "black", size=15),
    axis.text.y=element_text(colour = "black", size=15),
    # plot.margin=unit(c(2,2,2,2,2),"line"),
    legend.position="right",
    legend.text = element_text(size=15),
    legend.title = element_text(size = 15),
    panel.grid.minor = element_blank(),
    # panel.grid.major = element_blank(),
    legend.key.size = unit(1, "cm"),
    plot.title = element_text(face = "bold"),
    strip.background=element_rect(fill="white"),
    strip.text.x=element_text(size=12)
    # , text = element_text(family = "Lato")
  )
# +
#   annotate("text", x = 0.1, y = 2.05, label = "80% Power", size = 5) +
#   annotate("text", x = 0.1, y = 2.3, label = "90% Power", size = 5)

p