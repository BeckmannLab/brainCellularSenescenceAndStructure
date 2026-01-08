library(data.table)
library(qvalue)
# Simulation parameters

setwd("/bulk_smri/")

resFinal=fread("allSTRDEGresults_01.18.24.txt",data.table=FALSE)

subset_res = resFinal[resFinal$feature == "rh_superiorfrontal_area",]
pvals <- subset_res$P.Value

# Clean and ensure p-values are between 0 and 1 (remove NAs, etc.)
pvals <- pvals[!is.na(pvals)]
pvals <- pvals[pvals >= 0 & pvals <= 1]

# Estimate pi0 and thus pi1
qout    <- qvalue(pvals)
pi0_est <- qout$pi0
pi1_est <- 1 - pi0_est

# Compute histogram density values without plotting
hist_data <- hist(
  pvals,
  breaks = 50,
  plot = FALSE,
  freq = FALSE
)

# Determine max y value among histogram and reference lines
y_max <- max(max(hist_data$density), 1, pi0_est)

# Plot the histogram with proper y-axis limits
hist(
  pvals,
  breaks = 50,
  freq = FALSE,
  xlim = c(0, 1),
  ylim = c(0, y_max * 1.05),  # add a bit of headroom
  xlab = "P-Value",
  main = bquote("P-Value Distribution (aprox." ~ pi[1] ~ " = " ~ .(round(pi1_est, 2)) ~ ")")
)

# Add reference line for uniform null distribution (y=1)
abline(h = 1, col = "red", lwd = 2, lty = 2)

# Add dashed green line for estimated pi0
abline(h = pi0_est, col = "forestgreen", lwd = 3, lty = 2)

# Add legend with a hat over pi0
legend(
  "topright",
  legend = c(
    expression("Null Uniform Density (y=1)"),
    bquote(y == hat(pi)[0] ~ "=" ~ .(format(pi0_est, digits=2)))
  ),
  col = c("red", "forestgreen"),
  lty = c(2, 2),
  lwd = c(2, 3),
  bty = "n"
)

