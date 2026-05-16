# load the newest fields library
library(fields)
# make sure to set working directory to the location of this file by
# Session --> Set Working Directory --> To Source File Location
###########################

source("Func_InterpolationError.R")

range <- 0.25
#smoothness_grid <- c(0.5, 1, 1.5, 2.5, 3.5)
smoothness_grid <- c(0.5, 1, 1.5, 2.5)
sObs <- rbind( c( .5,.5))
grid_lengthout_tick <- c(50, 60, 80, 100, 200, 350, 500)
grid_lengthout <- c(40, 50, 60, 70, 80, 90, 120, 150, 200, 250, 300, 350, 400, 500)

#NN_grid <- seq(from = 2, to = 100, by = 2)
maxErrors <- matrix(NA, nrow = length(grid_lengthout), 
                    ncol = length(smoothness_grid))
fill_distance <- c()
for (j in 1:length(smoothness_grid)) {
  for (i in 1:length(grid_lengthout)){
    input_gridList <- list(x = seq(from = 0, to = 1, length.out = grid_lengthout[i]), 
                           y = seq(from = 0, to = 1, length.out = grid_lengthout[i]))
    temp_result <- interpolation_error(sObs,
                                       input_gridList,
                                       NNSize = 2, 
                                       aRange = range,
                                       smoothness = smoothness_grid[j])
    maxErrors[i, j] <- temp_result$maxErr
    NNGrid <- temp_result$NNGrid
    fineGrid <- temp_result$finePredGrid
    fill_distance[i] <- min(rdist(sObs, temp_result$NNGrid))
  }
}

plot(log10(fill_distance), log10(maxErrors[,1]), ylim = c(-10, 0),
     xlim = c(-10, 0))
abline(lm(log10(maxErrors[,1]) ~ log10(fill_distance)))
points(log10(fill_distance), log10(maxErrors[,2]), col = "red")
abline(lm(log10(maxErrors[,2]) ~ log10(fill_distance)), col = "red")
points(log10(fill_distance), log10(maxErrors[,3]), col = "blue")
abline(lm(log10(maxErrors[,3]) ~ log10(fill_distance)), col = "blue")
points(log10(fill_distance), log10(maxErrors[,4]), col = "green")
abline(lm(log10(maxErrors[,4]) ~ log10(fill_distance)), col = "green")



#### manually do log w. fill distance
pdf("../Plots/Fig05_MaxErrorFillDistance.pdf", width = 8, height = 4.5)
lwd_guides = 3
par(mar = c(4, 4.5, 3, 8))
x_vals <- 1 / fill_distance

#ylim_axis <- c(min(c(maxErrors)), max(c(maxErrors)))
ylim_axis <- c(min(c(maxErrors)), 0.1)
plot(x_vals, maxErrors_L2[, 1], pch = 15, 
     xlab = expression("Grid Size"), 
     ylab = "Maximum Absolute Error", 
     main = "",
     ylim = ylim_axis, col = lightGray, 
     log = 'xy', xaxt = "n")
points(x_vals, maxErrors[, 2], pch = 16,
       col = goldenTech)
points(x_vals, maxErrors[, 3], pch = 17,
       col = earthBlue)
points(x_vals, maxErrors[, 4], pch = 18,
       col = coloradoRed, cex = 1.3)
lines(x_vals, rep(mean(maxErrors[,1]), length(fill_distance)), 
      col = lightGray, lty = 2, lwd = lwd_guides)
lines(x_vals, 10^(log10(fill_distance^(1/2)) - 2.5), col = goldenTech,
      lty = 3, lwd = lwd_guides)
lines(x_vals, 10^(log10(fill_distance) - 2.8), col = earthBlue,
      lty = 4, lwd = lwd_guides)
lines(x_vals, 10^(log10(fill_distance^2) - 2), col = coloradoRed,
      lty = 5, lwd = lwd_guides)
## Custom x-axis: "50x50", "60x60", ...
axis(1, at = x_vals,
     labels = paste0(grid_lengthout, "x", grid_lengthout))
axis(3,
     at = x_vals,
     labels = format(signif(x_vals, 2), scientific = TRUE))
mtext(expression(1/delta), side = 3, line = 2)  # optional top-axis title
legend("topleft", title = "Observed Error",
       col = c(lightGray, goldenTech, earthBlue, coloradoRed),
       pch = c(15, 16, 17, 18), pt.cex = c(1, 1, 1, 1.3),
       legend = c(expression(nu == 1/2), expression(nu == 1), 
                  expression(nu == 3/2), expression(nu == 5/2)),
       inset=c(1.015, 0.05), bty = "n",
       xpd=TRUE, cex = 1)
legend("topleft", title = "Theoretical Order",
       col = c(lightGray, goldenTech, earthBlue, coloradoRed),
       lty = c(2, 3, 4, 5), lwd = lwd_guides,
       legend = c(expression(delta^{0}), expression(delta^{1/2}), 
                  expression(delta^{1}), expression(delta^{2})),
       inset=c(1.015, 0.5), bty = "n",
       xpd=TRUE)

dev.off()
