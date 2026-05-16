# load the newest fields library
library(fields)
# make sure to set working directory to the location of this file by
# Session --> Set Working Directory --> To Source File Location

###### School colors ######
# Webcolors
paleBlue = "#CFDCE9"; lightBlue = "#879EC3"; blasterBlue = "#09396C"
darkBlue = "#21314d"; coloradoRed = "#CC4628"
# Neutrals
lightGray = "#AEB3B8"; silver = "#81848A"; darkGray = "#75757D"
# Accents
goldenTech = "#F1B91A"; earthBlue = "#0272DE"; mutedBlue = "#57A2BD"
energyYellow = "#F0F600"; envrGreen = "#80C342"; redFlannel = "#B42024"
###########################

nObs <- c(200, 1500, 6500)
mList <- c(10, 20, 40, 60, 100, 140, 200, 260, 350, 500)
sim_predTiming_exact <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
sim_predTiming_NN2 <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
sim_predTiming_NN4 <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
sim_predTiming_NN8 <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
  
RUN = FALSE

if (RUN == TRUE){
  set.seed(123)
  for (i in 1:length(nObs)) {
    nObs_temp <- nObs[i]
    s_temp <- cbind(x = runif(nObs_temp), y = runif(nObs_temp))
    dist.mat_temp <- rdist(s_temp, s_temp)
    
    SigmaCov_temp <- Matern(dist.mat_temp, range = 0.05, smoothness = 1) 
    Sigma.c_temp <- chol(SigmaCov_temp)
    y0_temp <- t(Sigma.c_temp) %*% rnorm(nObs_temp) + rnorm(nObs_temp, mean = 0, sd = sqrt(0.2))
    
    # fitting the process
    krigObj_temp <- spatialProcess(s_temp, c(y0_temp), lambda = 0.2, 
                                   simpleKriging = TRUE,
                                   cov.args = list(Covariance = "Matern", aRange = 0.05, 
                                                   smoothness = 1, sigma2 = 1) )
    for (j in 1:length(mList)) {
      pred_grid_side <- seq(from = 0, to = 1, length.out = mList[j])
      
      # for exact
      t_temp <- c()
      print(paste("Starting nObs =", nObs_temp, 
                  "; exact method size = ", mList[j]))
      if (mList[j] < 350 & nObs_temp < 5000){
        for (k in 1:20) {
          t_temp[k] <- system.time(
            pred <- predictSurface(krigObj_temp, 
                                   gridList = list(x = pred_grid_side, y = pred_grid_side),
                                   extrap = TRUE)
          )[3]
        }
      } else if (nObs_temp > 5000 & mList[j] > 150) {
        t_temp[1] <- NA
      } else {
        t_temp[1] <- system.time(
          pred <- predictSurface(krigObj_temp, 
                                 gridList = list(x = pred_grid_side, y = pred_grid_side),
                                 extrap = TRUE)
        )[3]
      }
      sim_predTiming_exact[j, i] <- median(t_temp)
       
      # NN = 2
      t_temp <- c()
      print(paste("Starting nObs =", nObs_temp, 
                  "; NN2 size = ", mList[j]))
      for (k in 1:20) {
        t_temp[k] <- system.time(
          pred <- predictSurface(krigObj_temp, 
                               gridList = list(x = pred_grid_side, y = pred_grid_side), 
                               fast = TRUE, NNSize = 2, extrap = TRUE)
        )[3]
      }
      sim_predTiming_NN2[j, i] <- median(t_temp)
      
      # NN = 4
      t_temp <- c()
      for (k in 1:20) {
        t_temp[k] <- system.time(
          pred <- predictSurface(krigObj_temp, 
                                 gridList = list(x = pred_grid_side, y = pred_grid_side), 
                                 fast = TRUE, NNSize = 4, extrap = TRUE)
        )[3]
      }
      sim_predTiming_NN4[j, i] <- median(t_temp)
      
      # NN = 8
      t_temp <- c()
      print(paste("Starting nObs =", nObs_temp, 
                  "; NN8 size = ", mList[j]))
      for (k in 1:20) {
        t_temp[k] <- system.time(
          pred <- predictSurface(krigObj_temp, 
                                 gridList = list(x = pred_grid_side, y = pred_grid_side), 
                                 fast = TRUE, NNSize = 8, extrap = TRUE)
        )[3]
      }
      sim_predTiming_NN8[j, i] <- median(t_temp)
    }
  }

  save(sim_predTiming_exact, sim_predTiming_NN2, sim_predTiming_NN4,
       sim_predTiming_NN8, file = "SimPredTiming_New.RData")
} else {
  load("SimPredTiming.RData")
}

pdf("../Plots/Fig06_SimPredictionTiming.pdf", width=8, height=3)

color_kriging = goldenTech
color_ourMethod = blasterBlue

ylim_axis <- c(0.002, 61)

# set layout and outer margin (oma) to leave space on the right
par(mfrow = c(1, 3), oma = c(0, 0, 0, 7.5), mar = c(4, 4, 2, 1), xpd = NA)
plot(x = mList[4:10]^2, y = sim_predTiming_exact[4:10, 1], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, cex.lab = 1.3,
     ylab = "Time [s]", xlab = "", col = color_kriging, 
     ylim = ylim_axis)
     #ylim = c(0.001, max(c(sim_predTiming_exact[,1]))))
side_labels <- paste0(mList[4:10], " x ", mList[4:10])
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:10]^2, sim_predTiming_NN2[4:10, 1], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN4[4:10, 1], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN8[4:10, 1], type = 'o', pch = 19, 
      lwd = 3, lty = 2, col = color_ourMethod)
mtext("(a)", line=.5, adj=.05)

par(mar = c(4, 3, 2, 2))
plot(x = mList[4:10]^2, y = sim_predTiming_exact[4:10, 2], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, cex.lab = 1.3,
     xlab = "Prediction Grid Size", ylab = "", col = color_kriging,
     ylim = ylim_axis)
     #ylim = c(0.005, max(sim_predTiming_exact[,2])))
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:10]^2, sim_predTiming_NN2[4:10, 2], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN4[4:10, 2], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN8[4:10, 2], type = 'o', pch = 19, 
      lwd = 3, lty = 2, col = color_ourMethod)
mtext("(b)", line=.5, adj=.05)

#par(mar=c(4,4,2,10))
plot(x = mList[4:10]^2, y = sim_predTiming_exact[4:10, 3], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, 
     xlab = "", ylab = "", col = color_kriging,
     ylim = ylim_axis)
     #ylim = c(0.03, 16))
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:10]^2, sim_predTiming_NN2[4:10, 3], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN4[4:10, 3], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = color_ourMethod)
lines(mList[4:10]^2, sim_predTiming_NN8[4:10, 3], type = 'o', pch = 19, 
      lwd = 3, lty = 2, col = color_ourMethod)
mtext("(c)", line=.5, adj=.05)

legend("topleft", lty = c(1, 2, 1, 3),           # line types
       lwd = 3,                       # line width
       pch = c(19, 19, 19, 19),       # point symbols
       pt.cex = 1,                    # point size
       col = c(color_kriging, color_ourMethod, color_ourMethod, color_ourMethod),
       legend = c("Exact Method", "Rapid, L = 8", 
                  "Rapid, L = 4", "Rapid, L = 2"),
       inset = c(1.05, 0.3), xpd = NA, bty = "n", cex = 1.1)

dev.off()