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
sim_CSTiming_exact <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
sim_CSTiming_fast <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
sim_CSTiming_rapid <- data.frame(matrix(NA, nrow = length(mList), ncol = length(nObs)))
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
      predGrid <- make.surface.grid(list(x1 = pred_grid_side, x2 = pred_grid_side))
      # for exact
      t_temp <- c()
      print(paste("Starting nObs =", nObs_temp, 
                  "; exact method size = ", mList[j]))
      if (mList[j] < 130){
        for (k in 1:10) {
          t_temp[k] <- system.time(
            cs <- sim.spatialProcess(krigObj_temp, predGrid, M = 10)
          )[3]
        }
      } else {
        t_temp[1] <- NA
      } 
      sim_CSTiming_exact[j, i] <- median(t_temp)
       
      # Maggie's
      print(paste("Starting nObs =", nObs_temp, 
                  ";Maggie's size = ", mList[j]))
      t_temp <- c()
      if (mList[j] < 130){
        for (k in 1:20) {
        t_temp[k] <- system.time(
          cs <- simLocal.spatialProcess(krigObj_temp, 
                                        list(x = pred_grid_side, y = pred_grid_side),
                                        extrap = TRUE,
                                        M = 10, NNSize = 4)
        )[3]
        }
      } else if (mList[j] < 400 & nObs_temp < 3000) {
        t_temp[1] <- system.time(
          cs <- simLocal.spatialProcess(krigObj_temp, 
                                        list(x = pred_grid_side, y = pred_grid_side),
                                        extrap = TRUE,
                                        M = 10, NNSize = 4)
        )[3]
      } else {
        t_temp[1] <- NA
      }
      sim_CSTiming_fast[j, i] <- median(t_temp)
      
      # rapid
      print(paste("Starting nObs =", nObs_temp, 
                  ";Maggie's + Rapid size = ", mList[j]))
      t_temp <- c()
      for (k in 1:20) {
        t_temp[k] <- system.time(
          cs <- simLocal.spatialProcess(krigObj_temp, 
                                        list(x = pred_grid_side, y = pred_grid_side),
                                        M = 10, extrap = TRUE,
                                        fast = TRUE, NNSize = 4,
                                        NNSizePredict = 4)
        )[3]
      }
      sim_CSTiming_rapid[j, i] <- median(t_temp)
    }
  }

  save(sim_CSTiming_exact, sim_CSTiming_fast, sim_CSTiming_rapid, file = "SimCSTiming_New.RData")
} else {
  load("SimCSTiming.RData")
}



pdf("../Plots/Fig07_SimCSTiming.pdf", width=8, height=3)

ylim_axis <- c(2, 153)

# set layout and outer margin (oma) to leave space on the right
par(mfrow = c(1, 3), oma = c(0, 0, 0, 11), mar = c(4, 4, 2, 1), xpd = NA)
plot(x = mList[4:10]^2, y = sim_CSTiming_exact[4:10, 1], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, cex.lab = 1.3,
     ylab = "Time [s]", xlab = "", col = goldenTech, 
     ylim = ylim_axis)
     #ylim = c(0.001, max(c(sim_predTiming_exact[,1]))))
side_labels <- paste0(mList[4:10], " x ", mList[4:10])
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:10]^2, sim_CSTiming_fast[4:10, 1], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = lightGray)
lines(mList[4:10]^2, sim_CSTiming_rapid[4:10, 1], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = blasterBlue)
mtext("(a)", line=.5, adj=.05)

par(mar = c(4, 3, 2, 2))
plot(x = mList[4:10]^2, y = sim_CSTiming_exact[4:10, 2], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, cex.lab = 1.3,
     xlab = "Prediction Grid Size", ylab = "", col = goldenTech,
     ylim = ylim_axis)
     #ylim = c(0.005, max(sim_predTiming_exact[,2])))
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:7]^2, sim_CSTiming_fast[4:7, 2], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = lightGray)
lines(mList[4:10]^2, sim_CSTiming_rapid[4:10, 2], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = blasterBlue)
mtext("(b)", line=.5, adj=.05)

#par(mar=c(4,4,2,10))
plot(x = mList[4:10]^2, y = sim_CSTiming_exact[4:10, 3], log = 'xy',
     xaxt = "n", type = 'o', pch = 19, lwd = 3, 
     xlab = "", ylab = "", col = goldenTech,
     ylim = ylim_axis)
     #ylim = c(0.03, 16))
axis(1, at = mList[4:10]^2, labels = side_labels)
lines(mList[4:10]^2, sim_CSTiming_fast[4:10, 3], type = 'o', pch = 19, 
      lwd = 3, lty = 3, col = lightGray)
lines(mList[4:10]^2, sim_CSTiming_rapid[4:10, 3], type = 'o', pch = 19, 
      lwd = 3, lty = 1, col = blasterBlue)
mtext("(c)", line=.5, adj=.05)

legend("topleft", lty = 1,           # line types
       lwd = 3,                       # line width
       pch = c(19, 19, 19),       # point symbols
       pt.cex = 1,                    # point size
       col = c(goldenTech, lightGray, blasterBlue),
       legend = c("CS w. Exact Pred.", "Fast CS w. Exact Pred.", 
                  "Fast CS w. Rapid Pred."),
       inset = c(1.05, 0.3), xpd = NA, bty = "n", cex = 1.1)

dev.off()