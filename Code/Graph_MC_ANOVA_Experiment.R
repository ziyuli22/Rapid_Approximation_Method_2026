# load the newest fields library
library(fields)
# make sure to set working directory to the location of this file by
# Session --> Set Working Directory --> To Source File Location
###########################

n_obs <- c(200, 500, 1500) # n obs
corr_distance <- c(0.2, 0.4, 0.8) 
smoothness <- c(0.5, 1, 1.5) # actual smoothness fixed at 1
nugget <- c(0.01, 0.1, 0.5) # actual nugget = 0.45
NN <- c(2, 4, 8) # nearest neighbor
nGrid <- c(100, 350, 500) # grid

draws <- 50 # monte carlo draws

target <- rbind(
  c(0, 1), #corner
  c(0.25, 0.75), # mid-diag
  c(0.5, 0.5), # middle
  c(0, 0.5), # edge
  c(0.25, 0.5) # mid-left
)
# If RUN = TRUE, takes about 2 hrs
RUN = FALSE

if (RUN == TRUE){
  data_X <- expand.grid(
    corrDistance = corr_distance,
    smoothness = smoothness,
    nugget = nugget,
    NN = NN,
    nObs = n_obs,
    nGrid = nGrid
  )
  
  data <- data.frame(data_X)
  data$corner_mean <- NA
  data$mid_diag_mean <- NA
  data$center_mean <- NA
  data$edge_mean <- NA
  data$mid_center_mean <- NA
  data$corner_se <- NA
  data$mid_diag_se <- NA
  data$center_se <- NA
  data$edge_se <- NA
  data$mid_center_se <- NA
  data$corner_error_sd <- NA
  data$mid_error_sd <- NA
  data$center_error_sd <- NA
  data$edge_error_sd <- NA
  data$mid_center_error_sd <- NA
  
  L1_data <- data[, c(7:11, 17:ncol(data))]
  
  set.seed(123)
  
  for (i in 1:length(nGrid)) {
    nGrid_temp <- nGrid[i]
    pred_grid_side_temp <- seq(0, 1, length.out = nGrid_temp)
    gridList_temp <- list(x = pred_grid_side_temp, y = pred_grid_side_temp)
    
    grid_temp <- make.surface.grid(gridList_temp)
    d_to_target_temp <- rdist(grid_temp, target)
    id_target_temp <- apply(d_to_target_temp, 2, which.min) 
    target_temp <- grid_temp[id_target_temp, ]
    
    counter = 0
    for (j in 440:nrow(data)) {
      counter = counter + 1
      
      smoothness_temp <- data$smoothness[j]
      aRange_temp <- Matern.cor.to.range(data$corrDistance[j], smoothness_temp, 0.7)
      nugget_temp <- data$nugget[j]
      NN_temp <- data$NN[j]
      nObs_temp <- data$nObs[j]
      if (data$nGrid[j] == nGrid_temp){
        
        cases_temp <- matrix(NA, nrow = draws, ncol = 5)
        SE_temp <- matrix(NA, nrow = draws, ncol = 5)
        L1_temp <- matrix(NA, nrow = draws, ncol = 5)
        for(k in 1:draws){
          s0_temp = cbind(x = runif(nObs_temp), y = runif(nObs_temp))
          dist.mat_temp <- rdist(s0_temp,s0_temp)
          SigmaCov_temp <- Matern(dist.mat_temp, range = aRange_temp, smoothness = smoothness_temp) # Covariance matrix (distance dependent)
          Sigma.c_temp <- chol(SigmaCov_temp) # Cholesky factorization of the covariance matrix
          y0_temp <- t(Sigma.c_temp) %*% rnorm(nObs_temp) + rnorm(nObs_temp, mean = 0, sd = sqrt(nugget_temp))
          
          # fitting model
          fit_temp <- mKrig(s0_temp, y0_temp, 
                       Covariance = "Matern", cov.args=list(
                         aRange = aRange_temp, nu = smoothness_temp), 
                       lambda = nugget_temp,
                       sigma2 = 1.0, simpleKriging = TRUE
                       )
          exact_pred_temp <- predict(fit_temp, target_temp)
          exact_SE_temp <- predictSE.mKrig(fit_temp, target_temp)
          
          approx_pred_hold<- predictSurface(fit_temp, gridList = gridList_temp, 
                                            fast=TRUE, extrap=TRUE, NNSize = NN_temp)
          approx_pred_temp<- approx_pred_hold$z[id_target_temp]
          
          cases_temp[k, ] <- (exact_pred_temp - approx_pred_temp)
          L1_temp[k, ] <- abs(exact_pred_temp - approx_pred_temp)
          SE_temp[k, ] <- exact_SE_temp # exact SE
        }
        
        point_means <- apply(cases_temp, 2, mean)
        point_SE <- apply(SE_temp, 2, mean)
        point_error_sd <- apply(cases_temp, 2, sd)
        point_L1_means <- apply(L1_temp, 2, mean)
        point_L1_sd <- apply(L1_temp, 2, sd)
        
        data$corner_mean[j] <- point_means[1]
        data$mid_diag_mean[j] <- point_means[2]
        data$center_mean[j] <- point_means[3]
        data$edge_mean[j] <- point_means[4]
        data$mid_center_mean[j] <- point_means[5]
        data$corner_se[j] <- point_SE[1]
        data$mid_diag_se[j] <- point_SE[2]
        data$center_se[j] <- point_SE[3]
        data$edge_se[j] <- point_SE[4]
        data$mid_center_se[j] <- point_SE[5]
        data$corner_error_sd[j] <- point_error_sd[1]
        data$mid_error_sd[j] <- point_error_sd[2]
        data$center_error_sd[j] <- point_error_sd[3]
        data$edge_error_sd[j] <- point_error_sd[4]
        data$mid_center_error_sd[j] <- point_error_sd[5]
        
        L1_data$corner_mean[j] <- point_L1_means[1]
        L1_data$mid_diag_mean[j] <- point_L1_means[2]
        L1_data$center_mean[j] <- point_L1_means[3]
        L1_data$edge_mean[j] <- point_L1_means[4]
        L1_data$mid_center_mean[j] <- point_L1_means[5]
        L1_data$corner_error_sd[j] <- point_L1_sd[1]
        L1_data$mid_error_sd[j] <- point_L1_sd[2]
        L1_data$center_error_sd[j] <- point_L1_sd[3]
        L1_data$edge_error_sd[j] <- point_L1_sd[4]
        L1_data$mid_center_error_sd[j] <- point_L1_sd[5]
        
        if (counter %% 20 == 0) {
          timestamp()
          cat("nGrid =", nGrid_temp,
              "Processed", counter, "rows\n")
        }
      }
      
    }
    timestamp()
    cat("done nGrid =", nGrid_temp, "\n")
  }
  save(data, L1_data, file = "ANOVA_MC_Study_New.RData")
} else {
  load("ANOVA_MC_Study.RData")
}

data$corrDistance<- as.factor(data$corrDistance)
data$smoothness <- as.factor(data$smoothness)
data$nugget <- as.factor(data$nugget)
data$NN <- as.factor(data$NN)
data$nObs <- as.factor(data$nObs)
data$nGrid <- as.factor(data$nGrid)

### Corner L1 mean
lmod_corner <- lm(log10(L1_data$corner_mean) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
summary(lmod_corner)
anova(lmod_corner)

lmod_small_corner <- lm(log10(L1_data$corner_mean) ~ smoothness + NN + nGrid + nugget
           + smoothness:NN, data = data)
summary(lmod_small_corner)
anova(lmod_small_corner)

### Mid Diag L1 Mean
lmod_mid_diag <- lm(log10(L1_data$mid_diag_mean) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
summary(lmod_mid_diag)
anova(lmod_mid_diag)

### Center L1 Mean
lmod_center <- lm(log10(L1_data$center_mean) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
summary(lmod_center)
anova(lmod_center)

lmod_small_center <- lm(log10(L1_data$center_mean) ~ smoothness + NN + nGrid + nugget
                        + smoothness:NN, data = data)
summary(lmod_small_center)
anova(lmod_small_center)

### Edge L1 Mean
lmod_edge<- lm(log10(L1_data$edge_mean) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
anova(lmod_edge)

### Mid Center L1 Mean
lmod_mid_center <- lm(log10(L1_data$mid_center_mean) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
anova(lmod_mid_center)

####
### Corner L1 SD
lmod_corner_SD <- lm(log10(L1_data$corner_error_sd) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
anova(lmod_corner_SD)

### Mid Diag L1 SD
lmod_mid_diag_SD <- lm(log10(L1_data$mid_error_sd) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
anova(lmod_mid_diag_SD)

### Edge L1 Mean
lmod_edge_SD <- lm(log10(L1_data$edge_error_sd) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
anova(lmod_edge_SD)

# Combining 5 point's means
L1_mean_error <- apply(L1_data[, 1:5], 1, mean)
#lmod_averaged_error <- lm(log10(L1_mean_error) ~ (corrDistance + smoothness + nugget + NN + nGrid + nObs)^2, data = data)
lmod_averaged_error <- lm(log10(L1_mean_error) ~ corrDistance + smoothness + nugget + NN + nGrid + nObs, data = data)
summary(lmod_averaged_error)
anova(lmod_averaged_error)

lmod_averaged_error_small <- lm(log10(L1_mean_error) ~ smoothness + NN + nGrid + nugget
                                + smoothness:NN, data = data)
summary(lmod_averaged_error_small)

pdf("../Plots/Fig03_InteractionPlot.pdf", width=7, height=4.5)

y_lim = c(-8.5, -2.5)
legend_cex = 1.4
main_cex = 1.4
axis_cex = 1.2
tick_cex = 1.2

par(mfrow = c(2, 3))

# NN (x-axis) by corrDist (lines)
lev <- levels(factor(data$NN))
par(mar = c(4.5,4.5,3,0))  # extra right margin + allow drawing outside
interaction.plot(data$corrDistance, data$NN, log10(L1_data$center_mean),
                 main = expression("Correlation Dist. vs L"),
                 xlab = expression("Correlation Dist."), 
                 ylab = expression("Log Abs Error"),
                 ylim = y_lim,
                 trace.label = "",            # no internal legend title
                 legend = FALSE,              # <- turn off the built-in legend
                 lty = seq_along(lev), lwd = 2, 
                 cex.main = main_cex, cex.lab = axis_cex, cex.axis = tick_cex)


# nGrid  (x-axis) by corrDistance (lines)
interaction.plot(data$nGrid, data$NN, log10(L1_data$center_mean),
                 main = expression( "Grid Size vs L"), 
                 xlab = expression("Grid Size"), ylab = "",
                 ylim = y_lim, yaxt = "n",
                 trace.label = "",            # no internal legend title
                 legend = FALSE,              # <- turn off the built-in legend
                 lty = seq_along(lev), lwd = 2, 
                 cex.main = main_cex, cex.lab = axis_cex, cex.axis = tick_cex)

par(mar = c(2,0,0,0))
plot.new()
legend("center", title = expression("Nearest Neighbor (L)"),
       legend = lev, lty = seq_along(lev), bty = "n",
       cex = legend_cex, lwd = 2)

# corrDist (x-axis) by nGrid (lines)
lev <- levels(factor(data$smoothness))
par(mar = c(4.5,4.5,3,0))
interaction.plot(data$NN, data$smoothness,log10(L1_data$center_mean),
                 main = expression( "L vs " * nu), 
                 xlab = expression("Nearest Neighbor (L)"), 
                 ylab = "Log Abs Error",
                 ylim = y_lim,
                 trace.label = "",            # no internal legend title
                 legend = FALSE,              # <- turn off the built-in legend
                 lty = seq_along(lev), lwd = 2, 
                 cex.main = main_cex, cex.lab = axis_cex, cex.axis = tick_cex)


par(mar = c(2,0,0,0))
plot.new()
legend("center", title = expression("Smoothness (" * nu * ")"),
       legend = lev, lty = seq_along(lev), bty = "n",
       cex = legend_cex, lwd = 2)

dev.off()
