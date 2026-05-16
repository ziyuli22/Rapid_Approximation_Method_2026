library(fields)
library(LatticeKrig)
library(GpGp)

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

# if RUN = TRUE this takes a while to run
RUN = FALSE

### 
data("NorthAmericanRainfall2")
s0 <- cbind(NorthAmericanRainfall2$longitude,
            NorthAmericanRainfall2$latitude)
precip0 <- NorthAmericanRainfall2$precip/( 100*2.54)
subset100 <- (s0[, 1] > - 100- 5) & (s0[, 1] < - 100 + 8) &
  (s0[, 2] >= 27) & (s0[, 2] <= 55)
s <- s0[subset100, ]
precip <- precip0[subset100]

mList <- c(60, 100, 140, 200, 260, 350, 500)

if (RUN){
  
  # empty matrix to store timing info, first row is fitting
  methodPredictionTiming <- data.frame(matrix(NA, nrow = length(mList) + 1, 
                                        ncol = 6))
  names(methodPredictionTiming) <- c("Exact", "RapidNN2", "RapidNN4", "RapidNN8", "Lattice", "Vecchia")
  
  # fit the processes
  # for Rapid Method, takes about 6min if we don't restrict nu = 1.5
  predTiming <- c()
  for (i in 1:20) {
    predTiming[i] <- system.time(
      krigObj <- spatialProcess(s, precip, gridN = 15)
      )[3]
  }
  methodPredictionTiming$RapidNN4[1] <- median(predTiming)
  methodPredictionTiming$RapidNN2[1] <- median(predTiming)
  methodPredictionTiming$Exact[1] <- median(predTiming)
  methodPredictionTiming$RapidNN8[1] <- median(predTiming)
  # for Lkrig, just about 1 min
  predTiming <- c()
  for (i in 1:20) {
    predTiming[i] <- system.time(
      LKrig_Obj <- LatticeKrig( s, precip, NC=20,
                                findAwght = TRUE,
                                BasisFunction = "WendlandFunction",
                                BasisType = "Radial",
                                normalize=FALSE) 
      )[3]
  }
  methodPredictionTiming$Lattice[1] <- median(predTiming)
  # for Vecchia, a little more than lattice. 
  predTiming <- c()
  X <- cbind( 1, s)
  for (i in 1:20) {
    predTiming[i] <- system.time(
      GpGpfit <- fit_model(precip, s, X, "matern_isotropic")
    )[3]
  }
  methodPredictionTiming$Vecchia[1] <- median(predTiming)

  
  for (i in 1:length(mList)) {
    
    m <- mList[i]
    
    xGrid <- seq(from = min(s[, 1]), to = max(s[, 1]), length.out = m)
    yGrid <- seq(from = min(s[, 2]), to = max(s[, 2]), length.out = m)
    predGrid <- make.surface.grid(list(x1 = xGrid, x2 = yGrid))
    
    t_temp_exact <- c()
    for (j in 1:20) {
      t_temp_exact[j] <- system.time(
        pred_exact <- predictSurface(krigObj, 
                                     gridList = list(x1 = xGrid, x2 = yGrid),
                                     extrap = TRUE)
      )[3]
    }
    methodPredictionTiming$Exact[i + 1] <- median(t_temp_exact)
    
    t_temp_rapid2 <- c()
    for (j in 1:20) {
      t_temp_rapid2[j] <- system.time(
        pred_rapid <- predictSurface(krigObj, 
                                     gridList = list(x1 = xGrid, x2 = yGrid),
                                     extrap = TRUE, fast = TRUE, NNSize = 2)
      )[3]
    }
    methodPredictionTiming$RapidNN2[i + 1] <- median(t_temp_rapid2)
    
    t_temp_rapid4 <- c()
    for (j in 1:20) {
      t_temp_rapid4[j] <- system.time(
        pred_rapid <- predictSurface(krigObj, 
                                     gridList = list(x1 = xGrid, x2 = yGrid),
                                     extrap = TRUE, fast = TRUE, NNSize = 4)
      )[3]
    }
    methodPredictionTiming$RapidNN4[i + 1] <- median(t_temp_rapid4)
    
    t_temp_rapid8 <- c()
    for (j in 1:20) {
      t_temp_rapid8[j] <- system.time(
        pred_rapid <- predictSurface(krigObj, 
                                     gridList = list(x1 = xGrid, x2 = yGrid),
                                     extrap = TRUE, fast = TRUE, NNSize = 8)
      )[3]
    }
    methodPredictionTiming$RapidNN8[i + 1] <- median(t_temp_rapid8)
    
    t_temp_lattice <- c()
    for (j in 1:20){
      t_temp_lattice[j] <- system.time(
        pred_LKrig <- predictSurface(LKrig_Obj, 
                                     grid.list = list(x1 = xGrid, x2 = yGrid),
                                     extrap = TRUE)
      )[3]
    }
    methodPredictionTiming$Lattice[i + 1] <- median(t_temp_lattice)
    
    t_temp_vecchia <- c()
    for (j in 1:20) {
      X_pred <- cbind( 1, predGrid)
      t_temp_vecchia <- system.time(
        pred_vecchia <- predictions(GpGpfit, predGrid, X_pred)
      )[3]
    }
    methodPredictionTiming$Vecchia[i + 1] <- median(t_temp_vecchia)
    
    print(paste("The grid size is m =", mList[i]) )
    
  }

  save(methodPredictionTiming, mList, file = "RainfallPredictionMethodTiming_New.RData")
} else {
  load("RainfallPredictionMethodTiming.RData")
}


pdf("../Plots/Fig10_RainfallMethodTiming.pdf", width=9, height=4)

len_data <- nrow(methodPredictionTiming)
yLegend <- c(min(na.omit(methodPredictionTiming$RapidNN2)),
             max(na.omit(methodPredictionTiming$Vecchia)) + 20)
par( mar=c(4,4.5,2,11.5))
#yRange <- c(log10(min(c(rainfallCSTiming$original))), log10(max(c(rainfallCSTiming$Maggie))))
plot(x = mList^2, y = methodPredictionTiming$Vecchia[2:len_data], type = 'o', log = 'xy',
     pch = 19, lwd = 3, main = "Rainfall Prediction Timing",
     xlab = "Number of Grid Points", ylab = "Time [Seconds]",
     xaxt = "n", cex.main = 1.5, cex.lab = 1.4, ylim = yLegend, col = mutedBlue)
# Custom x-axis labels
side_labels <- paste0(mList, " x ", mList)
axis(1, at = mList^2, labels = side_labels)

lines(x = mList^2, y = methodPredictionTiming$Exact[2:len_data], type = 'o',
      pch = 19, lwd = 3,
      col = goldenTech)

lines(x = mList^2, y = methodPredictionTiming$Lattice[2:len_data], type = 'o',
      pch = 19, lwd = 3,
      col = silver)
lines(x = mList^2, y = methodPredictionTiming$RapidNN4[2:len_data], type = 'o', 
      pch = 19, lwd = 3,
      col = blasterBlue)
lines(x = mList^2, y = methodPredictionTiming$RapidNN2[2:len_data], type = 'o', 
      pch = 19, lwd = 3, lty = 2, 
      col = blasterBlue)
lines(x = mList^2, y = methodPredictionTiming$RapidNN8[2:len_data], type = 'o', 
      pch = 19, lwd = 3, lty = 3, 
      col = blasterBlue)

legend("topleft",lwd = 3,
       col = c( mutedBlue, goldenTech, silver, blasterBlue, blasterBlue, blasterBlue),
       lty = c(1, 1, 1, 3, 1, 2),
       legend = c("Vecchia (GpGp)", "Exact Prediction", "Lattice Krig",  "Rapid Pred. (L = 8)",
                  "Rapid Pred. (L = 4)", "Rapid Pred. (L = 2)"),
       inset=c(1.01,0.2), bty = "n",
       xpd=TRUE, cex = 1)

dev.off()
