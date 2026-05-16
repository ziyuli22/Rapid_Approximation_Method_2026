library(fields)
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

# if RUN = TRUE
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

if (RUN){
  
  mList <- c(20, 40, 60, 100, 140, 200, 260, 350, 500)
  
  # fit the process
  krigObj <- spatialProcess(s, precip, smoothness = 1.5, aRange = 0.5)
  
  X <- cbind( 1, s)
  GpGpfit <- fit_model(precip, s, X, "matern_isotropic")
  
  # empty matrix to store timing info
  rainfallCSTiming <- data.frame(matrix(NA, nrow = length(mList), 
                                        ncol = 4))
  
  
  names(rainfallCSTiming) <- c("original", "Maggie",
                               "RapidNN4", "Vecchia")
  
  for (i in 1:length(mList)) {
    
    m <- mList[i]
    
    xGrid <- seq(from = min(s[, 1]), to = max(s[, 1]), length.out = m)
    yGrid <- seq(from = min(s[, 2]), to = max(s[, 2]), length.out = m)
    predGrid <- make.surface.grid(list(x1 = xGrid, x2 = yGrid))
    
    print(paste("The grid size is m =", mList[i]) )
    
    # Original
    if (m < 150){
      temp_t <- c()
      for (j in 1:10) {
        temp_t[j] <- system.time(
          originalCS <- sim.spatialProcess(krigObj, predGrid, M = 10)
        )[3]
      }
      rainfallCSTiming$original[i] <- median(temp_t)
    }
    
    # Maggie's
    temp_t <- c()
    for (j in 1:10) {
      temp_t[j] <- system.time(
        fastCS <- simLocal.spatialProcess(krigObj, 
                                          list(x = xGrid, y = yGrid),
                                          extrap = TRUE,
                                          M = 10, NNSize = 4)
      )[3]
    }
    rainfallCSTiming$Maggie[i] <- median(temp_t)
    

    # rapid method (NN = 4)
    temp_t <- c()
    for (j in 1:10) {
      temp_t[j] <- system.time(
        rapidCS <- simLocal.spatialProcess(krigObj, 
                                           list(x = xGrid, y = yGrid),
                                           M = 10, extrap = TRUE,
                                           fast = TRUE, NNSize = 4,
                                           NNSizePredict = 4)
      )[3]
    }
    rainfallCSTiming$RapidNN4[i] <- median(temp_t)
    
    
    X_pred <- cbind(1, predGrid)
    temp_t <- c()
    for (j in 1:10) {
      temp_t[j] <- system.time(
        vecchiaCS <- cond_sim(GpGpfit, predGrid, X_pred,
                              nsims = 10)
      )[3]
    }
    rainfallCSTiming$Vecchia[i] <- median(temp_t)
  }

  save(rainfallCSTiming, mList, file = "RainfallCSTiming_New.RData")
} else {
  load("RainfallCSTiming.RData")
}

pdf("../Plots/Fig11_RainfallCSTiming.pdf", width=9, height=4)

yLegend <- c(min(na.omit(rainfallCSTiming$original)),
             max(na.omit(rainfallCSTiming$Maggie)) + 20)
par( mar=c(4,4.5,2,12.5))
#yRange <- c(log10(min(c(rainfallCSTiming$original))), log10(max(c(rainfallCSTiming$Maggie))))
plot(x = mList[3:9]^2, y = rainfallCSTiming$original[3:9], type = 'o', log = 'xy',
     pch = 19, lwd = 3, main = "Rainfall Conditional Simulation Timing",
     xlab = "Number of Grid Points", ylab = "Time [Seconds]",
     xaxt = "n", cex.main = 1.5, cex.lab = 1.4, ylim = yLegend, col = goldenTech)
# Custom x-axis labels
side_labels <- paste0(mList[3:9], " x ", mList[3:9])
axis(1, at = mList[3:9]^2, labels = side_labels)

lines(x = mList[3:9]^2, y = rainfallCSTiming$Vecchia[3:9], type = 'o', 
      pch = 19, lwd = 3,
      col = mutedBlue)
lines(x = mList[3:9]^2, y = rainfallCSTiming$Maggie[3:9], type = 'o',
      pch = 19, lwd = 3,
      col = silver)
lines(x = mList[3:9]^2, y = rainfallCSTiming$RapidNN4[3:9], type = 'o', 
      pch = 19, lwd = 3,
      col = blasterBlue)

legend("topleft",lwd = 3,
       col = c(mutedBlue, goldenTech, silver, blasterBlue),
       lty = c(1, 1, 1, 1),
       legend = c("Vecchia (GpGp)", "CS via Exact Pred.", 
                  "Fast CS via Exact Pred.", 
                  "Fast CS via Rapid Pred."),
       inset=c(1.01,0.35), bty = "n",
       xpd=TRUE, cex = 1)

dev.off()
