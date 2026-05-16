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

# if RUN = TRUE, this might take 5min to run at m = 500
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

m = 500
xGrid <- seq(from = min(s[, 1]), to = max(s[, 1]), length.out = floor(0.5 * m))
yGrid <- seq(from = min(s[, 2]), to = max(s[, 2]), length.out = m)
predGrid <- make.surface.grid(list(x1 = xGrid, x2 = yGrid))

if (RUN){
  # fit the process
  krigObj <- spatialProcess(s, precip, smoothness = 1.5)
  
  # kriging prediction
  krigingPreds <- predictSurface(krigObj, gridList = list(x1 = xGrid, x2 = yGrid),
                                 extrap = TRUE)
  # our method
  rapidMethod <- predictSurface(krigObj, gridList = list(x1 = xGrid, x2 = yGrid),
                              fast = TRUE, NNSize = 4, extrap = TRUE)
  
  
  
  save(krigingPreds, rapidMethod, file = "RainfallPredictionComparison_New.RData")
} else {
  load("RainfallPredictionComparison.RData")
}
  
# calculating relative error
relError <- matrix(NA, nrow = nrow(krigingPreds$z), 
                   ncol = ncol(krigingPreds$z))
for (i in 1:nrow(krigingPreds$z)) {
  for (j in 1:ncol(krigingPreds$z)){
    if (krigingPreds$z[i,j] == 0){
      relError[i,j] <- 0
    } 
    else{
      relError[i,j] <- (krigingPreds$z[i,j] - rapidMethod$z[i,j])/krigingPreds$z[i,j]
    }
  }
}

# absolute error
absError <- matrix(NA, nrow = nrow(krigingPreds$z), 
                   ncol = ncol(krigingPreds$z))
for (i in 1:nrow(krigingPreds$z)) {
  for (j in 1:ncol(krigingPreds$z)){
    absError[i, j] <- krigingPreds$z[i,j] - rapidMethod$z[i,j]
  }
}

pdf("../Plots/Fig08_RainfallPredictionComparison.pdf", width=7, height=3.5)

zLegend <- c(min(c(krigingPreds$z, rapidMethod$z, precip)),
             max(c(krigingPreds$z, rapidMethod$z, precip))
)

#par(mfrow = c(1, 4), oma = c(0, 0, 0, 9), mar = c(4, 4, 2, 1), xpd = NA)

layout(matrix(1:4, nrow = 1), widths = c(1.15, 0.75, 1.15, 1.3))
par(oma = c(0,0,0,1))

zLegend <- range(c(krigingPreds$z, rapidMethod$z, precip))

# Plot 1: Observations
par(mar = c(4, 5, 3, 0))
bubblePlot(s, precip, zlim = zLegend, main = "",
           noLegend = TRUE, xlab = "", ylab = "Lat",
           cex.lab = 1.2)
US(add = TRUE)
mtext("(a)", line=.5, adj=.05)

# kriging prediction
par(mar = c(4, 0.7, 3, 0))
bubblePlot(predGrid, krigingPreds$z, zlim = zLegend, main = "", 
           noLegend = TRUE, xlab = "Lon", ylab = "", yaxt = "n",
           cex.lab = 1.2)
US(add = TRUE)
# Convert predGrid and z to surface for contour
surfObj <- as.surface(predGrid, krigingPreds$z)
contour(surfObj$x, surfObj$y, surfObj$z,
        levels = c(9), add = TRUE, drawlabels = TRUE,
        col = "white", lwd = 1)
mtext("(b)", line=.5, adj=.05)

# rapid prediction
par(mar = c(4, 0.7, 3, 5))
bubblePlot(predGrid, rapidMethod$z, zlim = zLegend,main = "",
           xlab = "", ylab = "", yaxt = "n", noLegend = TRUE,
           cex.lab = 1.2)
US(add = TRUE)
# Convert predGrid and z to surface for contour
surfObj2 <- as.surface(predGrid, rapidMethod$z)
contour(surfObj2$x, surfObj2$y, surfObj2$z,
        levels = c(9), add = TRUE, drawlabels = TRUE,
        col = "white", lwd = 1)
mtext("(c)", line=.5, adj=.05)

# Add legend
image.plot( add = TRUE, legend.only = TRUE, zlim = zLegend,
            col = viridisLite::viridis(256),
            horizontal = FALSE, legend.width = 2.5,
            legend.args = list(text = "inches", col = "black", 
                               cex = 0.75, line = 0.5))

# Absolute Difference
col_legend <- two.colors(start = blasterBlue, middle = 'white', 
                         end = redFlannel)
zLegend_diff <- c(min(c(absError)), - min(c(absError)))
par(mar = c(4, 0, 3, 5))
bubblePlot(predGrid, absError, main = "",
           xlab = "", ylab = "", yaxt = "n", noLegend = TRUE,
           cex.main = 2, cex.lab = 2, zlim = zLegend_diff,
           col = col_legend)
US(add = TRUE)
mtext("(d)", line=.5, adj=.05)

# Add legend
image.plot( add = TRUE, legend.only = TRUE, zlim = zLegend_diff,
            col = col_legend,
            horizontal = FALSE, legend.width = 2.5,
            legend.args = list(text = "inches", col = 'black', 
                               cex = 0.75, line = 0.5))

dev.off()
