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

# grid for prediction
m <- 128
xGrid <- seq(from = min(s[, 1]), to = max(s[, 1]), length.out = 0.5*m)
yGrid <- seq(from = min(s[, 2]), to = max(s[, 2]), length.out = m)
predGrid <- make.surface.grid(list(x1 = xGrid, x2 = yGrid))

# fit the process
krigObj <- spatialProcess(s, precip, smoothness = 1.5)

# num of ensembles
ensemble <- 100

if (RUN){
  
  set.seed(123)
  # original CS
  originalCS <- sim.spatialProcess(krigObj, predGrid, 
                                   M = ensemble)
  
  # Maggie's CS
  set.seed(123)
  fastCS <- simLocal.spatialProcess(krigObj, 
                                    list(x = xGrid, y = yGrid),
                                    extrap = TRUE,
                                    M = ensemble, NNSize = 5)
  
  # fast CS
  set.seed(123)
  rapidCS <- simLocal.spatialProcess(krigObj, 
                                     list(x = xGrid, y = yGrid),
                                     M = ensemble, extrap = TRUE,
                                     fast = TRUE, NNSize = 5,
                                     NNSizePredict = 4)
  
  save(originalCS, fastCS, rapidCS, file = "RainfallCSComparison_New.RData")
} else {
  load("RainfallCSComparison.RData")
}

SE <- predictSurfaceSE(krigObj, grid.list = list(x = xGrid, y = yGrid),
                       extrap=TRUE)


# calculate empirical SE, for sanity check, should agree with real SE
SE_original <- apply(originalCS, 1, sd)

# maggies' and with rapid prediction
SE_fast <- apply(fastCS$z, c(1,2), sd)

# CE & rapid prediction
SE_rapid <- apply(rapidCS$z, c(1, 2), sd)


# CE & rapid prediction
meanSurface_rapid <- apply(rapidCS$z, c(1, 2), mean)

stats(c(SE_original/c(SE$z)) )
stats(c(SE_fast/SE$z) )
stats(c(SE_rapid/SE$z) )

pdf("../Plots/Fig09_RainfallCS9inContour.pdf", width=8, height=4)

col_legend <- two.colors(start = blasterBlue, middle = 'white', 
                         end = redFlannel)

layout(matrix(1:3, nrow = 1), widths = c(1.9, 1.5, 2.2))

par(mar = c(4, 4, 2, 4.5))
bubblePlot(predGrid, meanSurface_rapid, xlab = "Lon", ylab = "Lat", 
           main = "", cex.lab = 1.5,
           noLegend = TRUE)
US(add = TRUE)
for (i in 1:100) {
  surfObj_temp <- as.surface(predGrid, rapidCS$z[, , i])
  contour(surfObj_temp$x, surfObj_temp$y, surfObj_temp$z,
          levels = c(9), add = TRUE, drawlabels = FALSE,
          col = alpha("white", 0.1), lwd = 1)
}
mtext("(a)", line=.5, adj=.05)

image.plot( add = TRUE, legend.only = TRUE, zlim = range(c(meanSurface_rapid)),
            col = viridisLite::viridis(256),
            horizontal = FALSE, legend.width = 2.5,
            legend.args = list(text = "inches", col = "black", 
                               cex = 0.8, line = 0.5))

diff <- SE$z-SE_rapid
absBound <- max(abs(max(diff)), abs(min(diff)))
par(mar = c(4, 2.5, 2, 6))
bubblePlot(predGrid, diff, xlab = "Lon", yaxt = 'n', ylab = "",
           main = "", zlim = c(-absBound, absBound), col = col_legend,
           cex.lab = 1.5, noLegend = TRUE)
US(add = TRUE)
mtext("(b)", line=.5, adj=.05)
image.plot( add = TRUE, legend.only = TRUE, zlim = c(-absBound, absBound),
            col = col_legend,
            horizontal = FALSE, legend.width = 2.5,
            legend.args = list(text = "inches", col = "black", 
                               cex = 0.8, line = 0.5))

par(mar = c(8, 5, 8, 1))
plot(c(SE$z), c(SE_rapid), main = "",
     xlab = "Actual SE", ylab = "Empirical SE", pch = 19, 
     col = alpha(lightGray, 0.2), cex.lab = 1.5)
abline(0, 1, col = redFlannel, lwd = 2)
mtext("(c)", line=.5, adj=.05)

dev.off()

