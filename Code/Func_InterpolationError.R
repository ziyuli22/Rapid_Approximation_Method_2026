interpolation_error <- function(sObs,
                                input_gridList,
                                NNSize, 
                                aRange, 
                                smoothness){
  
  info<- augmentPredictionGrid( s=sObs, gridList=input_gridList, NNSize=NNSize)
  gridList<- info$predictionGrid
  obj<-approximateCovariance2D(sObs,
                               gridList=gridList,
                               NNSize=NNSize, 
                               Covariance="Matern", 
                               sigma2=1.0,
                               aRange=aRange,
                               covArgs= list( smoothness=smoothness),
                               verbose=FALSE
  )
  
  # this script only makes sense for one off grid point
  B<-spam2full( obj$B)
  sGridApprox<- make.surface.grid( obj$gridList)
  ind<- B!=0
  sGridApprox<- sGridApprox[ind,]
  
  sGrid <- list( x=seq( 0, 1, length.out = 427),
                 y=seq( 0, 1, length.out = 427))
  sGridPredict <- make.surface.grid( sGrid)
  
  Phi <- B[,ind]%*%Matern( rdist( sGridApprox,sGridPredict), aRange=aRange, smoothness=smoothness )
  truePhi <- Matern( rdist( sObs,sGridPredict), aRange=aRange, smoothness=smoothness )
  
  statsSummary <- max(abs(cbind( c( Phi-truePhi))))
  errorSurf <- Phi-truePhi
  return( list(maxErr = statsSummary, NNGrid = sGridApprox, finePredGrid = sGridPredict,
               error = errorSurf))
  
  # imagePlot( as.surface( sGrid, Phi-truePhi ))
  # points(sGridApprox, col="grey20", pch=16 )
  # points( sObs, col="magenta", pch=16)
  # points( make.surface.grid(input_gridList), col = "white", pch = 1)
  # 
}