patternCausality <- function(X,Y,E,tau,metric,h,weighted) {
  ###################################
  ### STEP 0: PREPARATORY ACTIONS ###
  ###################################
  NNSPAN = E+1 # Former NN | Reserves a minimum number of nearest neighbors
  CCSPAN = (E-1)*tau # This will remove the common coordinate NNs
  hashedpatterns <- patternHashing(E)
  #####################################
  ### STEP 1: THE SHADOW ATTRACTORS ###
  #####################################
  #= [A] =# State Space
  Mx <- stateSpace(X,E,tau)
  My <- stateSpace(Y,E,tau)
  #= [B] =# Signature Space
  SMx <- signatureSpace(Mx,E)
  SMy <- signatureSpace(My,E)
  #= [C] =# Pattern Space
  PSMx <- patternSpace(SMx,E)
  PSMy <- patternSpace(SMy,E)
  #= [D] =# Distance Matrix | First row corresponds to t=1
  Dx <- distanceMatrix(Mx,metric)
  Dy <- distanceMatrix(My,metric)
  #= Check whether time series length is sufficient
  FCP <- firstCausalityPoint(E,tau,h,X)
  #= Calculate the main loop duration of the algorithm
  al_loop_dur <- FCP:(length(X)-(E-1)*tau-h)
  #= Calculate the loop duration for out of sample forecasts
  out_of_sample_loop_dur <- ((length(X)-(E-1)*tau-h)+1):nrow(Mx)
  #= KEEPING THE PC MATRICES | Causality is considered only from FCP onwards
  predictedPCMatrix <- dataBank(type = "array",dimensions=c(3^(E-1),3^(E-1),length(Y)))
  #pb <- tkProgressBar(title = "Deploying PC Mk. II", min = 0,
  #                    max = length(al_loop_dur), width = 500)
  for(i in al_loop_dur) {
    if (!anyNA(c(Mx[i,],My[i+h,]))) {
      ###################################################################
      ### STEP 2: The Nearest Neighbours and their Future projections ###
      ###################################################################
      NNx <- pastNNsInfo(CCSPAN,NNSPAN,Mx,Dx,SMx,PSMx,i,h)
      if (!anyNA(Dy[i,NNx$times+h])) {
        projNNy <- projectedNNsInfo(My,Dy,SMy,PSMy,NNx$times,i,h)
        #######################################################################
        ### STEP 3: The affected variable's predicted pattern h steps ahead ###
        #######################################################################
        predictedSignatureY <- predictionY(E,projNNy,zeroTolerance=E-1)$predictedSignatureY
        predictedPatternY <- predictionY(E,projNNy,zeroTolerance=E-1)$predictedPatternY[1]
        #############################################
        ### STEP 4: The causal variable's pattern ###
        #############################################
        ####################signatureX <- signaE(E,SignX,i)
        signatureX <- SMx[i,]
        patternX <- PSMx[i,]
        ####################################################
        ### STEP 5: The affected variable's real pattern ###
        ####################################################
        #######realSignatureY <- signaE(E,SignY,(i+h))
        realSignatureY <- SMy[(i+h),]
        realPatternY <- PSMy[i+h]
        ##########################################################################
        ### STEP 6: The nature and intensity of causality at every time step t ###
        ##########################################################################
        pc <- fillPCMatrix(weighted,predictedPatternY,realPatternY,predictedSignatureY,realSignatureY,patternX,signatureX)
        predictedPCMatrix[which(hashedpatterns==patternX),which(hashedpatterns==predictedPatternY),i] <- pc$predicted
      }
    }
    #setTkProgressBar(pb, i, label=paste( i/al_loop_dur[length(al_loop_dur)], 0),"% PC Mk. II In-Sample Assignment Completion")
  }
  causality <- natureOfCausality(predictedPCMatrix,al_loop_dur,hashedpatterns,X)
  totalCausPercent <- 1-mean(causality$noCausality[al_loop_dur],na.rm = T)
  posiCausPercent <- mean(ifelse(causality$noCausality[al_loop_dur]!=1,causality$Positive[al_loop_dur],NA),na.rm = T)
  negaCausPercent <- mean(ifelse(causality$noCausality[al_loop_dur]!=1,causality$Negative[al_loop_dur],NA),na.rm = T)
  darkCausPercent <- mean(ifelse(causality$noCausality[al_loop_dur]!=1,causality$Dark[al_loop_dur],NA),na.rm = T)
  #return(list(causality,totalCausPercent,posiCausPercent,negaCausPercent,darkCausPercent))
  return(data.frame(total=totalCausPercent,positive=posiCausPercent,negative=negaCausPercent,dark=darkCausPercent))
}