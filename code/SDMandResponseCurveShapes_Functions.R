

##logit and inverse logit

logit <- function(x) log(x / (1 - x))
inv.logit <- function(x) (1 + exp(-x))^-1

#inverse of rho matrix
inv.rmat <- function(n,rho=0.9){
  rinv <- diag((1 + rho^2),n,n)
  rinv[1,1] <- 1
  rinv[n,n] <- 1
  rinv[row(rinv) == (col(rinv)-1)] <- -rho
  rinv[row(rinv) == (col(rinv)+1)] <- -rho
  return(rinv)
}  	

#calculate R from covariance matrix
rmat  <- function(n, rho = 0.9) {
  mat <- diag(rep(1,n))
  mat <- rho^abs(row(mat)-col(mat))
  ((1 - rho^2)^-1)*mat
}

## Functions to read specimen data
calc_ObservationsFromProbabilities <- function(probs, N, VAR, sigma,w){###need to incorporate long and Lat
  obs <- NULL
  
  #no residual variance, just a straight binomial draw
  if(VAR == "binom")
    obs <- t(sapply(probs, FUN=function(x) rbinom(n=N, size=1, prob=x))) #N Bernoulli draws
  
  #binomal draw with a residual variance in logit space such that at p = 0.5, sd = +/- 0.10
  if(VAR == "binom+res"){
    prob_new <-  t(sapply(probs, FUN=function(x) inv.logit(rnorm(N, logit(x), sigma))))#add variability
    obs <- matrix(rbinom(length(prob_new),size = 1, prob = prob_new),nrow=nrow(prob_new), ncol = ncol(prob_new))
  }
  
  if(VAR == "spatial"){
    prob_new <- inv.logit(logit(probs)+2 * ((w-min(w))/(max(w) - min(w)) - .5) )
    obs <- t(sapply(prob_new, FUN=function(x) rbinom(n=N, size=1, prob=x))) #N Bernoulli draws
  }
  
  if(VAR == "spatial+res"){
    prob_new <- inv.logit(logit(probs)+2 * ((w-min(w))/(max(w) - min(w)) - .5) )
    prob_new <-  t(sapply(prob_new, FUN=function(x) inv.logit(rnorm(N, logit(x), sigma))))#add variability
    obs <- matrix(rbinom(length(prob_new),size = 1, prob = prob_new),nrow=nrow(prob_new), ncol = ncol(prob_new))
  }
  
  return(obs=obs)
}


get_TypeData_FromFile <- function(type, center=TRUE, centerBasedOnRegionIDs=2){
  dat <- read.csv(file=file.path(dir.in, paste0(type, ".csv")))
  dat <- dat[, c("region", "x", "y", "LnP", "MinT", "prob")]
  
  if(center){###DMB revision: scaled variables need to be centered with the same mean for all regions
    mean_dat <- apply(dat[dat$region %in% centerBasedOnRegionIDs, c("LnP", "MinT")], 2, mean) 
    centerMeans <- c(LnPmean=mean(mean_dat["LnP"]), MinTmean=mean(mean_dat["MinT"]))
    dat$LnPScaled <- dat$LnP - centerMeans[1]
    dat$MinTScaled <- dat$MinT - centerMeans[2]
  }
  
  return(list(dat=dat, centerMeans=centerMeans))
}

get_GeographicRaster_FromFile <- function(region){
  raster(file.path(dir.gis, paste0("extent.raster.r", region, ".grd")))
}

set_Data <- function(type, error, obs, dat, samp, run){
  
  prepare_Data <- function(obs, dat, samp){
    return(list(resp.var=as.numeric(obs[samp]),
                resp.xy=dat[samp, c("x", "y")],
                expl.var=dat[samp, c("LnP", "LnPScaled", "MinT", "MinTScaled")]))
  }
  
  #calibration dataset
  blist <- prepare_Data(obs=obs, dat=dat, samp=samp)
  
  #evaluation dataset
  tmp.samp <- (1:length(obs))[-samp]
  if(equalSamples)
    tmp.samp <- get.balanced.sample(obsData[[paste(type,error,sep="_")]][ibase, run],tmp.samp)
  
  elist <- prepare_Data(obs=obs, dat=dat, samp=tmp.samp)
  
  bdat <- list(resp.name=type, 
               resp.var=blist$resp.var, eval.resp.var=elist$resp.var, #presence/absence
               resp.xy=blist$resp.xy,   eval.resp.xy=elist$resp.xy,   #xy
               expl.var=blist$expl.var, eval.expl.var=elist$expl.var) #covariates
  
  return(bdat = bdat)
}


## Functions to build the SDMs
set_options <- function(model, level=c("linear", "squared", "interaction")){
  level <- match.arg(arg=level, choices=c("linear", "squared", "interaction"), several.ok=TRUE)
  
  if("GLM" %in% model){ #GLM options
    term <- NULL
    if(any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
    if(any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
    if(any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formGLM <- formula(paste0("eval(parse(text=type)) ~ ", term))
   
    bopt <- list(myFormula=formGLM,
				   test="none",
				   family=binomial(link = 'logit'),
				   mustart=0.5,
				   control=glm.control(epsilon = 1e-06, trace = FALSE, maxit = 100))
  }
  
  if("GAM" %in% model){ #GAM options
    term <- NULL
    if(any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "s(LnP) + s(MinT)")
    if(any(level == "interaction")) term <- paste0("","ti(MinTScaled, LnPScaled)")
    formGAM <- formula(paste0("eval(parse(text=type)) ~ ", term))
   
    bopt <- list(myFormula=formGAM,
				   algo='GAM_mgcv',
				   k=NULL,                    
				   family=binomial(link = 'logit'),
				   control=gam.control(epsilon = 1e-06, 
									   trace = FALSE, 
									   maxit = 100,
									   keepData=TRUE))
  }
  
  if ("MaxEnt" %in% model) { # This is the Tsuruoka and not the Phillips implementation of MaxEnt!
	term <- NULL
	if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
#    if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
#    if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formMET <- formula(paste0("eval(parse(text=type)) ~ ", term))

 	# http://www.nactem.ac.uk/tsuruoka/maxent/
	#    supporting L1/L2 regularization [1,2]
	#    supporting fast parameter estimation algorithms (LBFGS [3], OWLQN [4], and SGD [5])
	#    supporting real-valued features
    #[1] Jun'ichi Kazama and Jun'ichi Tsujii. 2003. Evaluation and Extension of Maximum Entropy Models with Inequality Constraints, In Proceedings of EMNLP, pp. 137-144.
	#[2] Stanley F. Chen and Ronald Rosenfeld. 1999. A Gaussian Prior for Smoothing Maximum Entropy Models, Technical Report CMU-CS-99-108, Computer Science Department, Carnegie Mellon University.
	#[3] Jorge Nocedal. 1980. Updating Quasi-Newton Matrices with Limited Storage, Mathematics of Computation, Vol. 35, No. 151, pp. 773-782.
	#[4] Galen Andrew and Jianfeng Gao. 2007. Scalable training of L1-regularized log-linear models, In Proceedings of ICML.
	#[5] Yoshimasa Tsuruoka, Jun'ichi Tsujii, and Sophia Ananiadou. 2009. Stochastic Gradient Descent Training for L1-regularized Log-linear Models with Cumulative Penalty, In Proceedings of ACL-IJCNLP, pp. 477-485
	
	#Yoshimasa Tsuruoka recommends using one of following three methods if you see overfitting.
	#1. Set the l1_regularizer parameter to 1.0, leaving l2_regularizer and set_heldout as default.
	#2. Set the l2_regularizer parameter to 1.0, leaving l1_regularizer and set_heldout as default.
	#3. Set the set_heldout parameter to hold-out a portion of your data, leaving l1_regularizer and l2_regularizer as default.
	#If you are using a large number of training samples, try setting the use_sgd parameter to TRUE.

	bopt <- list(myFormula = formMET,
				l1_regularizer = 1,
				l2_regularizer = 0,
				use_sgd = TRUE, # SGD = stochastic gradient descent
				set_heldout = 0)
  }

  if ("RF" %in% model) {
    term <- NULL
    if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
#    if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
#    if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formRF <- formula(paste0("eval(parse(text=type)) ~ ", term))

	# Breiman, L (2002), “Manual On Setting Up, Using, And Understanding Random Forests V3.1”, https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf
	# Cutler, D. R., T. C. Edwards, K. H. Beard, A. Cutler, K. T. Hess, J. Gibson, and J. J. Lawler. 2007. Random forests for classification in ecology. Ecology 88:2783-2792.
	bopt <- list(myFormula = formRF,
                      ntree = 501, # Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.; Cutler et al. 2007: even ntree = 50 produced quite stable results; odd number so that there will be no ties in predictions (which would be broken at random)
                      mtry = 'default', # 'default' for classification (sqrt(p) where p is number of variables in x); Cutler et al. 2007: RF is insensitive to mtry
                      nodesize = 5, #NOTE: randomForest's default for classification is 1 (which Breiman 2002 recommends); biomod2 sets it to 5. Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
                      maxnodes = NULL)
  }

  if ("BRT" %in% model) {
    term <- NULL
    if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
    if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
#    if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formBRT <- formula(paste0("eval(parse(text=type)) ~ ", term))

	bopt <- list(myFormula = formBRT,
						n.trees = 2500,
						interaction.depth = if (any(level == "interaction")) 2L else 1L,
						n.minobsinnode = 5,
						shrinkage = 0.001,
						bag.fraction = 0.5,
						train.fraction = 1,
						cv.folds = 1, #NOTE: biomod2 has this set to 3
						keep.data = FALSE,
						verbose = FALSE)
  }
  
  bopt
}

#calculate cutoff based on maximizing TSS
get.cutoff <- function(pred, obs, pred.eval, obs.eval, method){
  
  p.cut <- seq(0,1,by=.001)
  stat <- list()
  if('TSS' %in% method) stat$TSS <- cbind(p.cut,p.cut,p.cut,p.cut)*NA
  if('KAPPA' %in% method) stat$KAPPA <- cbind(p.cut,p.cut,p.cut,p.cut)*NA
  if('ROC' %in% method) stat$ROC <- cbind(p.cut,p.cut,p.cut,p.cut)*NA
  
  for(pc in 1:length(p.cut)){
    
    p.obs <- as.integer(pred >= p.cut[pc])
    
    a <- length(which(obs == 1 & p.obs == 1))
    b <- length(which(obs == 0 & p.obs == 1))
    c <- length(which(obs == 1 & p.obs == 0))
    d <- length(which(obs == 0 & p.obs == 0))
    n <- length(obs)
    
    p.obs.eval <- as.integer(pred.eval >= p.cut[pc])
    
    ae <- length(which(obs.eval == 1 & p.obs.eval == 1))
    be <- length(which(obs.eval == 0 & p.obs.eval == 1))
    ce <- length(which(obs.eval == 1 & p.obs.eval == 0))
    de <- length(which(obs.eval == 0 & p.obs.eval == 0))
    ne <- length(obs.eval)
    
    if('TSS' %in% method){
      stat$TSS[pc,1] <-  a / (a + c) + d / (b + d) - 1		#Testing Data -- sensitivity + specificity - 1
      stat$TSS[pc,2] <-  ae / (ae + ce) + de / (be + de) - 1	#Evaluating Data -- sensitivity + specificity - 1
      stat$TSS[pc,3] <-  a / (a + c)					    	#sensitivity
      stat$TSS[pc,4] <-  d / (b + d) 							#specificity
    }
    
    if('KAPPA' %in% method){
      stat$KAPPA[pc,1] <-  ((a + d) / n - ((a + b) * (a + c) + (c + d) * (b + d)) / (n^2)) /
        (1 - ((a + b) * (a + c) + (c + d) * (b + d)) / (n^2))
      stat$KAPPA[pc,2] <-  ((ae + de) / ne - ((ae + be) * (ae + ce) + (ce + de) * (be + de)) / (ne^2)) /
        (1 - ((ae + be) * (ae + ce) + (ce + de) * (be + de)) / (ne^2))
      stat$KAPPA[pc,3] <-  a / (a + c)					    	#sensitivity
      stat$KAPPA[pc,4] <-  d / (b + d) 							#specificity
    }
    
    if('ROC' %in% method){
      stat$ROC[pc,1] <-  a / (a + c)
      stat$ROC[pc,2] <-  d / (b + d)
      stat$ROC[pc,3] <-  ae / (ae + ce)					    	#sensitivity
      stat$ROC[pc,4] <-  de / (be + de) 							#specificity
    }
    
    
  }
  
  cutoff <- matrix(NA,nrow=length(method),ncol=5)
  rownames(cutoff) <- method
  colnames(cutoff) <- c('Testing.data','Evaluating.data','Cutoff','Sensitivity','Specificity')
  
  if('TSS' %in% method) cutoff['TSS',] <- c(stat$TSS[which.max(stat$TSS[,1]),1],
                                            stat$TSS[which.max(stat$TSS[,1]),2],
                                            p.cut[which.max(stat$TSS[,1])],
                                            stat$TSS[which.max(stat$TSS[,1]),3],
                                            stat$TSS[which.max(stat$TSS[,1]),4])
  if('KAPPA' %in% method) cutoff['KAPPA',] <- c(stat$KAPPA[which.max(stat$KAPPA[,1]),1],
                                                stat$KAPPA[which.max(stat$KAPPA[,1]),2],
                                                p.cut[which.max(stat$KAPPA[,1])],
                                                stat$KAPPA[which.max(stat$KAPPA[,1]),3],
                                                stat$KAPPA[which.max(stat$KAPPA[,1]),4])
  if('ROC' %in% method) cutoff['ROC',] <- c(sum(stat$ROC[order(stat$ROC[,2]),2][-1]*diff(1 - stat$ROC[order(stat$ROC[,2]),1])),
                                            sum(stat$ROC[order(stat$ROC[,2]),4][-1]*diff(1 - stat$ROC[order(stat$ROC[,2]),3])),
                                            p.cut[which.min(abs(stat$ROC[,1] - stat$ROC[,2]))],
                                            stat$ROC[which.min(abs(stat$ROC[,1] - stat$ROC[,2])),1],
                                            stat$ROC[which.min(abs(stat$ROC[,1] - stat$ROC[,2])),2])
  
  return(cutoff)
  
}

#fit model, get cutoff, and evaluation statistics
calc_sdms <- function(type, error, model, bdat, bopt,eval.methods){
  type <<- type	 #Need to pass 'type' because of model formulae
  error <<- error
  data.tmp <- cbind(bdat$resp.var,bdat$expl.var)
  colnames(data.tmp)[1] <- type
  eval.tmp <- cbind(bdat$eval.resp.var,bdat$eval.expl.var)
  colnames(eval.tmp)[1] <- type
  
  bsdms <- list()
  
  if(model == 'GLM'){
    bsdms$m <- stats::glm(formula = bopt$myFormula,
                        family = bopt$family,
                        data = data.tmp,
                        #mustart = rep(bopt$mustart,nrow(data.tmp)),
                        control = bopt$control,
                        x = FALSE, y = FALSE)
    bsdms$DEV <- deviance(bsdms)
  }
  
  if(model == 'GAM'){
    bsdms$m <- mgcv::gam(formula = bopt$myFormula,
                       family = bopt$family,
                       data = data.tmp,
                       control = bopt$control)
    bsdms$DEV <- deviance(bsdms)
  }

  
  if (model == "MaxEnt") {
  	# http://www.nactem.ac.uk/tsuruoka/maxent/
  	
  	warning("This is the Tsuruoka and not the Phillips implementation of MaxEnt!")
  	
  	bsdms$m <- maxent::maxent(feature_matrix = data.tmp[, -1],
							code_vector = as.factor(data.tmp[, 1]),
							l1_regularizer = bopt$l1_regularizer,
							l2_regularizer = bopt$l2_regularizer,
							use_sgd = bopt$use_sgd,
							set_heldout = bopt$set_heldout,
							verbose = FALSE)
  	
    bsdms$DEV <- NA
    
    pred.obs <- as.numeric(maxent::predict.maxent(bsdms$m, feature_matrix = data.tmp[, -1])[, "1"]) # fitted
    pred.eval <- as.numeric(maxent::predict.maxent(bsdms$m, feature_matrix = eval.tmp[, -1])[, "1"])
  }

  if (model == "RF") {
  	bsdms$m <- randomForest::randomForest(x = data.tmp[, -1], y = as.factor(data.tmp[, 1]), #classification random.forest: y must be a factor
# 										formula = bopt$myFormula, #For large data sets, especially those with large number of variables, calling randomForest via the formula interface is not advised: There may be too much overhead in handling the formula.
#  										data = data.tmp,
  										ntree = bopt$ntree,
  										importance = FALSE,
  										norm.votes = TRUE,
#										strata = factor(c(0, 1)), # NOTE: biomod2 sets this, but it has no influence on result
  										nodesize = bopt$nodesize,
  										maxnodes = bopt$maxnodes)
  	
  	bsdms$DEV <- NA
    pred.obs <- as.numeric(predict(bsdms$m, newdata = data.tmp[, -1], type = "prob")[, "1"]) # fitted
    pred.eval <- as.numeric(predict(bsdms$m, newdata = eval.tmp[, -1], type = "prob")[, "1"])
  }

  if (model == "BRT") {
stop("TODO(drs): the gbm call gives a weird error message about a missing 'p'???")
	bsdms$m <- gbm::gbm(formula = bopt$myFormula,
  						distribution = "bernoulli",
  						data = data.tmp,
  						weights = rep(1, nrow(data.tmp)),
            			n.trees = bopt$n.trees, 
            			cv.folds = bopt$cv.folds,
            			interaction.depth = bopt$interaction.depth, 
            			n.minobsinnode = bopt$n.minobsinnode,
            			shrinkage = bopt$shrinkage, 
            			bag.fraction = bopt$bag.fraction,
            			train.fraction = bopt$train.fraction, 
            			verbose = FALSE)
  	bsdms$DEV <- NA
  }

  
  ##evaluate models and get cutoffs
  if (model %in% c("GLM", "GAM")) {
  	pred.eval <- predict(object = bsdms$m, newdata = eval.tmp[,-1], type = 'response', se.fit = FALSE)
  	pred.obs <- fitted(bsdms$m)
  }

  
  bsdms$cutoff <- get.cutoff(pred = pred.obs, obs = data.tmp[,1],
                             pred.eval = pred.eval, obs.eval = eval.tmp[,1],
                             method = eval.methods)	
  
  
  ##return list
  return(bsdms)
}

make_projection <- function(bsdm, newData, projName){ #Project the SDMs based on one region to the others
  
  bproj <- list(projName = projName,
                pred = as.integer(1000*round(predict(object = bsdm,
                                                     newdata = newData,
                                                     type = 'response',
                                                     se.fit = FALSE),3))
  )
  
  return(bproj)
}

get.balanced.sample <- function(obs,samp){
  
  if(length(which(obs[samp]==1)) < .5*length(obs)){
    
    samp <- c(samp[obs[samp]==1],
              samp[sample(x=which(obs[samp]==0),size=length(samp[obs[samp]==1]),replace=FALSE)])
    
  }
  
  if(length(which(obs[samp]==1)) >= .5*length(obs)){
    
    samp <- c(samp[obs[samp]==0],
              samp[sample(x=which(obs[samp]==1),size=length(samp[obs[samp]==0]),replace=FALSE)])
    
  }
  
  
  return(samp)
}

#previously used biomod 2 and now involves direct modeling and evaluation
make.SDM <- function(i){
  type <- runRequests[i, "types"]
  mlevel <- runRequests[i, "mlevels"]
  model <- runRequests[i, "models"]
  error <- runRequests[i, "errors"]
  
  #we fixed: predictorsN == 10
  DSplit <- 100 * (1 - 1/(1 + sqrt(predictorsN - 1))) #Fielding, A. H., and J. F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24:38-49.
  
  #Get data
  ibase <- (climData[, "region"] == baseRegion)
  N <- sum(ibase)
  samp <- sample(x=1:N, size=trunc(N*DSplit/100), replace=FALSE)
  
  if(equalSamples)
    samp <- get.balanced.sample(obsData[[paste(type,error,sep="_")]][ibase, runRequests[i, 'realizations']],samp)
  
  bdat <- set_Data(type = type,
                   error = error,
                   obs = obsData[[paste(type,error,sep="_")]][ibase, runRequests[i, 'realizations']],
                   dat = climData[ibase, ],
                   samp = samp,
                   run = runRequests[i, 'realizations'])
  
  #Build models -- full data
  bresM <- vector(mode="list", length=6)
  names(bresM) <- c("runID", "centerMeans", "SDMs", "Proj", "samp", "ResponseCurvePreds")
  bresM$runID <- runRequests[i, ]
  bresM$centerMeans <- centerMeansData
  bresM$samp <- samp
  
  bopt <- set_options(model = model, level=  mlevels[[mlevel]])
  
  bresM$SDMs <- calc_sdms(type = type, error = error, model = model, bdat = bdat, bopt = bopt, eval.methods = eval.methods) 
  
  #Project model onto all regions and for response curve plots
  bresM$Proj <- bresM$ResponseCurvePreds <- vector(mode="list", length=length(regions))
  orig.variables <- names(climData[,-(1:3)])
  if(length(temp <- grep("Scaled", orig.variables)) > 0){
    scaled.variables <- orig.variables[temp]
    orig.variables <- orig.variables[-temp]
  } else {
    scaled.variables <- NULL
  }
  for(ir in seq_along(regions)){
    iregion <- (climData[, "region"] == regions[ir])
    newData <- climData[iregion, c("LnP", "LnPScaled", "MinT", "MinTScaled")]
    #Project model onto region
    bresM$Proj[[ir]]$Proj <- make_projection(bsdm=bresM$SDMs$m,
                                             newData=newData,
                                             projName=paste0(model, "_region", regions[ir]))		
    
    #Prepare response curve plot predictions
    bresM$ResponseCurvePreds[[ir]] <- try(our.response.plot2(
      modelObj=bresM$SDMs$m, modelName=model,
      Data=newData,
      orig.variables=orig.variables,
      scaled.variables=scaled.variables,
      centerMeans=bresM$centerMeans,
      data_species=obsData[[paste(type,error,sep="_")]][iregion, runRequests[i, 'realizations']],
      fixed.var.metric="mean") #mean was proposed in the 'evaluation strip' by Elith et al. 2005
      , silent=TRUE)
  }
  
  
  
  #clean fitted model objects, i.e., it will NOT work with predict.glm resp. predict.gam
warning("TODO(drs): fix this; it is now SDMs$m")
  bresM$SDMs <- bresM$SDMs[c('coefficients','DEV','family','cutoff','edf')] #'model'
  
  
  return(bresM)
}


## Functions to evaluate models
rmse <- function(obs, pred, na.rm=FALSE) sqrt(mean((obs - pred) ^ 2, na.rm=na.rm))
mae <- function(obs, pred, na.rm=FALSE) mean(abs(obs - pred), na.rm=na.rm)


eval2.SDMs <- function(i, runEval, type, error, mlevel, model, runID, bsub){
  if(length(bsub) == 0) return(NULL)
  ids <- matrix(unlist(t(sapply(bsub, FUN=function(l) l$runID))[, c("realizations", "run")]), ncol=2, byrow=FALSE, dimnames=list(NULL, c("realizations", "run")))
  
  quantile.probs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
  
  #Prepare result container
  bevalM <- vector(mode="list", length=4)
  names(bevalM) <- c("evalID", "Eval", "Deviance", "Proj") #removed variable importance for now
  bevalM$evalID <- runEval
  
  #Evaluate models based on evaluation datasplits and realizations
  stat.methods <- c('Testing.data','Evaluating.data','Cutoff','Sensitivity','Specificity')
  temp.Eval <- array(NA,dim=c(length(eval.methods), length(stat.methods), temp <- apply(ids, 2, max)),
                     dimnames=list(eval.methods, stat.methods, eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))
  temp.Dev <- array(NA, dim=temp)
  
  for(j in 1:nrow(ids)){
    temp.Eval[,, ids[j, "realizations"], ids[j, "run"]] <- bsub[[j]]$SDMs$cutoff[,stat.methods]
    temp.Dev[ids[j, "realizations"], ids[j, "run"]] <- bsub[[j]]$SDMs$DEV
  }
  
  bevalM$Eval$mean <- apply(temp.Eval, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
  bevalM$Eval$sd   <- apply(temp.Eval, MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
  bevalM$Eval$quantiles   <- apply(temp.Eval, MARGIN=c(1,2),FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
  bevalM$Deviance <- c(mean=mean(temp.Dev, na.rm=TRUE), sd=sd(temp.Dev, na.rm=TRUE))
  bevalM$Deviance$quantiles   <- quantile(temp.Dev, probs=quantile.probs, type=8, na.rm=TRUE)
  
  #Project models onto all regions and take differences between projected and base region
  bevalM$Proj <- vector(mode="list", length=length(regions))
  stat.methods <- c("Testing.data", "Cutoff", "Sensitivity", "Specificity")
  
  for(ir in c(baseRegion, seq_along(regions)[-baseRegion])){
    temp.Eval <- array(NA,dim=c(length(eval.methods), length(stat.methods), temp <- apply(ids, 2, max)),
                       dimnames=list(eval.methods, stat.methods, eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))
    temp.Prob <- array(NA,dim=c(2, temp),
                       dimnames=list(c('rmse','mae'), eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))
    
    #Evaluate projections based on complete dataset
    data.probs <- probData[[paste(type,error,sep="_")]][(climData[, "region"] == regions[ir])]
    
    for(j in 1:nrow(ids)){
      data.obs <- obsData[[paste(type,error,sep="_")]][(climData[, "region"] == regions[ir]), ids[j, "realizations"]]
      FittedData <- bsub[[j]]$Proj[[ir]]$Proj$pred/1000
      temp.Eval[,, ids[j, "realizations"], ids[j, "run"]] <- 
        get.cutoff(pred = FittedData,
                   obs = data.obs,
                   pred.eval = FittedData,
                   obs.eval = data.obs,
                   method = eval.methods)[,stat.methods]
      
      temp.Prob['rmse', ids[j, "realizations"], ids[j, "run"]] <- rmse(obs=data.probs, pred = FittedData)
      temp.Prob['mae', ids[j, "realizations"], ids[j, "run"]]  <- mae(obs=data.probs, pred = FittedData)
      
    }
    
    #Differences to base region
    if(ir == baseRegion){#Code assumes that ir takes as first value the value of base region
      base.Eval <- temp.Eval
      base.Prob <- temp.Prob
    }
    diff.Eval <- temp.Eval - base.Eval
    diff.Prob <- temp.Prob - base.Prob
    
    #Aggregated evaluations
    bevalM$Proj[[ir]]$Eval$mean <- apply(temp.Eval, MARGIN=c(1, 2), FUN=mean, na.rm=TRUE)
    bevalM$Proj[[ir]]$Eval$sd   <- apply(temp.Eval, MARGIN=c(1, 2), FUN=sd, na.rm=TRUE)
    bevalM$Proj[[ir]]$Eval$quantiles <- apply(temp.Eval, MARGIN=c(1, 2), FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalDiffToBase$mean <- apply(diff.Eval, MARGIN=c(1, 2), FUN=mean, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalDiffToBase$sd   <- apply(diff.Eval, MARGIN=c(1, 2), FUN=sd, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalDiffToBase$quantiles <- apply(diff.Eval, MARGIN=c(1, 2), FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
    
    #Evaluate projections against underlying probabilities
    bevalM$Proj[[ir]]$EvalProb$mean <- apply(temp.Prob, MARGIN=1, FUN=mean, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalProb$sd <- apply(temp.Prob, MARGIN=1, FUN=sd, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalProb$quantiles <- apply(temp.Prob, MARGIN=1, FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalProbDiffToBase$mean <- apply(diff.Prob, MARGIN=1, FUN=mean, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalProbDiffToBase$sd <- apply(diff.Prob, MARGIN=1, FUN=sd, na.rm=TRUE)
    bevalM$Proj[[ir]]$EvalProbDiffToBase$quantiles <- apply(diff.Prob, MARGIN=1, FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
  }
  
  return(bevalM)
}


## Functions to plot figures and maps
map_distributions <- function(Obs, Fit, XY, model, fun = mean, maxPred = 1000, figname,save=TRUE){#TODO: make raster maps
  #Obs: dims=sites, NbRunEval*realizations
  #Fit: dims=sites, NbRunEval*realizations
  zmax <- max(Fit, maxPred)/1000
  
  preds <- apply(Fit/maxPred, MARGIN=1, FUN=function(x) fun(x))
  obs <- apply(Obs, MARGIN=1, FUN=function(x) mean(x))
  cols <- rev(terrain.colors(n=maxPred))
  xseq <- seq(min(preds),max(preds),length=maxPred)
  
  if(save){
    jpeg(width=12/2.54, height=10/2.54, file=figname,units='in',res=600)
    layout(matrix(c(1,2), nrow=2, ncol=1, byrow=TRUE), respect=TRUE, width=1, height=c(.8,.2))
    op <- par(mgp=c(1, 0.5, 0), mar=c(1,2,1,0))
  }
  
  plot(XY, pch=15, cex=0.43, col=cols[findInterval(preds,xseq)], asp=1, xlab="", ylab="", axes=FALSE)
  axis(side=1)
  axis(side=2)
  points(XY, pch=16, cex=ifelse(obs > 0, 0.06 + 0.25*obs, 0), col=col2alpha("black", 0.4))
  mtext(text=model)
  
  #Legend
  if(save) par(mar=c(2,0,0.1,0))
  plot(x=seq(from=0, to=zmax, length=maxPred), y=rep(0, maxPred), xlim=c(0, zmax), col=cols, xlab="Suitability score", axes=FALSE)
  axis(side=1, pos=-0.1, cex=0.6)
  
  if(save){
    par(op)
    dev.off()
  }
}

plot_scatterPredvsTrueProbs <- function(Prob, Fit, model, maxPred = 1000, figname){
  Fit <- Fit / maxPred
  binsize <- 0.03
  
  pdf(width=5, height=5, file=figname)
  op <- par(mgp=c(1.5, 0.5, 0), mar=c(3,3,2,1))
  if(ncol(Fit) > 10){
    xorder <- order(Prob)
    Prob <- Prob[xorder]
    Fit <- Fit[xorder, ]
    xt <- seq(from=0, to=1, by=binsize)
    binorder <- findInterval(Prob, xt)
    xProb <- aggregate(Prob, by=list(binorder), FUN=mean)$x
    quantile.probs <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    quantileFit <- apply(Fit, MARGIN=1, FUN=quantile, probs=quantile.probs, type=8)
    xquantileFit <- sapply(1:nrow(quantileFit), FUN=function(i) aggregate(quantileFit[i, ], by=list(binorder), FUN=mean)$x)
    plot(x=xProb, y=xquantileFit[, 4], type="l", xlim=c(0, 1), ylim=c(0, 1), xlab="'True' probability", ylab="Predicted probability", axes=FALSE)
    polygon(x=c(xProb, rev(xProb)), y=c(xquantileFit[, 1], rev(xquantileFit[, 7])), border=NA, col="blue")
    polygon(x=c(xProb, rev(xProb)), y=c(xquantileFit[, 2], rev(xquantileFit[, 6])), border=NA, col="red")
    polygon(x=c(xProb, rev(xProb)), y=c(xquantileFit[, 3], rev(xquantileFit[, 5])), border=NA, col="orange")
    lines(x=xProb, y=xquantileFit[, 4], lwd=2, col="black")
    legend(x="bottomright", inset=0.05, legend=c("median pred", "50% of pred", "90% of pred", "95% of pred"), bty="n", fill=c("black", "orange", "red", "blue"), border=NA)
  } else {
    matplot(x=Prob, y=Fit, xlim=c(0, 1), ylim=c(0, 1), pch=".", axes=FALSE)
  }
  axis(side=1, pos=0)
  axis(side=2, pos=0)
  lines(x=c(0, 1), y=c(0, 1), col="black", lwd=2)
  title(main=model)
  par(op)
  dev.off()
}

plot_responseCurves2 <- function(newdata, Prob, Obs, Fit, respCurvePreds, model, maxPred = 1000, centerMeans, figname){
  #Obs: vector of length sites
  #Fit: dims=sites, NbRunEval
  repN <- ncol(Obs)
  
  Fit <- Fit / maxPred
  
  orig.variables <- names(climData[,-(1:3)])
  scaled.variables <- NULL
  isScaled <- FALSE
  if(length(temp <- grep("Scaled", orig.variables)) > 0){
    isScaled <- TRUE
    scaled.variables <- orig.variables[temp]
    orig.variables <- orig.variables[-temp]
    #DataRe <- cbind(newdata$LnPScaled + centerMeans["LnPmean"], newdata$MinTScaled + centerMeans["MinTmean"])
  }
  DataRe <- newdata[, orig.variables]
  stopifnot(length(orig.variables) == 2)
  
  goodRuns <- sapply(respCurvePreds, FUN=function(l) !inherits(l, "try-error"))
  if(all(goodRuns)){
    responseAgg_TF <- (repN > 10)
    if(responseAgg_TF){
      temp <- array(unlist(respCurvePreds), dim=c(dim(respCurvePreds[[1]][[1]]), length(respCurvePreds[[1]]), repN), dimnames=list(NULL, c("Var", "Pred"), c("Var1", "Var2"), idruns=1:repN))						
      aggR2s <- array(NA, dim=c(dim(temp)[-4], 5))
      aggR2s[, , , 1] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.05, na.rm=TRUE)
      aggR2s[, , , 2] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.33, na.rm=TRUE)
      aggR2s[, , , 3] <- apply(temp, MARGIN=1:3, FUN=median, na.rm=TRUE)
      aggR2s[, , , 4] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.66, na.rm=TRUE)
      aggR2s[, , , 5] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.95, na.rm=TRUE)
    }
    
    pdf(width=12/2.54, height=12/2.54, file=figname)
    layout(matrix(c(2,0,1,3), nrow=2, ncol=2, byrow=TRUE), width=c(lcm(9), lcm(3)), height=c(lcm(3), lcm(9)))
    op <- par(mgp=c(1.5, 0.5, 0), mar=c(3,3,0,0))
    par(mar=c(3,3,0,0))
    plot(x=DataRe[, 1], y=DataRe[, 2], pch=20, cex=0.5, col="black", xlab=orig.variables[1], ylab=orig.variables[2])
    mtext(text=model, line=-1.5, adj=1)
    
    par(mar = c(0,3,0,0))
    plot(x=DataRe[, 1], y=Prob, pch=20, cex=0.25, col="darkblue", ylim=c(0, 1), xlab="", ylab="", axes=FALSE)
    axis(side=2, at=seq(0, 1, length=6), labels=FALSE)
    if(!responseAgg_TF){
      for(m in 1:ncol(Obs)){
        points(x=DataRe[, 1], y=Obs[, m], pch=20, cex=0.25, col="gray")
        if(!inherits(respCurvePreds[[m]], "try-error")){
          matlines(x=respCurvePreds[[m]][[1]][, 1] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[m]][[1]])[1]), names(centerMeans))] else 0, y=respCurvePreds[[m]][[1]][, -1], col=rainbow(n=round(1.2*repN))[m], lty=1)
        }
      }
    } else {
      x <- aggR2s[, 1, 1, 3] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[1]][[1]])[1]), names(centerMeans))] else 0
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, 1, 1], rev(aggR2s[, 2, 1, 5])), 
              border=NA, col=col2alpha("blue", 0.3))
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, 1, 2], rev(aggR2s[, 2, 1, 4])), 
              border=NA, col=col2alpha("orange", 0.9))
      lines(	x=x,
             y=aggR2s[, 2, 1, 3], lwd=1, col="black")
    }
    
    par(mar = c(3,0,0,0))
    plot(x=Prob, y=DataRe[, 2], pch=20, cex=0.25, col="darkblue", xlim=c(0, 1), xlab="", ylab="", axes=FALSE)
    axis(side=3, at=seq(0, 1, length=6), labels=FALSE)
    if(!responseAgg_TF){
      for(m in 1:ncol(Obs)){
        points(x=Obs[, m], y=DataRe[, 2], pch=20, cex=0.25, col="gray")
        if(!inherits(respCurvePreds[[m]], "try-error")){
          matlines(x=respCurvePreds[[m]][[2]][, -1], y=respCurvePreds[[m]][[2]][, 1] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[m]][[2]])[1]), names(centerMeans))] else 0, col=rainbow(n=round(1.2*repN))[m], lty=1)
        }
      }
    } else {
      x <- aggR2s[, 1, 2, 3] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[1]][[2]])[1]), names(centerMeans))] else 0
      polygon(x=c(aggR2s[, 2, 2, 1], rev(aggR2s[, 2, 2, 5])), 
              y=c(x, rev(x)), 
              border=NA, col=col2alpha("blue", 0.3))
      polygon(x=c(aggR2s[, 2, 2, 2], rev(aggR2s[, 2, 2, 4])), 
              y=c(x, rev(x)), 
              border=NA, col=col2alpha("orange", 0.9))
      lines(	x=aggR2s[, 2, 2, 3],
             y=x, lwd=1, col="black")
    }
    par(op)
    dev.off()
    res <- 0
  } else {
    res <- 1
  }
  return(res)
}

plot_responseCurves_CurvesOnly <- function(newdata, Prob, Obs, Fit, respCurvePreds, model, maxPred = 1000, centerMeans,env){
  #Obs: vector of length sites
  #Fit: dims=sites, NbRunEval
  repN <- ncol(Obs)
  
  Fit <- Fit / maxPred
  
  orig.variables <- names(climData[,-(1:3)])
  scaled.variables <- NULL
  isScaled <- FALSE
  if(length(temp <- grep("Scaled", orig.variables)) > 0){
    isScaled <- TRUE
    scaled.variables <- orig.variables[temp]
    orig.variables <- orig.variables[-temp]
    #DataRe <- cbind(newdata$LnPScaled + centerMeans["LnPmean"], newdata$MinTScaled + centerMeans["MinTmean"])
  }
  DataRe <- newdata[, orig.variables]
  stopifnot(length(orig.variables) == 2)
  
  goodRuns <- sapply(respCurvePreds, FUN=function(l) !inherits(l, "try-error"))
  if(all(goodRuns)){
    responseAgg_TF <- (repN > 10)
    if(responseAgg_TF){
      temp <- array(unlist(respCurvePreds), dim=c(dim(respCurvePreds[[1]][[1]]), length(respCurvePreds[[1]]), repN), dimnames=list(NULL, c("Var", "Pred"), c("Var1", "Var2"), idruns=1:repN))  					
      aggR2s <- array(NA, dim=c(dim(temp)[-4], 5))
      aggR2s[, , , 1] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.05, na.rm=TRUE)
      aggR2s[, , , 2] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.33, na.rm=TRUE)
      aggR2s[, , , 3] <- apply(temp, MARGIN=1:3, FUN=median, na.rm=TRUE)
      aggR2s[, , , 4] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.66, na.rm=TRUE)
      aggR2s[, , , 5] <- apply(temp, MARGIN=1:3, FUN=quantile, probs=0.95, na.rm=TRUE)
    }
    
    env.ind <- which(orig.variables == env)
    
    plot(x=DataRe[, env.ind], y=Prob, pch=20, cex=0.25, col="darkred", ylim=c(0, 1), xlab="", ylab="", axes=FALSE,frame.plot=TRUE)
    axis(side=2, at=seq(0, 1, length=6), labels=FALSE,tck=.05)
    
    xrange <- round(c(min(DataRe[,env.ind])+.1*diff(range(DataRe[,env.ind])),    #minimum + 10% of range
                      max(DataRe[,env.ind])-.1*diff(range(DataRe[,env.ind]))),2) #maximum - 10% of range
    
    print(seq(xrange[1],xrange[2],length=3))
    
    axis(side=1, at=seq(xrange[1],xrange[2],length=3), labels=FALSE,tck=.05)
    
    if(!responseAgg_TF){
      for(m in 1:ncol(Obs)){
        points(x=DataRe[, env.ind], y=Obs[, m], pch=20, cex=0.25, col="darkred")
        if(!inherits(respCurvePreds[[m]], "try-error")){
          matlines(x=respCurvePreds[[m]][[1]][, 1] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[m]][[1]])[1]), names(centerMeans))] else 0, y=respCurvePreds[[m]][[1]][, -1], col=rainbow(n=round(1.2*repN))[m], lty=1)
        }
      }
    } else {
      x <- aggR2s[, 1, env.ind, 3] + if(isScaled) centerMeans[grep(sub("Scaled", "", colnames(respCurvePreds[[1]][[env.ind]])[1]), names(centerMeans))] else 0
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, env.ind, 1], rev(aggR2s[, 2, env.ind, 5])), 
              border=NA, col=col2alpha("blue", 0.3))
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, env.ind, 2], rev(aggR2s[, 2, env.ind, 4])), 
              border=NA, col=col2alpha("orange", 0.9))
      lines(	x=x,
             y=aggR2s[, 2, env.ind, 3], lwd=1, col="black")
    }
    
    res <- 0
  } else {
    res <- 1
  }
  return(res)
}

make2.figures <- function(i, type, error, mlevel, model, runID, bsub){
  if(length(bsub) == 0) return(NULL)
  ids <- matrix(unlist(t(sapply(bsub, FUN=function(l) l$runID))[, c("realizations", "run")]), ncol=2, byrow=FALSE, dimnames=list(NULL, c("realizations", "run")))
  
  success <- list()
  for(ir in seq_along(regions)){
    coordsXY <- climData[(climData[, "region"] == regions[ir]), c("x", "y")]
    
    FittedData <- ObsData <- matrix(NA, nrow=length(bsub[[1]]$Proj[[ir]]$Proj$pred), ncol=nrow(ids))
    respCurvePreds <- vector(mode="list", length=nrow(ids))
    for(rr in 1:nrow(ids)){
      FittedData[, rr] <- bsub[[rr]]$Proj[[ir]]$Proj$pred
      ObsData[, rr] <- obsData[[paste(type,error,sep="_")]][(climData[, "region"] == regions[ir]), ids[rr, 'realizations']]
      respCurvePreds[[rr]] <- bsub[[rr]]$ResponseCurvePreds[[ir]]
    }
    
    #Figures
    #Map: 'true' vs. predicted distribution #TODO: maps of binary and/or continous predictions?
    success[[1]] <- try(map_distributions(Obs=ObsData, Fit=FittedData, XY=coordsXY, model=model, fun = mean, maxPred = 1000, figname=file.path(dir.maps, paste0("Map_TrueVsPredicted_", type,"_",error, "_", mlevel, "_", model, paste0("_region", regions[ir]), ".pdf"))), silent=TRUE)
    success[[2]] <- try(map_distributions(Obs=ObsData, Fit=FittedData, XY=coordsXY, model=model, fun = sd, maxPred = 500, figname=file.path(dir.maps, paste0("Map_UncertaintyVsPredicted_", type,"_",error, "_", mlevel, "_", model, paste0("_region", regions[ir]), ".pdf"))), silent=TRUE)
    
    #Predicted vs true probabilities
    success[[3]] <- try(plot_scatterPredvsTrueProbs(Prob=probData[[paste(type,error,sep="_")]][(climData[, "region"] == regions[ir])],
                                                    Fit=FittedData,
                                                    model=model, maxPred=1000,
                                                    figname=file.path(dir.figs, paste0("Fig_TrueVsPredictedProbs_", type,"_",error, "_", mlevel, "_", model, paste0("_region", regions[ir]), ".pdf"))), silent=TRUE)
    
    #Response curves: 'true' vs predicted -- not yet modified for new data structure
    success[[4]] <- try(plot_responseCurves2(
      newdata=climData[(climData[, "region"] == regions[ir]), c("LnP", "LnPScaled", "MinT", "MinTScaled")],
      Prob=probData[[paste(type,error,sep="_")]][(climData[, "region"] == regions[ir])],
      Obs=ObsData, Fit=FittedData,
      respCurvePreds=respCurvePreds,
      model=model, maxPred=1000,
      centerMeans=centerMeansData,
      figname=file.path(dir.figs, paste0("Fig_ResponseCurves2_", type,"_",error, "_", mlevel, "_", model, paste0("_region", regions[ir]), ".pdf"))), silent=TRUE)
    
    if(any(temp <- sapply(success, FUN=function(x) inherits(x, "try-error")))){
      success <- success[temp][[1]]
      break
    }
  }
  return(success)
}

work <- function() {
  # Note the use of the tag for sent messages: 1=ready_for_task, 2=done_task, 3=exiting
  # Note the use of the tag for received messages: 1=task, 2=done_tasks
  junk <- 0
  done <- 0
  while (done != 1) {
    # Signal being ready to receive a new task
    mpi.send.Robj(junk, 0, 1)
    
    # Receive a task
    dataForRun <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
    task_info <- mpi.get.sourcetag()
    tag <- task_info[2]
    
    if (tag == 1) { #1=ready_for_task
      print(paste("Worker processes", dataForRun$fun, dataForRun$i))
      if(dataForRun$fun == "eval"){
        result <- try(eval2.SDMs(i=dataForRun$i, runEval=dataForRun$runEval, type=dataForRun$type, error=dataForRun$error,mlevel=dataForRun$mlevel, model=dataForRun$model, runID=dataForRun$runID, bsub=dataForRun$bsub), silent=TRUE)
      } else if(dataForRun$fun == "figures"){
        result <- try(make2.figures(i=dataForRun$i, type=dataForRun$type, error=dataForRun$error, mlevel=dataForRun$mlevel, model=dataForRun$model, runID=dataForRun$runID, bsub=dataForRun$bsub), silent=TRUE)
        result <- as.integer(inherits(result, "try-error"))
      } else {
        result <- 1
      }
      mpi.send.Robj(list(i=dataForRun$i, r=result), 0, 2) # Send result back to the master
    } else if (tag == 2) { #2=done_task
      done <- 1
    }
    # We'll just ignore any unknown messages
  }
  mpi.send.Robj(junk, 0, 3)
}

