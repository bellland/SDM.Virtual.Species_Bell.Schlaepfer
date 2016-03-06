###############################################################################
#
# Bell, D. M., and D. R. Schlaepfer. Impacts of the data-generating processes and species distribution model complexity on ecological fidelity and global change predictions.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
###############################################################################



##logit and inverse logit

logit <- function(x) log(x / (1 - x))
inv.logit <- function(x) (1 + exp(-x))^-1

#inverse of rho matrix
inv.rmat <- compiler::cmpfun(function(n,rho=0.9){
  rinv <- diag((1 + rho^2),n,n)
  rinv[1,1] <- 1
  rinv[n,n] <- 1
  rinv[row(rinv) == (col(rinv)-1)] <- -rho
  rinv[row(rinv) == (col(rinv)+1)] <- -rho
  return(rinv)
})  	

#calculate R from covariance matrix
rmat  <- compiler::cmpfun(function(n, rho = 0.9) {
  mat <- diag(rep(1,n))
  mat <- rho^abs(row(mat)-col(mat))
  ((1 - rho^2)^-1)*mat
})


get_temp_fname <- function(x, base) {
	file.path(dir.sdm, paste(x[c("models", "types")], collapse = "_"), paste0(base, ".rds"))
} 



## Functions to read specimen data
calc_ObservationsFromProbabilities <- compiler::cmpfun(function(probs, N, VAR, sigma,w){###need to incorporate long and Lat
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
})


get_TypeData_FromFile <- compiler::cmpfun(function(type, center=TRUE, centerBasedOnRegionIDs=2){
  dat <- read.csv(file=file.path(dir.in, paste0(type, ".csv")))
  dat <- dat[, c("region", "x", "y", "LnP", "MinT", "prob")]
  
  if(center){###DMB revision: scaled variables need to be centered with the same mean for all regions
    mean_dat <- apply(dat[dat$region %in% centerBasedOnRegionIDs, c("LnP", "MinT")], 2, mean) 
    centerMeans <- c(LnPmean=mean(mean_dat["LnP"]), MinTmean=mean(mean_dat["MinT"]))
    dat$LnPScaled <- dat$LnP - centerMeans[1]
    dat$MinTScaled <- dat$MinT - centerMeans[2]
  }
  
  return(list(dat=dat, centerMeans=centerMeans))
})

get_GeographicRaster_FromFile <- function(region){
  raster(file.path(dir.gis, paste0("extent.raster.r", region, ".grd")))
}

set_Data <- compiler::cmpfun(function(type, error, obs, dat, samp, run){
  
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
})


## Functions to build the SDMs
set_options <- compiler::cmpfun(function(model, level=c("linear", "squared", "interaction")){
  level <- match.arg(arg=level, choices=c("linear", "squared", "interaction"), several.ok=TRUE)
  
  bopt <- list()
  
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
				   k=NULL,                    
				   family=binomial(link = 'logit'),
				   control=mgcv::gam.control(epsilon = 1e-06, 
									   trace = FALSE, 
									   maxit = 100,
									   keepData=TRUE))
  }
  
  if ("MaxEnt" %in% model) { # This is the Tsuruoka and not the Phillips implementation of MaxEnt!
	term <- NULL
	if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
    if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
    if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
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
	#1. Set the l1_regularizer (lasso) parameter to 1.0, leaving l2_regularizer (ridge) and set_heldout as default.
	#2. Set the l2_regularizer parameter to 1.0, leaving l1_regularizer and set_heldout as default.
	#3. Set the set_heldout parameter to hold-out a portion of your data, leaving l1_regularizer and l2_regularizer as default.
	#If you are using a large number of training samples, try setting the use_sgd parameter to TRUE.

	bopt <- list(myFormula = formMET,
				l1_regularizer = 1,
				l2_regularizer = 0,
				use_sgd = TRUE, # SGD = stochastic gradient descent
				set_heldout = 0)
  }

  if ("MaxEntP" %in% model) { # This is the Phillips implementation of MaxEnt!
  	# Requires their java application: download from http://www.cs.princeton.edu/~schapire/maxent
	
	bopt <- list(path_to_maxent.jar = path_to_MaxEntP,
				memory_allocated = 512,
				maximumiterations = 500, # default = 500; biomod2 = 200 (but not used?)
				linear = TRUE,
				quadratic = TRUE,
				product = any(level == "interaction"), #interactions
				threshold = TRUE,
				hinge = TRUE,
				lq2lqptthreshold = 80,
				l2lqthreshold = 10,
				hingethreshold = 15,
				beta_threshold = -1.0,
				beta_categorical = -1.0,
				beta_lqp = -1.0,
				beta_hinge = -1.0,
				betamultiplier = 1, # regularization? 0, no regularization
				defaultprevalence = 0.5)
	
  }

  if ("RF" %in% model) {
    term <- NULL
    if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
    #if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
    #if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formRF <- formula(paste0("eval(parse(text=type)) ~ ", term))

	# Breiman, L (2002), “Manual On Setting Up, Using, And Understanding Random Forests V3.1”, https://www.stat.berkeley.edu/~breiman/Using_random_forests_V3.1.pdf
	# Cutler, D. R., T. C. Edwards, K. H. Beard, A. Cutler, K. T. Hess, J. Gibson, and J. J. Lawler. 2007. Random forests for classification in ecology. Ecology 88:2783-2792.
	bopt <- list(myFormula = formRF,
                      ntree = 501, # Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.; Cutler et al. 2007: even ntree = 50 produced quite stable results; odd number so that there will be no ties in predictions (which would be broken at random)
                      mtry = if (any(level == "interaction")) 2L else 1L, # 'default' for classification (sqrt(p) where p is number of variables in x); Cutler et al. 2007: RF is insensitive to mtry
                      nodesize = if (any(level == "interaction")) 1L else 5L, #NOTE: randomForest's default for classification is 1 (which Breiman 2002 recommends); biomod2 sets it to 5. Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
                      maxnodes = NULL)
  }

  if ("BRT" %in% model) {
    term <- NULL
    if (any(level == "linear")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "LnP + MinT")
    if (any(level == "squared")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(LnPScaled * LnPScaled) + I(MinTScaled * MinTScaled)")
   #if (any(level == "interaction")) term <- paste0(term, ifelse(length(term) > 0, " + ", ""), "I(MinTScaled * LnPScaled)")
    formBRT <- formula(paste0("eval(parse(text=type)) ~ ", term))
	
	# gbm vignette: "I usually aim for 3,000 to 10,000 iterations with shrinkage rates between 0.01 and 0.001"
	# 5/BRT/woInt/AIF/binom/1/1: 	shrinkage	n.trees.optim
	#								0.001		17534
	#								0.005		3425
	#								0.01		2588
	bopt <- list(myFormula = formBRT,
						n.trees = 3500, # biomod2 has 2500, but also shrinkage 0.001 (and doesn't use gbm.more), they likely use too few trees; according to vignette: it may need 10,000 trees to get optimal performance with 0.001 shrinkage
						interaction.depth = if (any(level == "interaction")) 4L else 1L,
						n.minobsinnode = 5,
						shrinkage = 0.005, #biomod2 has 0.001
						bag.fraction = 0.5, # gbm vignette: "(0.5 is recommended), gbm computes an out-of-bag estimate of the improvement in predictive performance"
						train.fraction = 1,
						perf.method = "cv",
						cv.folds = 5, #NOTE: biomod2 has this set to 3; 1 expresses a bug in gbm::gbm; gbm vignette: "My recommendation is to use 5- or 10-fold cross validation if you can afford the computing time"
						keep.data = FALSE,
						verbose = FALSE)
	if (!(bopt$perf.method == "cv")) bopt$cv.folds <- 0
  }
  
  bopt
})

#calculate cutoff based on maximizing TSS
get.cutoff_slow <- function(pred, obs, pred.eval, obs.eval, method){
  
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

# the new get.cutoff() is 25% faster than the older version get.cutoff_slow
get.cutoff <- compiler::cmpfun(function(pred, obs, pred.eval, obs.eval, method){
  p.cut <- seq(0,1,by=.001)
  responses1 <- c('Testing.data', 'Evaluating.data', 'Sensitivity', 'Specificity')
  stat <- array(NA, dim = c(length(method), length(p.cut), length(responses1)), dimnames = list(method, NULL, responses1))
  responses2 <- c('Testing.data', 'Evaluating.data', 'Cutoff', 'Sensitivity', 'Specificity')
  cutoff <- array(NA, dim = c(length(method), length(responses2)), dimnames = list(method, responses2))

  for (pc in seq_along(p.cut)) {
    p.obs <- as.integer(pred >= p.cut[pc])
    lobs <- list(temp <- obs == 0, !temp)
    lpobs <- list(temp <- p.obs == 0, !temp)
    a <- sum(lobs[[2]] & lpobs[[2]])
    b <- sum(lobs[[1]] & lpobs[[2]])
    c <- sum(lobs[[2]] & lpobs[[1]])
    d <- sum(lobs[[1]] & lpobs[[1]])
    n <- length(obs)

    p.obs.eval <- as.integer(pred.eval >= p.cut[pc])
    lobse <- list(temp <- obs.eval == 0, !temp)
    lpobse <- list(temp <- p.obs.eval == 0, !temp)
    ae <- sum(lobse[[2]] & lpobse[[2]])
    be <- sum(lobse[[1]] & lpobse[[2]])
    ce <- sum(lobse[[2]] & lpobse[[1]])
    de <- sum(lobse[[1]] & lpobse[[1]])
    ne <- length(obs.eval)

	sensitivity <- a / (a + c)
	specificity <- d / (b + d)
    
    if('TSS' %in% method){
      stat["TSS", pc, "Testing.data"] <-  a / (a + c) + d / (b + d) - 1		#Testing Data -- sensitivity + specificity - 1
      stat["TSS", pc, "Evaluating.data"] <-  ae / (ae + ce) + de / (be + de) - 1	#Evaluating Data -- sensitivity + specificity - 1
      stat["TSS", pc, "Sensitivity"] <- sensitivity
      stat["TSS", pc, "Specificity"] <- specificity
    }
    
    if('KAPPA' %in% method){
      stat["KAPPA", pc, "Testing.data"] <-  ((a + d) / n - ((a + b) * (a + c) + (c + d) * (b + d)) / (n^2)) /
        (1 - ((a + b) * (a + c) + (c + d) * (b + d)) / (n^2))
      stat["KAPPA", pc, "Evaluating.data"] <-  ((ae + de) / ne - ((ae + be) * (ae + ce) + (ce + de) * (be + de)) / (ne^2)) /
        (1 - ((ae + be) * (ae + ce) + (ce + de) * (be + de)) / (ne^2))
      stat["KAPPA", pc, "Sensitivity"] <- sensitivity
      stat["KAPPA", pc, "Specificity"] <- specificity
    }
    
    if('ROC' %in% method){
      stat["ROC", pc, "Testing.data"] <-  sensitivity
      stat["ROC", pc, "Evaluating.data"] <-  specificity
      stat["ROC", pc, "Sensitivity"] <-  ae / (ae + ce)
      stat["ROC", pc, "Specificity"] <-  de / (be + de)
    }
  }


  if ('TSS' %in% method) {
  	imax <- which.max(stat["TSS", , 'Testing.data'])
  	cutoff['TSS', responses1] <- stat["TSS", imax, responses1]
  	cutoff['TSS', 'Cutoff'] <- p.cut[imax]
  }
  if ('KAPPA' %in% method) {
  	imax <- which.max(stat["KAPPA", , "Testing.data"])
  	cutoff['KAPPA', responses1] <- stat["KAPPA", imax, responses1]
  	cutoff['KAPPA', 'Cutoff'] <- p.cut[imax]
  }
  if ('ROC' %in% method) {
  	ROCsorted <- stat["ROC", order(stat["ROC", , "Evaluating.data"]), ]
  	cutoff['ROC', "Testing.data"] <- sum(ROCsorted[-1, "Evaluating.data"] * diff(1 - ROCsorted[, "Testing.data"]), na.rm = TRUE)
  	cutoff['ROC', "Evaluating.data"] <- sum(ROCsorted[-1, "Specificity"] * diff(1 - ROCsorted[, "Sensitivity"]), na.rm = TRUE)
	imax <- which.min(abs(stat["ROC", , "Testing.data"] - stat["ROC", , "Evaluating.data"]))
  	cutoff['ROC', 'Cutoff'] <- p.cut[imax]
	cutoff['ROC', "Sensitivity"] <- stat["ROC", imax, "Testing.data"]
	cutoff['ROC', "Specificity"] <- stat["ROC", imax, "Evaluating.data"]
  }
  
  cutoff
})
	


make_prediction <- compiler::cmpfun(function(bsdm, newData) {
	if (inherits(bsdm, "glm") || inherits(bsdm, "gam")) { #gam objects inherit from classes glm and lm
		preds <- predict(bsdm, newdata = newData, type = 'response', se.fit = FALSE)

	} else if (inherits(bsdm, "maxent")) {
		preds <- predict(bsdm, feature_matrix = newData[, c("LnP", "MinT")])[, "1"]
	
	} else if (inherits(bsdm, "MaxEntP")) {
		stop("predict.MaxEntP not implemented yet")
	
	} else if (inherits(bsdm, "randomForest")) {
		preds <- predict(bsdm, newdata = newData, type = "prob")[, "1"]
    
    } else if (inherits(bsdm, "gbm")) {
    	preds <- predict(bsdm, newdata = newData, type = "response", n.trees = bsdm[["n.trees.opt"]])

	} else {
		preds <- NULL
	}	
	
	as.numeric(preds)
})

make_integer <- compiler::cmpfun(function(x) as.integer(1000 * round(x, 3)))


calc_sdms_and_predict <- compiler::cmpfun(function(runID, bdat, bopt, new_regions, new_data, scaled.variables, orig.variables, centerMeans, eval_disc.methods, dir_temp = NULL){
	# Data for model fitting
	vars <- c(orig.variables, scaled.variables)
	type <<- runID["types"]	 #Need to pass 'type' because of model formulae
	data.tmp <- cbind(bdat$resp.var, bdat$expl.var[, vars])
	colnames(data.tmp)[1] <- runID["types"]

	# Data for model evaluation
	eval.tmp <- cbind(bdat$eval.resp.var, bdat$eval.expl.var[, vars])
	colnames(eval.tmp)[1] <- runID["types"]
  
	# Data for predictions
	ptemp_obs <- cbind(species = rep("obs", nrow(data.tmp)), bdat[["resp.xy"]][, c("x", "y")], data.tmp[, vars])
	ptemp_eval <- cbind(species = rep("eval", nrow(eval.tmp)), bdat[["eval.resp.xy"]][, c("x", "y")], eval.tmp[, vars])
	ptemp_regions <- do.call(rbind, lapply(seq_along(new_regions), function(ir)
						cbind(species = names(new_regions)[ir], new_data[new_regions[[ir]], c("x", "y", vars)])))
	dat_curves <- lapply(seq_along(new_regions), function(ir)
							pred.response.plot2(data_species = new_data[new_regions[[ir]], "obs"],
										data_env = ptemp_regions[ptemp_regions[, "species"] == names(new_regions)[ir], vars],
										orig.variables, scaled.variables, centerMeans,
										fixed.var.metric = 'mean')) #mean was proposed in the 'evaluation strip' by Elith et al. 2005
	ptemp_curves <- do.call(rbind, do.call(rbind, lapply(seq_along(new_regions), function(ir) lapply(seq_along(dat_curves[[ir]]), function(ic) {
						x <- dat_curves[[ir]][[ic]][["new_data"]]
						cbind(species = paste0("region", ir, "_rc_", names(dat_curves[[ir]])[ic]), x = 0, y = 0, x[, vars])}))))
	
	dat_pred <- rbind(ptemp_obs, ptemp_eval, ptemp_regions, ptemp_curves)
	id_pred <- dat_pred[, "species"]


	# Model fitting
	bsdms <- list()
	bsdms$DEV <- NA
  
	if (runID["models"] == 'GLM') {
		bsdms$comp_time <- system.time(
			temp <- stats::glm(formula = bopt$myFormula,
								family = bopt$family,
								data = data.tmp,
								#mustart = rep(bopt$mustart,nrow(data.tmp)),
								control = bopt$control,
								x = FALSE, y = FALSE))["elapsed"]
		bsdms$m <- temp
		bsdms$DEV <- deviance(bsdms$m)
	} else
  
	if (runID["models"] == 'GAM') {
		bsdms$comp_time <- system.time(
			temp <- mgcv::gam(formula = bopt$myFormula,
						   family = bopt$family,
						   data = data.tmp,
						   control = bopt$control))["elapsed"]
		bsdms$m <- temp
		bsdms$DEV <- deviance(bsdms$m)
	} else 
  
	if (runID["models"] == "MaxEnt") {
#		# http://www.nactem.ac.uk/tsuruoka/maxent/
#		warning("This is the Tsuruoka and not the Phillips implementation of MaxEnt!")
#	
#		if (any(grepl("I(MinTScaled * LnPScaled)", as.character(attr(terms(bopt$myFormula), "variables")), fixed = TRUE))) {
#			# with interaction
#			feature_matrix <- with(data.tmp, data.frame(LnP = LnP, MinT = MinT, LnP_x_MinT = MinTScaled * LnPScaled))
#		} else { #no interaction
#			feature_matrix <- with(data.tmp, data.frame(LnP = LnP, MinT = MinT))
#		}
  	
		bsdms$comp_time <- system.time(
			temp <- maxent::maxent(feature_matrix = data.tmp[, c("LnP", "MinT")],
								code_vector = as.factor(data.tmp[, 1]),
								l1_regularizer = bopt$l1_regularizer,
								l2_regularizer = bopt$l2_regularizer,
								use_sgd = bopt$use_sgd,
								set_heldout = bopt$set_heldout,
								verbose = FALSE))["elapsed"]
		bsdms$m <- temp
	} else 
  
	if (runID["models"] == "MaxEntP") {		
		# Setup Phillips-Maxent
		#	- directories
		dir_temp_MEP <- file.path(dir_temp, paste(runID, collapse = "_"))
		dir_work_MEP <- file.path(dir_temp_MEP, "MaxEnt_Phillips_work")
		dir_out_MEP <- file.path(dir_temp_MEP, "MaxEnt_Phillips_output")
		temp <- lapply(c(dir_work_MEP, dir_out_MEP), function(x) dir.create(path = x, showWarnings = FALSE, recursive = TRUE))
		
		#	- observational data
		sp_name <- c(absence = "background", presence = as.character(runID["types"]))
		iobs <- lapply(c(0, 1), function(x) data.tmp[, 1] == x)
		names(iobs) <- names(sp_name)
		dat_obs_MEP <- lapply(seq_along(sp_name), function(x)
								cbind(species = rep(sp_name[[x]], sum(iobs[[x]])), 
								bdat[["resp.xy"]][iobs[[x]], c("x", "y")],
								data.tmp[iobs[[x]], c("LnP", "MinT")]))
		names(dat_obs_MEP) <- names(sp_name)
		temp <- lapply(names(dat_obs_MEP), function(x) write.table(dat_obs_MEP[[x]],
									file = file.path(dir_work_MEP, paste0(x, ".csv")),
									quote = FALSE, row.names = FALSE, sep = ","))
		
		#	- data for predictions
		write.table(dat_pred[, c("species", "x", "y", "LnP", "MinT")], file = file.path(dir_work_MEP, "predict.csv"), quote = FALSE, row.names = FALSE, sep = ",")

  		# Run Phillips-Maxent
		comp_time2 <- system.time(
			temp <- system2(command = "java", args = paste0(
							 "-mx", bopt$memory_allocated, "m",
							 " -jar ", file.path(bopt$path_to_maxent.jar, "maxent.jar"), 
							 " environmentallayers=\"", file.path(dir_work_MEP, "absence.csv"),
							 "\" samplesfile=\"", file.path(dir_work_MEP, "presence.csv"),
							 "\" projectionlayers=\"", file.path(dir_work_MEP, "predict.csv"), 
							 "\" outputdirectory=\"", dir_out_MEP, "\"",
							 " outputformat=logistic ",
							 " maximumiterations=", bopt$maximumiterations,
							 " linear=", bopt$linear,
							 " quadratic=", bopt$quadratic,
							 " product=", bopt$product,
							 " threshold=", bopt$threshold,
							 " hinge=", bopt$hinge,
							 " lq2lqptthreshold=", bopt$lq2lqptthreshold,
							 " l2lqthreshold=", bopt$l2lqthreshold,
							 " hingethreshold=", bopt$hingethreshold,
							 " beta_threshold=", bopt$beta_threshold,
							 " beta_categorical=", bopt$beta_categorical,
							 " beta_lqp=", bopt$beta_lqp,
							 " beta_hinge=", bopt$beta_hinge,
							 " betamultiplier=", bopt$betamultiplier,
							 " defaultprevalence=", bopt$defaultprevalence,
							 " autorun redoifexists", 
							 " plots nodoclamp novisible nowarnings notooltips noaddsamplestobackground",
							 " nowritebackgroundpredictions"),
					stdout = TRUE, stderr = TRUE)
		)["elapsed"]
	
		# computational time as reported by MaxEntP
		temp <- readLines(con = file.path(dir_out_MEP, "maxent.log"), n = 38)[38]
		bsdms$comp_time <- as.numeric(strsplit(temp, split = " ", fixed = TRUE)[[1]][4])
		bsdms$comp_time_total <- comp_time2
	
		# model
		m <- list(	lambdas = readLines(con = file.path(dir_out_MEP, paste0(runID["types"], ".lambdas"))),
					evals = read.csv(file.path(dir_out_MEP, "maxentResults.csv")))
		class(m) <- "MaxEntP"
		bsdms$m <- m
	} else 

	if (runID["models"] == "RF") {
		bsdms$comp_time <- system.time(
			temp <- randomForest::randomForest(x = data.tmp[, c("LnP", "MinT")], y = as.factor(data.tmp[, 1]), #classification random.forest: y must be a factor
#	 										formula = bopt$myFormula, #For large data sets, especially those with large number of variables, calling randomForest via the formula interface is not advised: There may be too much overhead in handling the formula.
#	  										data = data.tmp,
											ntree = bopt$ntree,
											mtry = bopt$mtry,
											importance = FALSE,
											norm.votes = TRUE,
#											strata = factor(c(0, 1)), # NOTE: biomod2 sets this, but it has no influence on result
											nodesize = bopt$nodesize,
											maxnodes = bopt$maxnodes))["elapsed"]

# balancing classes seems to worsen evaluation performance
#temp2 <- rfUtilities::rf.classBalance(xdata = data.tmp[, c("LnP", "MinT")], ydata = as.factor(data.tmp[, 1]), #classification random.forest: y must be a factor
##	 										formula = bopt$myFormula, #For large data sets, especially those with large number of variables, calling randomForest via the formula interface is not advised: There may be too much overhead in handling the formula.
##	  										data = data.tmp,
#											ntree = bopt$ntree,
#											mtry = bopt$mtry,
#											importance = FALSE,
#											norm.votes = TRUE,
##											strata = factor(c(0, 1)), # NOTE: biomod2 sets this, but it has no influence on result
#											nodesize = bopt$nodesize,
#											maxnodes = bopt$maxnodes)
	
		temp[["OOB"]] <- temp$err.rate[temp$ntree, "OOB"] #OOB estimate of error rate
	
		bsdms$m <- temp
	} else 

	if (runID["models"] == "BRT") {
		its <- comp_time2 <- comp_time1 <- n.trees.opt <- 0
		temp <- list()
	
		# gbm vignette: "I usually aim for 3,000 to 10,000 iterations with shrinkage rates between 0.01 and 0.001"
		comp_time2 <- try(system.time(
			while (bopt$shrinkage <= 0.1 && bopt$shrinkage >= 0.0001 && its < 10) {
				its <- its + 1
				comp_time1 <- try(system.time(
					temp <- gbm::gbm(formula = bopt$myFormula,
									distribution = "bernoulli",
									data = data.tmp,
									n.trees = bopt$n.trees, 
									cv.folds = bopt$cv.folds,
									interaction.depth = bopt$interaction.depth, 
									n.minobsinnode = bopt$n.minobsinnode,
									shrinkage = bopt$shrinkage, 
									bag.fraction = bopt$bag.fraction,
									train.fraction = bopt$train.fraction, 
									verbose = FALSE,
									n.cores = 0) # we don't want this to run parallel; the calls to make.SDMs() are parallelized
				)["elapsed"], silent = TRUE)
				if (inherits(comp_time1, "try-error")) stop(comp_time1)
				
				n.trees.opt <- gbm::gbm.perf(temp, method = bopt$perf.method, plot.it = FALSE)

				if (n.trees.opt <= 3000) {
					bopt$shrinkage <- bopt$shrinkage / 5
				} else if (n.trees.opt >= 10000) {
					bopt$shrinkage <- bopt$shrinkage * 2
				} else {
					if (n.trees.opt >= bopt$n.trees) {
						bopt$n.trees <- 2 * bopt$n.trees
					} else {
						break # success!
					}
				}
			}
		)["elapsed"], silent = TRUE)
		
		if (inherits(comp_time2, "try-error")) stop(runID, " failed during call to 'gbm': ", comp_time2)
		if (its >= 10) stop(runID, " did not reach convergence")

		bsdms$comp_time <- comp_time1
		bsdms$comp_time_total <- comp_time2
		
		bsdms$shrinkage <- bopt$shrinkage
		bsdms$n.trees <- bopt$n.trees

		temp[["n.trees.opt"]] <- n.trees.opt
		temp[["its"]] <- its
		bsdms$m <- temp
	}

	# Make the predictions
	if (!is.null(bsdms$m)) {
		if (runID["models"] == "MaxEntP") {
			# predictions
			preds <- read.csv(file.path(dir_out_MEP, paste0(runID["types"], "_predict.csv")))[, 3]
	
			# remove tmp dir
			unlink(dir_temp_MEP, recursive = TRUE)
		} else {
			preds <- make_prediction(bsdms$m, dat_pred)
		}
		
		# Evaluate models and get cutoffs
		pred.eval <- preds[id_pred == "eval"]
		pred.obs <- preds[id_pred == "obs"]
		bsdms$cutoff <- get.cutoff(pred = pred.obs, obs = data.tmp[,1],
								 pred.eval = pred.eval, obs.eval = eval.tmp[,1],
								 method = eval_disc.methods)	
		
		# Predict across regions
		bsdms$Proj <- lapply(names(new_regions), function(ir) preds[id_pred == ir])
		names(bsdms$Proj) <- names(new_regions)
		
		# Get values for response curves
		bsdms$ResponseCurvePreds <- lapply(seq_along(new_regions), function(ir) {
										temp <- lapply(paste0(names(new_regions)[ir], "_rc_", names(dat_curves[[ir]])), function(x) preds[id_pred == x]) 
										arrange.response.plot2(show.variables = names(dat_curves[[ir]]),
																modelName = runID[["models"]],
																xdat = dat_curves[[ir]],
																ydat = temp)
									})
		names(bsdms$ResponseCurvePreds) <- names(new_regions)
	}
	

	# Return result
	bsdms
})




get.balanced.sample <- compiler::cmpfun(function(obs,samp){
  
  if(length(which(obs[samp]==1)) < .5*length(obs)){
    
    samp <- c(samp[obs[samp]==1],
              samp[sample(x=which(obs[samp]==0),size=length(samp[obs[samp]==1]),replace=FALSE)])
    
  }
  
  if(length(which(obs[samp]==1)) >= .5*length(obs)){
    
    samp <- c(samp[obs[samp]==0],
              samp[sample(x=which(obs[samp]==1),size=length(samp[obs[samp]==0]),replace=FALSE)])
    
  }
  
  
  return(samp)
})

#previously used biomod 2 and now involves direct modeling and evaluation
make.SDM <- compiler::cmpfun(function(i){
	ftemp <- get_temp_fname(runRequests[i, ], runRequestIDs[i])
	if (!file.exists(ftemp)) {
		print(paste(Sys.time(), "make.SDM:", i, paste(runRequests[i, ], collapse = "-")))

		type <- runRequests[i, "types"]
		mlevel <- runRequests[i, "mlevels"]
		model <- runRequests[i, "models"]
		error <- runRequests[i, "errors"]

		orig.variables <- names(climData[,-(1:3)])
		if(length(temp <- grep("Scaled", orig.variables)) > 0){
			scaled.variables <- orig.variables[temp]
			orig.variables <- orig.variables[-temp]
		} else {
			scaled.variables <- NULL
		}
		
		#we fixed: predictorsN == 10
		DSplit <- 100 * (1 - 1/(1 + sqrt(predictorsN - 1))) #Fielding, A. H., and J. F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24:38-49.

		#Get data
		i_regions <- lapply(regions, function(x) climData[, "region"] == x)
		names(i_regions) <- paste0("region", regions)
		obs <- obsData[[paste(type,error,sep="_")]][, runRequests[i, 'realizations']]
		
		N <- sum(i_regions[[baseRegion]])
		samp <- sample(x=1:N, size=trunc(N*DSplit/100), replace=FALSE)

		if (equalSamples) samp <- get.balanced.sample(obs[i_regions[[baseRegion]]],samp)
		
		bdat <- set_Data(type = type,
					   error = error,
					   obs = obs[i_regions[[baseRegion]]],
					   dat = climData[i_regions[[baseRegion]], ],
					   samp = samp,
					   run = runRequests[i, 'realizations'])
  
		#Build models -- full data
		bresM <- list()
		bresM[["runID"]] <- runRequests[i, ]
		bresM[["centerMeans"]] <- centerMeansData
		bresM[["samp"]] <- samp

		bopt <- set_options(model = model, level = mlevels[[mlevel]])

		temp <- try(calc_sdms_and_predict(runID = bresM[["runID"]], bdat = bdat, bopt = bopt,
					new_regions = i_regions,
					new_data = cbind(obs, climData),
					scaled.variables = scaled.variables, orig.variables = orig.variables,
					centerMeans = bresM[["centerMeans"]],
					eval_disc.methods = eval_disc.methods,
					dir_temp = dirname(ftemp)), silent = TRUE)

		if (inherits(temp, "try-error")) {
			res <- temp
		} else {
			bresM <- modifyList(bresM, temp)

			#clean fitted model objects, i.e., it will NOT work with predict.glm etc.
			bresM[["m"]] <- if (identical(model, "GLM")) {
								bresM[["m"]]$family <- bresM[["m"]]$family[c("family", "link")]
								bresM[["m"]][c('coefficients', 'family', 'df.null', 'df.residual')]
							} else if (identical(model, "GAM")) { # gam objects inherit from classes glm and lm
								bresM[["m"]]$family <- bresM[["m"]]$family[c("family", "link")]
								bresM[["m"]][c('coefficients', 'family', 'df.null', 'df.residual', 'edf')]
							} else if (identical(model, "MaxEnt")) {
								bresM[["m"]]
							} else if (identical(model, "MaxEntP")) {
								bresM[["m"]]
							} else if (identical(model, "RF")) {
								bresM[["m"]]$class <- "randomForest"
								bresM[["m"]][c("type", "importance", "ntree", "mtry", "confusion", "OOB")]
							} else if (identical(model, "BRT")) {
								bresM[["m"]][c("initF", "n.trees.opt", "distribution", "interaction.depth", "its")]
							}
			bresM[["class"]] <- switch(EXPR = model, GLM = "glm", GAM = "gam", MaxEnt = "maxent", MaxEntP = "MaxEntP", RF = "randomForest", BRT = "gbm")
				
			# Memoize the results
			saveRDS(bresM, ftemp)
			
			res <- i
		}
	}  
  
	res
})


## Functions to evaluate models
rmse <- compiler::cmpfun(function(obs, pred, na.rm=FALSE) sqrt(mean((obs - pred) ^ 2, na.rm=na.rm)))
mae <- compiler::cmpfun(function(obs, pred, na.rm=FALSE) mean(abs(obs - pred), na.rm=na.rm))


calc_eval_regions <- compiler::cmpfun(function(i) {
	stat.methods <- c("Testing.data", "Evaluating.data", "Cutoff", "Sensitivity", "Specificity")
	
	res <- list(discrete = array(NA, dim = c(length(regions), length(eval_disc.methods), length(stat.methods)), dimnames = list(regions, eval_disc.methods, stat.methods)),
				cont = array(NA, dim = c(length(regions), length(eval_cont.methods)), dimnames = list(regions, eval_cont.methods)))

  	ftemp <- get_temp_fname(runRequests[i, ], runRequestIDs[i])
  	if (!file.exists(ftemp)) {
  		return(res)
  	} else {
  		bresM <- try(readRDS(file = ftemp), silent = TRUE)
		
		if (!inherits(bresM, "try-error") && !is.null(bresM$evals)) {
			return(bresM$evals)
		} else {
			for (ir in seq_along(regions)) {
				obs  <- obsData[[paste(runRequests[i,'types'], runRequests[i,'errors'], sep="_")]][, runRequests[i,'realizations']][climData[,'region'] == ir]
				pred <- bresM$Proj[[ir]] / 1000

				res[["discrete"]][ir, eval_disc.methods, stat.methods] <- get.cutoff(pred = pred, obs = obs,
																					   pred.eval = pred, obs.eval = obs,
																					   method = eval_disc.methods)

				prob <- probData[[paste(runRequests[i,'types'], runRequests[i,'errors'], sep="_")]][climData[,'region'] == ir]
				res[["cont"]][ir, 'RMSE'] <- rmse(obs = prob, pred = pred)
				res[["cont"]][ir, 'MAE'] <- mae(obs = prob, pred = pred)
			}

			bresM$evals <- res
			saveRDS(bresM, file = ftemp)
			print(paste(Sys.time(), i, runRequestIDs[i], ": evaluation for each region saved to disk"))
		}
	}
	
	res
})


get_region_eval <- compiler::cmpfun(function(i, ir, stat.methods) {
	outs <- c(eval_disc.methods, eval_cont.methods)
	res <- matrix(NA, nrow = length(stat.methods), ncol = length(outs), dimnames = list(stat.methods, outs))

  	ftemp <- get_temp_fname(runRequests[i, ], runRequestIDs[i])
  	if (file.exists(ftemp)) {
  		bresM <- try(readRDS(file = ftemp), silent = TRUE)

		if (!inherits(bresM, "try-error")) {
			temp <- if (is.null(bresM$evals)) calc_eval_regions(i) else bresM$evals
			res[, eval_disc.methods] <- t(temp[["discrete"]][ir, eval_disc.methods, stat.methods])
			res[, eval_cont.methods] <- rep(temp[["cont"]][ir, eval_cont.methods], each = length(stat.methods))
		}
	}
		
	res
})

calc_region_partition <- compiler::cmpfun(function(j) {
	ftemp2 <- file.path(dir.tables, paste0("Partition_Props_region", ir, "_var", j, "_temp2.rds"))
	if (action == "continue" && file.exists(ftemp2)) {
		prop <- readRDS(file = ftemp2)
	} else {
		print(paste(Sys.time(), ": Partition of region:", ir, "; variable:", j, "of", length(variables)))

		tmp <- part.mat[,c(variables[j], colnames(runRequests))]
		colnames(tmp) <- c('y',colnames(runRequests))

		tmp$realizations <- apply(tmp[,factors],1,paste,collapse=".")
		
		# with 320,000 cases, this uses about 31 GB of memory
		mfit <- aov(y ~ factor(types)*factor(errors)*factor(models)*factor(mlevels) + Error(realizations), data=tmp) 
		sfit <- summary(mfit)

		if (length(sfit) == 2) {
			fnames <- rownames(sfit[["Error: realizations"]][[1]])
			fnames <- c("realizations", gsub("factor(", "", gsub(")", "", trimws(fnames), fixed = TRUE), fixed = TRUE))

			prop <- as.data.frame(matrix(NA,nrow=length(fnames),ncol=3))
			colnames(prop) <- c('factor','SS','prop')
			prop[,1] <- fnames

			#extract sum of squares
			prop[1, "SS"]  <- sfit[["Error: Within"]][[1]][, "Sum Sq"]
			prop[-1, "SS"] <- sfit[["Error: realizations"]][[1]][, "Sum Sq"]

			#calculate proportion of variation explained
			prop[, 'prop'] <- as.numeric(prop[, "SS"]) / sum(as.numeric(prop[, "SS"]))
		} else {
			prop <- NULL
		}
		
		saveRDS(prop, file = ftemp2)
	}
	prop
})

calc_region_partition2 <- compiler::cmpfun(function(j, ir) {
	ftemp2 <- file.path(dir.tables, "Partition", paste0("Partition_Props_region", ir, "_var", j, "_temp2.rds"))
	if (action == "continue" && file.exists(ftemp2)) {
		prop <- readRDS(file = ftemp2)
	} else {
		print(paste(Sys.time(), ": Partition of region:", ir, "; variable:", j, "of", length(variables)))

		# with 320,000 cases, calling aov() uses about 31 GB of memory
		# instead: repeat random subsample of runs and realizations
		
		tmp_dat <- mdata[, c(variables[j], colnames(runRequests))]
		colnames(tmp_dat) <- c('y', colnames(runRequests))
		irun <- tmp_dat[, "run"]
		ireals <- tmp_dat[, "realizations"]
		tmp_dat$realizations <- apply(tmp_dat[,factors], 1, paste,collapse=".")

		aov_formula1 <- y ~ factor(types)*factor(errors)*factor(models)*factor(mlevels)
		fnames <- attr(terms.formula(aov_formula1), "term.labels")
		fnames <- c("realizations", gsub("factor(", "", gsub(")", "", trimws(fnames), fixed = TRUE), fixed = TRUE), "Residuals")
		prop <- array(NA, dim = c(repeatsN_partition_subsample, length(fnames), 2),
							dimnames = list(NULL, fnames, c('SS', 'prop')))

		set.seed(127)

		for (im in seq_len(repeatsN_partition_subsample)) {
# sample(., size = x)
#	- x = 1: subset size of 14240 and uses 2.1 GB during the aov() call, but within-SS is not estimable
# 	- x = 3: subset size of 41760 and uses 4.6 GB during the aov() call
#	- x1 = 4, x2 = 4: subset size of 55040 and uses 5.7 GB during the aov() call
#	- x1 = 5, x2 = 4: subset size of 60800 and uses 8.0 GB during the aov() call

			tmp_sub <- irun %in% sample(x = evaluationRepeatsN, size = 4) |
						ireals %in% sample(x = presenceRealizationsN, size = 4)
#			tmp_sub <- irun %in% sample(x = evaluationRepeatsN, size = 0) |
#						ireals %in% sample(x = presenceRealizationsN, size = 1)

			is_balanced <- all(diff(table(tmp_dat[tmp_sub, head(factors, n = -1)])) == 0)
			if (!is_balanced) next

			aov_sum <- summary(aov(as.formula(paste("y ~", as.character(aov_formula1)[3], "+ Error(realizations)")),
									data = tmp_dat, subset = tmp_sub))		

			if (length(aov_sum) == 2) {
				#extract sum of squares
				prop[im, 1, "SS"]  <- aov_sum[["Error: Within"]][[1]][, "Sum Sq"]
				prop[im, -1, "SS"] <- aov_sum[["Error: realizations"]][[1]][, "Sum Sq"]

				#calculate proportion of variation explained
				prop[im, , 'prop'] <- as.numeric(prop[im, , "SS"]) / sum(as.numeric(prop[im, , "SS"]))
			} else {
				stop(str(prop))
			}
		}
		
		#coverage <- sum(tmp_sub) / length(tmp_sub)
		
		saveRDS(prop, file = ftemp2)
	}
	prop
})



eval3.SDMs <- compiler::cmpfun(function(i){
	ftemp <- get_temp_fname(runEvals[i, ], runEvalIDs[i])
	if (!file.exists(ftemp)) {
		print(paste(Sys.time(), "eval3.SDMs:", i, paste(runEvals[i, ], collapse = "-")))

		type <- runEvals[i, "types"]
		error <- runEvals[i, "errors"]
		xt <- with(runRequests, which(types == type & errors == error & mlevels == runEvals[i, "mlevels"] & models == runEvals[i, "models"]))
		xt <- xt[xt %in% runIDs]
		
		bsub <- lapply(xt, function(x) try(readRDS(file = get_temp_fname(runRequests[x, ], runRequestIDs[x])), silent = TRUE))
	  	ibad <- sapply(bsub, function(x) inherits(x, "try-error"))
		if (sum(ibad) > 0) bsub <- bsub[!ibad]
		if (length(bsub) == 0) stop("bsub does not contain data")

		ids <- matrix(unlist(t(sapply(bsub, FUN=function(l) l$runID))[, c("realizations", "run")]), ncol=2, byrow=FALSE, dimnames=list(NULL, c("realizations", "run")))

		#Prepare result container
		bevalM <- vector(mode="list", length=4)
		names(bevalM) <- c("evalID", "Eval", "Deviance", "Proj") #removed variable importance for now
		bevalM$evalID <- runEvals[i, ]

		#Evaluate models based on evaluation datasplits and realizations
		stat.methods <- c('Testing.data', 'Evaluating.data', 'Cutoff', 'Sensitivity', 'Specificity')
		temp.Eval <- array(NA,dim=c(length(eval_disc.methods), length(stat.methods), temp <- apply(ids, 2, max)),
						 dimnames=list(eval_disc.methods, stat.methods, eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))
		temp.Dev <- array(NA, dim=temp)

		for (j in 1:nrow(ids)){
			temp.Eval[,, ids[j, "realizations"], ids[j, "run"]] <- bsub[[j]]$cutoff[eval_disc.methods, stat.methods]
			temp <- bsub[[j]]$DEV
			temp.Dev[ids[j, "realizations"], ids[j, "run"]] <- if (is.null(temp)) NA else temp
		}

		bevalM$Eval$mean <- apply(temp.Eval, MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
		bevalM$Eval$sd   <- apply(temp.Eval, MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
		bevalM$Eval$quantiles   <- apply(temp.Eval, MARGIN=c(1,2),FUN=quantile, probs=quantile.probs, type=8, na.rm=TRUE)
		bevalM$Deviance$mean <- mean(temp.Dev, na.rm=TRUE)
		bevalM$Deviance$sd <- sd(temp.Dev, na.rm=TRUE)
		bevalM$Deviance$quantiles <- quantile(temp.Dev, probs=quantile.probs, type=8, na.rm=TRUE)

		#Project models onto all regions and take differences between projected and base region
		bevalM$Proj <- vector(mode="list", length=length(regions))
		stat.methods <- c("Testing.data", "Cutoff", "Sensitivity", "Specificity")

		for(ir in c(baseRegion, seq_along(regions)[-baseRegion])){
			temp.Eval <- array(NA,dim=c(length(eval_disc.methods), length(stat.methods), temp <- apply(ids, 2, max)),
							   dimnames=list(eval_disc.methods, stat.methods, eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))
			temp.Prob <- array(NA,dim=c(length(eval_cont.methods), temp),
							   dimnames=list(eval_cont.methods, eval(parse(text=paste(names(temp)[1], "=1:", temp[1]))), eval(parse(text=paste(names(temp)[2], "=1:", temp[2])))))

			#Evaluate projections based on complete dataset
			for (j in 1:nrow(ids)) {
			  etemp <- get_region_eval(xt[j], ir, stat.methods)
			  temp.Eval[,, ids[j, "realizations"], ids[j, "run"]] <- t(etemp[, eval_disc.methods])
			  temp.Prob[, ids[j, "realizations"], ids[j, "run"]] <- etemp[1, eval_cont.methods]
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
	
		
		saveRDS(bevalM, file = ftemp)
	}
  	
	i
})



## Functions to estimate model complexity
get_complexity <- compiler::cmpfun(function(i) {	
	bresM <- try(readRDS(file = get_temp_fname(runRequests[i, ], runRequestIDs[i])), silent = TRUE)
	res <- rep(NA, 3)
	
	if (!inherits(bresM, "try-error")) {
		res[1] <- 	if(runRequests[i, "models"] == "GAM") {
						sum(bresM$m$edf)
					} else if(runRequests[i, "models"] == "GLM") {
						with(bresM$m, df.null - df.residual + 1)
					} else {
						NA
					}
		res[2] <- bresM$comp_time
		res[3] <- if (is.null(bresM$comp_time_total)) res[2] else bresM$comp_time_total
	}

	res
})

plot_complexity <- function(preds = NULL, y, ylab, fname) {
	dir.create(dir_temp <- file.path(dir.figs, "Complexity"), showWarnings = FALSE)
	ftemp <- file.path(dir_temp, fname)
	
	if (is.null(preds)) preds <- colnames(runEvals)
	ylim <- c(0, max(y, na.rm = TRUE))
	cats <- unique(runRequests[, preds])
	ncats <- nrow(cats)
	
	# determine figure and axis label size
	pwidth <- max(5, 10 / 40 * ncats)
	pheight_target <- 7
	png(file = ftemp, width = pwidth, height = pheight_target, units="in", res=600)
	xcex <- 0.85
	nchar_xlab <- max(strwidth(apply(cats, 1, paste, collapse = "_"), units = "inches", cex = xcex))
	pheight_ratio <- min(6, max(3, pheight_target / nchar_xlab))
	pheight <- round(nchar_xlab * pheight_ratio * 2) / 2
	dev.off()
	
	# plot
	png(file = ftemp, width = pwidth, height = pheight, units="in", res=600)
	op <- par(mar = c(0.5 + ceiling(nchar_xlab / par("cin")[2]), 4, 1, 1))
	tmp <- boxplot(y ~ apply(runRequests[, preds],1,paste,collapse="_"), ylim = ylim, axes=FALSE,frame.plot=TRUE)
	axis(2)
	axis(1,at = 1:length(tmp$names), labels = tmp$names,las=3, cex.axis = xcex)
	mtext(ylab,side=2,line=2.5)
	par(op)
	dev.off()
}

calc_NT1 <- function(refdat, prodat) {
	stopifnot(identical(colnames(refdat), colnames(prodat)))
	
	#——————————————————————–#
	# NT1 – UNIVARIATE EXTRAPOLATION
	#——————————————————————–#
	# Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014. Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity and Distributions 20:1147-1159.
	# code based on comment by Matthew Bayly to https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/

	range_ref <- t(matrixStats::colRanges(refdat))
	#dimnames(range_ref) <- list(c("min", "max"), colnames(refdat))
	range_ref_arr <- array(range_ref, dim = c(dim(range_ref), nrow(prodat)), dimnames = list(c("min", "max"), colnames(refdat), NULL))

	diffs_ref <- matrixStats::colDiffs(range_ref)
	#colnames(diffs_ref) <- colnames(refdat)
	diffs_ref_arr <- matrix(diffs_ref, nrow = nrow(prodat), ncol = ncol(prodat), byrow = TRUE)

	iud <- array(0, dim = c(dim(prodat), 3))
	iud[, , 2] <- prodat - t(range_ref_arr["min", ,])
	iud[, , 3] <- t(range_ref_arr["max", ,]) - prodat

	UDs <- apply(iud, 1:2, min) / diffs_ref_arr
	NT1 <- rowSums(UDs)
}

calc_NT2 <- function(refdat, prodat) {
	stopifnot(identical(colnames(refdat), colnames(prodat)))

	#——————————————————————–#
	# NT2 – MULTIVARIATE EXTRAPOLATION
	#——————————————————————–#
	# Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. 2014. Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity and Distributions 20:1147-1159.
	# code modified from on https://pvanb.wordpress.com/2014/05/13/a-new-method-and-tool-exdet-to-evaluate-novelty-environmental-conditions/

	# Calculate the average and covariance matrix of the variables 
	# in the reference set
	ref.av  <- colMeans(refdat, na.rm=TRUE)
	ref.cov <- var(refdat, na.rm=TRUE)
 
	# Calculate the mahalanobis distance of each raster 
	# cell to the environmental center of the reference 
	# set for both the reference and the projection data 
	# set and calculate the ratio between the two.
	mah.ref <- mahalanobis(x = refdat, center = ref.av, cov = ref.cov)
	mah.pro <- mahalanobis(x = prodat, center = ref.av, cov = ref.cov)
	mah.max <- max(mah.ref[is.finite(mah.ref)])
	NT2 <- mah.pro / mah.max
}

plot_extrapolation <- function(NT1rast, NT2rast, file) {
#——————————————————————–#
# Plot the extrapolation rasters
#——————————————————————–#

	nt2_max <- ceiling(cellStats(NT2rast, max))
	n <- 0:255
	d <- nt2_max / length(n)
	n1 <- round(1 / d)
	zlim1 <- c(cellStats(NT1rast, min), 0)
	zlim2 <- c(0, nt2_max)
	blegend <- c(-95, -93.5, 30, 48)

	cols_sim <- colorRampPalette(c("gray", "darkslategray1"))
	cols_diss <- colorRampPalette(c("gold1", "orange", "red", "darkred", "purple"))
	cols_NT1 <- cols_NT2 <- list()
	cols_NT1[["added_below"]] <- cols_NT1[["added_above"]] <- FALSE
	cols_NT1[["colors_label"]] <- rev(c("gray", cols_diss(n = length(n) - 1)))
	cols_NT2[["added_below"]] <- cols_NT2[["added_above"]] <- FALSE
	cols_NT2[["colors_label"]] <- c(cols_sim(n = n1), cols_diss(n = length(n) - n1))


	cex <- 1
	npanelsX <- 2; npanelsY <- 1
	h.panel <- 2.5; h.edge_lo <- 0.2; h.edge_up <- 0.05
	w.panel <- 1.8; w.edge_left <- 0.1; w.edge_right <- 0.0

	png(height = h.edge_lo + h.panel * npanelsX +  h.edge_up,
		width = w.edge_left + w.panel * (1 + npanelsY) + w.edge_right, units = "in",
		res = 600, file = file)

	panels <- matrix(0, nrow = 1 + npanelsX + 1, ncol = 1 + npanelsY + 1, byrow=FALSE)
	panels[-c(1, nrow(panels)), -c(1, ncol(panels))] <- 1:(npanelsX * npanelsY)
	layout(panels,
			heights = c(h.edge_up, rep(h.panel, times = npanelsX), h.edge_lo),
			widths = c(w.edge_left, rep(w.panel, times = npanelsY), w.edge_right))
	op <- par(mgp = c(1, 0, 0), mar = c(0.5, 0.5, 0.5, 0.2), tcl = 0.3, cex = cex, xpd = NA)

		raster::image(NT1rast, asp = 1, zlim = zlim1, col = cols_NT1[["colors_label"]], xlab = "", ylab = "", axes = FALSE)
		add_legend(zlim = zlim1, zextreme = zlim1, col_desc = cols_NT1, grid = NT1rast, box = blegend, cex = 0.85)
		raster::image(rbaseRegion, col = "darkgray", add = TRUE)
		plot(allborders, add = TRUE)
		axis(2)
		mtext(text = "(a)", line = -1, adj = 0.05, font = 2)

		raster::image(NT2rast, asp = 1, zlim = zlim2, col = cols_NT2[["colors_label"]], xlab = "", ylab = "", axes = FALSE)
		add_legend(zlim = zlim2, zextreme = zlim2, col_desc = cols_NT2, grid = NT2rast, box = blegend, cex = 0.85)
		image(rbaseRegion, col = "darkgray", add = TRUE)
		plot(allborders, add = TRUE)
		axis(1); axis(2)
		mtext(text = "(b)", line = -1, adj = 0.05, font = 2)
	par(op)
	dev.off()
}

# Create and plot a continuous color legend
add_legend <- function(zlim, zextreme, col_desc, grid, box=c(-100, -97, -50, -10), whitebox=TRUE, horiz=FALSE, signed=1, fun_inv_ens=NULL, srt=90, cex=1){
	if(is.null(fun_inv_ens)) fun_inv_ens <- function(x) x

	# Color ramp
	zr <- raster(xmn=box[1], xmx=box[2], ymn=box[3], ymx=box[4], crs=projection(grid), resolution=res(grid), vals=NULL)
	zr[] <- if(horiz) rep(1:dim(zr)[2], times=dim(zr)[1]) else rep(dim(zr)[1]:1, each=dim(zr)[2])
	if(whitebox) raster::image(zr, col="white", add=TRUE)
	raster::image(zr, col=col_desc$colors_label, add=TRUE)

	# Labels	
	atz <- pretty(zlim, n=6) #generate default labels	
	if((temp1z <- sum(temp2z <- atz <= zlim[1])) > 0){#adjust lowest label to limit value
		atz[tail(which(temp2z), n=1)] <- zlim[1]
		if(temp1z > 1) atz <- atz[-head(which(temp2z), n=-1)]
	}
	if((temp1z <- sum(temp2z <- atz >= zlim[2])) > 0){#adjust highest label to limit value
		atz[head(which(temp2z), n=1)] <- zlim[2]
		if(temp1z > 1) atz <- atz[-tail(which(temp2z), n=-1)]
	}
	datz <- diff(atz)
	temp <- 0.8 / (if(abs(srt) > 45) cex else 1)
	if(length(datz) > 0 && length(temp <- which(datz < temp * max(datz))) > 0){#remove labels if too close together, but not limits and not 0
		id_remove <- findInterval(x=1+temp, vec=1:length(atz), all.inside=TRUE)
		if(length(temp <- which(0 == atz[id_remove])) > 0) id_remove <- id_remove[-temp]
		if(length(id_remove) > 0) atz <- atz[-id_remove]
	}
	if(length(datz) == 0) atz <- zextreme
	ltxt <- prettyNum(signif(signed * fun_inv_ens(atz), 2))
	ltext_extreme <- prettyNum(signif(signed * fun_inv_ens(zextreme), 2))
	
	# Tick position
	ext <- extent(zr)
	if(horiz){
		xmin_orig <- ext@xmin
		if(col_desc$added_below) xmin_orig <- xmin_orig + col_desc$ncol_label_added / length(col_desc$colors_label) * (ext@xmax - ext@xmin)
		xmax_orig <- ext@xmax
		if(col_desc$added_above) xmax_orig <- xmax_orig - col_desc$ncol_label_added / length(col_desc$colors_label) * (ext@xmax - ext@xmin)
		xs_orig <- (temp <- atz / (max(atz) - min(atz)) * (xmax_orig - xmin_orig)) + (xmin_orig - min(temp))
	} else {
		ymin_orig <- ext@ymin
		if(col_desc$added_below) ymin_orig <- ymin_orig + col_desc$ncol_label_added / length(col_desc$colors_label) * (ext@ymax - ext@ymin)
		ymax_orig <- ext@ymax
		if(col_desc$added_above) ymax_orig <- ymax_orig - col_desc$ncol_label_added / length(col_desc$colors_label) * (ext@ymax - ext@ymin)
		ys_orig <- (temp <- atz / (max(atz) - min(atz)) * (ymax_orig - ymin_orig)) + (ymin_orig - min(temp))
	}
	
	# Draw ticks
	lwd_seg <- max(0.5, min(1, cex)) * par("lwd")
	if(horiz){
		segments(x0=xs_orig, x1=xs_orig, y0=ext@ymax - (ext@ymax - ext@ymin) / 3, y1=ext@ymax, lwd=lwd_seg)
		if(col_desc$added_below) segments(x0=ext@xmin, x1=ext@xmin, y0=ext@ymin, y1=ext@ymax, lwd=lwd_seg)
		if(col_desc$added_above) segments(x0=ext@xmax, x1=ext@xmax, y0=ext@ymin, y1=ext@ymax, lwd=lwd_seg)
	} else {
		segments(x0=ext@xmax - (ext@xmax - ext@xmin) / 3, x1=ext@xmax, y0=ys_orig, y1=ys_orig, lwd=lwd_seg)
		if(col_desc$added_below) segments(x0=ext@xmin, x1=ext@xmax, y0=ext@ymin, y1=ext@ymin, lwd=lwd_seg)
		if(col_desc$added_above) segments(x0=ext@xmin, x1=ext@xmax, y0=ext@ymax, y1=ext@ymax, lwd=lwd_seg)
	}
	
	# Write tick labels
	if(horiz){
		if(abs(srt) > 45){
			adj <- c(0.5, 1.3)
			ly <- ext@ymax
		} else {
			adj <- c(0.5, NA)
			ly <- ext@ymax + strheight(ltxt, units="user", cex=cex * 1.05)
		}
		text(x=xs_orig, y=ly, labels=ltxt, srt=srt, adj=adj, cex=cex, xpd=TRUE)
		if(col_desc$added_below) text(x=ext@xmin, y=ext@ymax, labels=ltext_extreme[1], srt=90, adj=c(1.3, 1), cex=cex, xpd=TRUE)
		if(col_desc$added_above) text(x=ext@xmax, y=ext@ymax, labels=ltext_extreme[2], srt=90, adj=c(-0.5, 1), cex=cex, xpd=TRUE)
	} else {
		if(abs(srt) > 45){
			adj <- c(0.5, 1.3)
			lx <- ext@xmax
		} else {
			adj <- c(1, NA)
			lx <- ext@xmax + (if(ext@xmax > 0) -1 else +1) * max(strwidth(ltxt, units="user", cex=cex * if(cex < 0.5) 1.5 else 1.05))
		}
		text(x=lx, y=ys_orig, labels=ltxt, srt=srt, adj=adj, cex=cex, xpd=TRUE)
		if(col_desc$added_below) text(x=ext@xmax, y=ext@ymin, labels=ltext_extreme[1], srt=0, adj=c(1, 1.3), cex=cex, xpd=TRUE)
		if(col_desc$added_above) text(x=ext@xmax, y=ext@ymax, labels=ltext_extreme[2], srt=0, adj=c(1, -0.5), cex=cex, xpd=TRUE)
	}
		
	invisible(0)
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
  points(XY, pch=16, cex=ifelse(obs > 0, 0.06 + 0.25*obs, 0), col = adjustcolor("black", alpha.f = 0.4))
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


arrange.response.plot2 <- compiler::cmpfun(function(show.variables, modelName, xdat, ydat) {
	list.out <- NULL
	for (i in seq_along(show.variables)) {
		# 5. Storing results
		list.out[[i]] <- data.frame(xdat[[i]][["plot_data"]], ydat[[i]])
		colnames(list.out[[i]]) <- c(show.variables[i], modelName)
	}
	names(list.out) <- show.variables
	
	list.out
})



pred.response.plot2 <- compiler::cmpfun(function(data_species, data_env, orig.variables, scaled.variables=NULL, centerMeans, fixed.var.metric = 'mean'){
##biomod2::response.plot2:
#	- doesn't find loaded models even though they were loaded properly
#	--> our own version  loads the model instead of assuming that the models are loaded
#	- doesn't account for 'coupled' variables like MAT and MATScaled
#	--> our own version co-varies 'coupled' variables properly
##deleted arguments: do.bivariate = FALSE, save.file="no", name="response_curve", ImageSize=480, plot=FALSE

	# 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
	if(is.null(scaled.variables)){
		show.variables <- orig.variables
		isScaled <- FALSE
	} else {
		show.variables <- scaled.variables
		isScaled <- TRUE
	}

	nb.pts <- 100
	if(is.null(data_species)){
		data_species <- rep(1,nrow(data_env))
	} else {
		data_species[data_species!=1 | is.na(data_species)] <- 0
	}


	# 2. build function outputs
	factor_id <- which(sapply(data_env,is.factor))
	list.out <- list()

	# Create a ranged data table
	ref_table <- data_env[1,,drop=F]
	rownames(ref_table) <- NULL

	for(i in 1:ncol(data_env)){
		temp <- data_env[data_species==1,i]
		if(is.numeric(data_env[,i])){
			ref_table[,i] <- switch(fixed.var.metric,
								  mean = mean(temp),
								  median = median(temp),
								  min = min(temp),
								  max = max(temp))
		} else{
			# return the majoritary class
			sum_level <- summary(temp)
			ref_table[,i] <- names(sum_level)[sum_level==max(sum_level)]
		}
	}

	for (i in seq_along(show.variables)) {
		# creating Tmp data
		if (is.factor(data_env[, show.variables[i]])) {
			pts.tmp <- as.factor(levels(data_env[, show.variables[i]]))
			pts.tmp.orig <- NULL
		} else {
			pts.tmp <- seq(min(data_env[, show.variables[i]]), max(data_env[, show.variables[i]]), length.out = nb.pts)
			temp <- grep(iorig <- sub("Scaled", "", show.variables[i]), names(centerMeans))
			pts.tmp.orig <- if(isScaled) pts.tmp + centerMeans[temp] else NULL
		}

		Data.r.tmp <- eval(parse(text=paste("cbind(", show.variables[i], "=pts.tmp,ref_table[,-which(colnames(ref_table)==show.variables[i]),drop=F])",sep="")))
		Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
		if(length(factor_id)){
			for(f in factor_id){
				Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(data_env[,f]))
			}
		}
		if(!is.null(pts.tmp.orig)) Data.r.tmp[, iorig] <- pts.tmp.orig

		list.out[[i]] <- list(new_data = Data.r.tmp, plot_data = pts.tmp)
	}
	names(list.out) <- show.variables

	list.out
})


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
              border=NA, col = adjustcolor("blue", alpha.f = 0.3))
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, 1, 2], rev(aggR2s[, 2, 1, 4])), 
              border=NA, col = adjustcolor("orange", alpha.f = 0.9))
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
              border=NA, col = adjustcolor("blue", alpha.f = 0.3))
      polygon(x=c(aggR2s[, 2, 2, 2], rev(aggR2s[, 2, 2, 4])), 
              y=c(x, rev(x)), 
              border=NA, col = adjustcolor("orange", alpha.f = 0.9))
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
              border=NA, col = adjustcolor("blue", alpha.f = 0.3))
      polygon(x=c(x, rev(x)), 
              y=c(aggR2s[, 2, env.ind, 2], rev(aggR2s[, 2, env.ind, 4])), 
              border=NA, col = adjustcolor("orange", alpha.f = 0.9))
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
    
    FittedData <- ObsData <- matrix(NA, nrow=length(bsub[[1]]$Proj[[ir]]), ncol=nrow(ids))
    respCurvePreds <- vector(mode="list", length=nrow(ids))
    for(rr in 1:nrow(ids)){
      FittedData[, rr] <- bsub[[rr]]$Proj[[ir]]
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

work <- compiler::cmpfun(function() {
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
})

