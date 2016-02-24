###############################################################################
#
# Bell, D. M., and D. R. Schlaepfer. Impacts of the data-generating processes and species distribution model complexity on ecological fidelity and global change predictions.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
###############################################################################


## Actions
action <- "continue" # "new", restart all computations; "continue", attempts account for already completed simulations
do.ExampleForDropbox <- FALSE #set this to TRUE to access/write to 'Example' subset on dropbox
do.SDMs <- FALSE					# '20160223' with 320,000 runs: 2016/02/24 completed in 237 core-hours
do.Partition <- TRUE				
do.Complexity <- FALSE
do.Evaluation <- FALSE
do.EvaluationSummary <- FALSE
do.Figures <- FALSE

##
## Load libraries
libraries <- c("reshape2", "mgcv", "raster", "parallel", "mvtnorm", "randomForest", "gbm", "maxent")
temp <- lapply(libraries, FUN=require, character.only=TRUE)

date.run <- "20160223" #label for output folders (20140228; 20140304; 20140314; 20140320; 20140321; 20140627; 20150130)
if (do.ExampleForDropbox) date.run <- paste0(date.run, "_Example")

## Directories
computer <- "Daniel" # "Dave"

if (computer == "Dave") {
	dir.prj <- "E:/Work_Share_Dave_SDM&ResponseCurves"
	path.functions <-"E:/Dropbox/Work_Share_Dave_SDM&ResponseCurves/BIOMOD2" 
	dir.big <- dir.prj
} else if (computer == "Daniel") {
	dir.prj <- "~/Dropbox/Work_Stuff/2_Research/200907_UofWyoming_PostDoc/Product21_SDM_AsymmetryResponseCurve/3_Simulations"
	path.functions <- file.path(dir.prj, "code")
	dir.big <- "~/Downloads/Product21_SDM_AsymmetryResponseCurve"
} else stop(computer, " is not implemented")

dir.dat <- file.path(normalizePath(dir.prj), "inst", "extdata")
dir.in <- file.path(dir.dat, "csv.files")
dir.gis <- file.path(dir.dat, "rasters")
dir.res <- file.path(normalizePath(dir.big), "Output")
dir.res2 <- file.path(normalizePath(dir.prj), "Output")
dir.sdm <- file.path(dir.res, paste0("SDMs_", date.run))
dir.maps <- file.path(dir.res2, paste0("Maps_", date.run))
dir.figs <- file.path(dir.res2, paste0("Plots_", date.run))
dir.tables <- file.path(dir.res2, paste0("Tables_", date.run))
temp <- lapply(c(dir.sdm, dir.maps, dir.figs, dir.tables), function(x) dir.create(path = x, showWarnings=FALSE, recursive=TRUE))


## Settings
num_cores <- 22
parallel_backend <- "parallel" # "parallel" (here, uses sockets on windows and forks on unix systems) or "mpi" (requiring a MPI installed)
fflag <- paste0("v2_", date.run)
filename.runRequests <- paste0("runSetup_", fflag, ".RData")
filename.saveSDMs <- paste0("SDMs_", fflag, ".rds")
filename.saveRunIDs <- paste0("runIDs_", fflag, ".rds")
filename.saveEvals <- paste0("SDMs_Evaluations_", fflag, ".RData")
filename.saveParts <- paste0("SDMs_VarPartition_", fflag, ".rds")
filename.saveComplexity <- paste0("SDMs_ModelComplexity_", fflag, ".RData")
filename.saveTypeData <- paste0("TypeData_", fflag, ".RData")

baseRegion <- 2 #NorthWest
regions <- 1:4
types <- c("AIF", "SCT", "SIF", "SIT")
mlevels <- list(woInt=c("linear", "squared"), wInt=c("linear", "squared", "interaction"))
sdm.models <- c("GLM", "GAM", "MaxEnt", "RF", "BRT")
eval.methods <- c('TSS','ROC','KAPPA')
errors <- c("binom","binom+res","spatial","spatial+res")
evaluationRepeatsN <- switch(EXPR=paste0("v_", date.run),
                            v_20140228=50,
                            v_20140304=10,
                            v_20140314=10,
                            v_20140320=2,
                            v_20140321=50,
                            v_20140627=25,
                            v_20150130=25,
                            v_20160209=4,
                            v_20160223=50)
presenceRealizationsN <- switch(EXPR=paste0("v_", date.run),
                                v_20140228=40,
                                v_20140304=5,
                                v_20140314=40,
                                v_20140320=2,
                                v_20140321=40,
                                v_20140627=25,
                                v_20150130=25,
                                v_20160209=4,
                            	v_20160223=40)
predictorsN <- 10
equalSamples <- FALSE #if TRUE, ensures that subsamples have the same number of presences as absences


## Define set of SDM runs
if (action == "continue" && file.exists(ftemp <- file.path(dir.res, filename.runRequests))) {
	load(file = ftemp) #runRequests, runRequestIDs, runEvals, runEvalIDs
} else {
	if (do.ExampleForDropbox) {
		#runRequests <- runRequests[sample(x=nrow(runRequests), size=5), ]
		runRequests <- rbind(runRequests[runRequests$models == "GLM" & runRequests$mlevels == "wInt" & runRequests$types=="AIF" & runRequests$run == 1, ][1:5, ],
							 runRequests[runRequests$models == "GAM" & runRequests$mlevels == "woInt" & runRequests$types=="AIF" & runRequests$run == 1, ][1:5, ])
	} else {
		runRequests <- expand.grid(sdm.models, names(mlevels), types, errors, 1:presenceRealizationsN, 1:evaluationRepeatsN, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)
		colnames(runRequests) <- c("models", "mlevels", "types", "errors", "realizations", "run")
	}
	runRequestIDs <- apply(runRequests, MARGIN = 1, function(x) paste0("SDM_", paste(trimws(x), collapse = "_")))
	
	runEvals <- unique(runRequests[, c("models", "mlevels", "types", "errors")]) #get the unique combinations of model, type, and mlevel
	runEvalIDs <- apply(runEvals, MARGIN = 1, function(x) paste0("Eval_", paste(trimws(x), collapse = "_")))
	
	save(runRequests, runRequestIDs, runEvals, runEvalIDs, file = ftemp)
}

runFolders <- apply(unique(runEvals[, c("models", "types")]), 1, paste, collapse = "_")
# nrow(runRequests) / length(runFolders) # we shouldn't have too many files per folder
temp <- lapply(runFolders, function(x) dir.create(path = file.path(dir.sdm, x), showWarnings=FALSE, recursive=TRUE))


######
print(Sys.time())
print(sessionInfo())

####### Set up parallel environment
#if(!interactive()){
num_cores <- min(num_cores, parallel::detectCores() - 2)

if(identical(parallel_backend, "mpi")){
  stopifnot(require("Rmpi"))
  libraries <- c(libraries, "Rmpi")
  
  .Last <- function() { #Properly clean up mpi before quitting R (e.g., at a crash)
    if(is.loaded("mpi_initialize") && exists("mpi.comm.size")){
      if (mpi.comm.size(1) > 0) mpi.close.Rslaves()
      .Call("mpi_finalize")
    }
  }
  
  mpi.spawn.Rslaves(nslaves=num_cores)
  
  exportObjects <- function(allObjects) {
    for(obj in 1:length(allObjects)) {
      bcast.tempString <- allObjects[obj]
      bcast.tempValue <- try(eval(as.name(allObjects[obj])))
      if(!inherits(bcast.tempValue, "try-error")){
        mpi.bcast.Robj2slave(bcast.tempString)
        mpi.bcast.Robj2slave(bcast.tempValue)
        mpi.bcast.cmd(cmd=try(assign(bcast.tempString, bcast.tempValue)))
      } else {
        print(paste(obj, bcast.tempString, "not successful"))
      }
    }
  }
  
} else if(identical(parallel_backend, "parallel")){
  
	cl  <- if (.Platform$OS.type == "unix") {
				makeCluster(num_cores, type = "FORK", outfile = "log_sdm.txt")
			} else if (.Platform$OS.type == "windows") {
				makeCluster(num_cores, type = "PSOCK", outfile = "log_sdm.txt")
			} else {
				stop("Running this code on this type of platform is currently not implemented.")
			}
  
  .Last <- function() { #Properly clean up cluster before quitting R (e.g., at a crash)
    if(exists("stopCluster") && exists("cl")){
      stopCluster(cl)	#clean up cluster
    }
  }
} else {
  stop("No available parallel backend specified.")
}	
#}

####### Functions

source(file.path(path.functions, "ResponseExtrapolation_Functions_Bell_and_Schlaepfer.R"))

################
## Read data from files once and generate observations
ftemp <- file.path(dir.res, filename.saveTypeData)
if (action == "continue" && file.exists(ftemp)) {
  print(paste(Sys.time(), ": Loading generated data"))
  load(ftemp)
  
} else {
  print(paste(Sys.time(), ": Data generation started"))

  tname <- paste(rep(types,times=length(errors)),rep(errors,each=length(types)),sep="_")
  tmat <- cbind(type=rep(types,times=length(errors)),var=rep(errors,each=length(types)))
  
  typeData <- varData <- obsData <- probData <- vector("list", length=length(tname))
  names(typeData) <- names(varData) <- names(obsData) <- names(probData) <- tname
  
  sigma <- logit(.5) - logit(.4) #standard deviation set so that sd = +/- 0.10 when p = 0.5
  
  for(tp in 1:length(tname)){
    typeData[[tname[tp]]] <- get_TypeData_FromFile(type=tmat[tp,"type"], center=TRUE, centerBasedOnRegionIDs=baseRegion)
    varData[[tname[tp]]]  <- tmat[tp,"var"]
    probData[[tname[tp]]] <- typeData[[tname[tp]]]$dat[,'prob']
    
    if(tp == 1) 
      dtmp <- as.matrix(dist(typeData[[tname[tp]]]$dat[,c("x","y")]),nrow = nrow(typeData[[tname[tp]]]$dat[,c("x","y")]))
    
      w <- probData[[tname[tp]]]*0
      for(j in 1:length(w)) {
       
        kk <- which(dtmp[,j] != 0 & dtmp[,j] < 1)
        
        w[j] <- sum(dtmp[kk,j] * probData[[tname[tp]]][kk]) / sum(dtmp[kk,j])
      }

    
    obsData[[tname[tp]]]  <- calc_ObservationsFromProbabilities(probs=probData[[tname[tp]]], 
                                                          N=presenceRealizationsN,
                                                          VAR=varData[[tname[tp]]],
                                                          sigma = sigma,
                                                          w = w)
    if (is.null(obsData[[tname[tp]]]))
    	stop("Call to 'calc_ObservationsFromProbabilities' with variance type = ", varData[[tname[tp]]], " failed for tname[", tp, "] = ", tname[tp])
  }
  
  rm(list = c("dtmp", "kk"))
  
  # Climate Data -- one matrix with x, y, region, and climate data
  climData <- typeData[['AIF_binom']]$dat[,-grep('prob',colnames(typeData[['AIF_binom']]$dat))]
  # centerMeans Data
  centerMeansData <- typeData[['AIF_binom']]$centerMeans
  
  rm(typeData)
  save(probData, obsData, climData, centerMeansData, file=ftemp)

  print(paste(Sys.time(), ": Data generation ended"))
}



## Build SDMs
if (do.SDMs) {
  print(paste(Sys.time(), ": SDMs started"))
    
  list.export <- c("libraries", "climData", "obsData", "probData", "centerMeansData", 
                   "runRequests", "runRequestIDs", "baseRegion", "regions", "mlevels", "sdm.models", 
                   "eval.methods", "predictorsN", "make.SDM", "set_Data", "calc_sdms", 
                   "set_options", "make_prediction", "make_projection","dir.sdm", "dir.in","equalSamples",
                   "get.balanced.sample","get.cutoff", "our.response.plot2", "get_temp_fname")


	# determine which SDMs have already been calculated and stored to disk
	xt <- seq_len(nrow(runRequests))

	if (action == "continue") { 
		fdone <- list.files(dir.sdm, pattern = "SDM_", recursive = TRUE)
		
		if (length(fdone) > 0) {
			wd <- setwd(dir.sdm)
			fsize <- file.info(fdone, extra_cols = FALSE)$size

			fgood <- fsize > 1024
			fbad <- !fgood
			if (any(fbad)){
				unlink(fbad)
				fdone <- fdone[fgood]
			}
			setwd(wd)
			
			fbase <- sub(".rds", "", basename(fdone))
			ifb <- match(fbase, table = runRequestIDs, nomatch = 0)
			xt <- xt[-ifb]
		}
	}

  # call make.SDM
  if(length(xt) > 0){
	  if (identical(parallel_backend, "mpi")) {
		exportObjects(list.export)
		mpi.bcast.cmd(lapply(libraries, FUN=require, character.only=TRUE))
	
		idones <- mpi.applyLB(x=xt, fun=make.SDM)
	
		mpi.bcast.cmd(rm(list=ls()))
		mpi.bcast.cmd(gc())
	  }
	  if (identical(parallel_backend, "parallel")) {
		  clusterExport(cl, list.export)
		  clusterEvalQ(cl, lapply(libraries, FUN=require, character.only=TRUE))

		  idones <- parLapply(cl, xt, function(i) try(make.SDM(i), silent=TRUE))

		  clusterEvalQ(cl, rm(list=ls()))
		  clusterEvalQ(cl, gc())
	  }

	  badRuns <- sapply(idones, function(x) inherits(x, "try-error"))
	  if (sum(badRuns) > 0) warning("There were ", sum(badRuns), " SDM runs that threw an error")
  }
  

  # Get data from disk
  bres <- list()
  for (i in seq_len(nrow(runRequests))) {
  	bres[[i]] <- try(readRDS(file = get_temp_fname(runRequests[i, ], runRequestIDs[i])), silent = TRUE)
  }
  
  #Identify runs
  goodRuns <- sapply(bres, FUN=function(l) !inherits(l, "try-error"))
  if(any(!goodRuns)) print(bres[[which(!goodRuns)[1]]])
  if(all(!goodRuns)) stop(paste(Sys.time(), ": No SDM successful"))
  print(paste(Sys.time(), ":", sum(goodRuns), "out of", length(goodRuns), "SDMs successful"))
  bres <- bres[goodRuns]
  temp <- t(sapply(bres, FUN=function(l) l$runID))
  #runIDs = index for which row of runRequests corresponds to elements of bres
  runIDs <- na.exclude(match(apply(temp, 1, paste, collapse="_"), table=apply(runRequests, MARGIN=1, FUN=function(x) paste(trimws(x), collapse="_"))))
  
  # print(object.size(bres), units = "GB") # 25.6 Gb for length(bres) == 320000
  saveRDS(runIDs, file = file.path(dir.sdm, filename.saveRunIDs))
  temp <- try(saveRDS(bres, file = file.path(dir.sdm, filename.saveSDMs)), silent = TRUE)
  if (inherits(temp, "try-error")) warning("Saving the object 'bres' to disk failed likely because of insufficient memory: the size of 'bres' is ", print(object.size(bres), units = "GB"), " and up to twice as much memory is required: ", temp)
  
  # Model object sizes
  if (FALSE) {
	  msize <- data.frame(model = sapply(bres, function(x) x$runID$model),
							size_MB = sapply(bres, function(x) object.size(x)))
	  msize[, "size_MB"] <- msize[, "size_MB"] / 1024^2 # convert bytes -> MB
	  with(msize, boxplot(size_MB ~ model))
  
		print_elem_size <- function(x) {
			if (inherits(x, "maxent")) {
				for (it in seq_along(slotNames(x))) cat(slotNames(x)[it], format(object.size(slot(x, slotNames(x)[it])), units = "KB"), "\n")
			} else {	
				for (it in seq_along(x)) cat(names(x)[it], format(object.size(x[[it]]), units = "KB"), "\n")
			}
		}

		for (im in seq_along(sdm.models)) {
			cat(sdm.models[im], "our object: \n")
			print_elem_size(bres[[im]])
			cat(sdm.models[im], "model object: \n")
			print_elem_size(bres[[im]]$SDMs$m)
			cat("\n")
		}
  }

  rm(bres)
  print(paste(Sys.time(), ": SDMs done"))
}


if (do.Partition) {
	print(paste(Sys.time(), ": Partition started"))
  
	if (!exists("runIDs")) runIDs <- readRDS(file = file.path(dir.sdm, filename.saveRunIDs))
  
	stat.methods <- 'Testing.data'
	variables <- c('TSS','KAPPA','ROC','RMSE','MAE')
	factors <- c('types',"errors",'models','mlevels','realizations')
	
	pfile <- file.path(dir.res, filename.saveParts)
	if (action == "continue" && file.exists(pfile)) {
		part.region <- readRDS(file = pfile)
	} else {
		list.export <- c("climData", "obsData", "probData", "centerMeansData", "action",
					   "runRequests", "runRequestIDs", "baseRegion", "regions", "factors", "variables",
					   "get_region_eval", "get.cutoff", "eval.methods", "stat.methods", "rmse", "mae",
					   "get_temp_fname", "dir.tables")
		if (identical(parallel_backend, "parallel")) {
			clusterExport(cl, list.export)
		}
  
		part.region <- list()	
		for(ir in regions){
			print(paste(Sys.time(), ": Partition of region:", ir))

			#get evaluation statistics (this takes about 1.24 s per iteration)
			ftemp1 <- file.path(dir.tables, paste0("Partition_Evals_region", ir, "_temp1.rds"))
			if (action == "continue" && file.exists(ftemp1)) {
				part.mat <- readRDS(file = ftemp1)
			} else {
				part.mat <- cbind(runRequests,matrix(NA,nrow=nrow(runRequests),ncol=length(variables)))
				colnames(part.mat) <- c(colnames(runRequests),variables)
	
				if (identical(parallel_backend, "parallel")) {
					clusterExport(cl, "ir")

					idones <- parSapply(cl, runIDs, function(i) get_region_eval(i), USE.NAMES = FALSE)

					clusterEvalQ(cl, gc())
				} else stop("this is not implemented 1")
		
				temp <- t(idones)
				temp <- temp[order(temp[, 1]), ]
				part.mat[temp[, 1], variables] <- temp[, -1]

				saveRDS(part.mat, file = ftemp1)
			}
	
			#run ANOVAs to estimate partitioning of variation
			if (identical(parallel_backend, "parallel")) {
				clusterExport(cl, c("ir", "part.mat"))

				part.out <- parLapply(cl, seq_along(variables), function(j) calc_region_partition(j))

				clusterEvalQ(cl, gc())
			} else stop("this is not implemented 1")
		
			names(part.out) <- variables
			part.region[[ir]]<- part.out
		}
  
		saveRDS(part.region, file = pfile)

		if (identical(parallel_backend, "parallel")) {
			clusterEvalQ(cl, rm(list=ls()))
			clusterEvalQ(cl, gc())
		}
	}
  
	var.sort <- c(2:(length(fnames)-1),1,length(fnames))
	reg.sort <- c(2,1,3,4)
  
	part.mat <- matrix(NA,nrow=length(reg.sort)*length(variables),ncol=length(var.sort)+2)
    colnames(part.mat) <- c('region','metric',(part.region[[1]][[2]]$prop[var.sort,1]))
    part.mat[,'region'] <- rep(c('NR','SR','SW','GP'),each=length(variables))
    part.mat[,'metric'] <- rep(variables,times=length(regions))
  
	for(j in 1:length(reg.sort))
		for(i in 1:length(variables))
			part.mat[which(part.mat[,'metric'] == variables[i])[j],2+1:length(var.sort)] <- (part.region[[reg.sort[j]]][[i]]$prop[var.sort,3])
  
	write.csv(part.mat,file.path(dir.tables,"var.part.csv"),quote=FALSE)
  
	print(paste(Sys.time(), ": Partition done"))
}


##Get effective degrees of freedom
if (do.Complexity) {
	print(paste(Sys.time(), ": Model complexity started"))
	
	if (!exists("bres")) bres <- readRDS(file = file.path(dir.sdm, filename.saveSDMs))
	if (!exists("runIDs")) runIDs <- readRDS(file = file.path(dir.sdm, filename.saveRunIDs))
	
	edf <- rep(NA,length(bres))

	stop("TODO(drs): fix this: what are EDFs for RF, MaxEnt, and BRT?")
	for(ee in 1:length(edf)){
		edf[ee] <- if(runRequests[ee,"models"] == "GAM") {
						sum(bres[[ee]]$SDMs$m$edf)
					} else if(runRequests[ee,"models"] == "GLM") {
						with(bresM$SDMs$m, df.null - df.residual + 1)
					} else if (inherits(bresM$SDMs$m, "maxent")) {
						NA
					} else if (inherits(bresM$SDMs$m, "randomForest")) {
						NA
					} else if (inherits(bresM$SDMs$m, "gbm")) {
						NA
					}
	}

	save(edf, file=file.path(dir.res, filename.saveComplexity))

	png(paste(dir.figs,"DF.png",sep="/"),width=6,height=4,units="in",res=600)
	par(mar=c(8,4,1,1))
	tmp <- boxplot(edf ~ apply(runRequests[, colnames(runEvals)],1,paste,collapse="_"),axes=FALSE,frame.plot=TRUE)
	axis(2)
	axis(1,at = 1:length(tmp$names), labels = tmp$names,las=3,cex.axis=.85)
	mtext("Degrees of Freedom",side=2,line=2.5)
	dev.off()

	print(paste(Sys.time(), ": Model complexity done"))
}
stop("here done")

## Evaluate SDMs
if(do.Evaluation){
  print(paste(Sys.time(), ": Evaluation started"))
  
	if (!exists("bres")) bres <- readRDS(file = file.path(dir.sdm, filename.saveSDMs))
	if (!exists("runIDs")) runIDs <- readRDS(file = file.path(dir.sdm, filename.saveRunIDs))
  
  list.export <- c("libraries", "obsData","probData","climData","centerMeansData", "baseRegion", 
                   "eval2.SDMs", "regions", "eval.methods", "sdm.models", "rmse", "mae","get.cutoff",
                   "dir.sdm", "dir.in", "get_temp_fname")
  if(identical(parallel_backend, "mpi")){
    exportObjects(c("work", list.export))
    mpi.bcast.cmd(lapply(libraries, FUN=require, character.only=TRUE))
    mpi.bcast.cmd(work())
    
    jobs <- 1:nrow(runEvals)
    idone <- vector("list", length=length(jobs)) #Result container
    workersN <- (mpi.comm.size() - 1)
    junk <- 0
    closed_slaves <- 0
    runs.completed <- 1
    
    while(closed_slaves < workersN) {
      complete <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
      complete_info <- mpi.get.sourcetag()
      slave_id <- complete_info[1]
      tag <- complete_info[2]
      
      if(tag == 1) {# slave is ready for a task. Give it the next task, or tell it tasks are done if there are none.
        if(runs.completed <= length(jobs)){# Send a task, and then remove it from the task list
          dataForRun <- list(fun="eval",
                             i=i <- jobs[runs.completed],
                             runEval=runEvals[i, ],
                             type=type <- runEvals[i, "types"],
                             error=error <- runEvals[i, "errors"],
                             mlevel=mlevel <- runEvals[i, "mlevels"],
                             model=model <- runEvals[i, "models"],
                             runID=runID <- with(runRequests[runIDs, ], which(types == type & errors == error & mlevels == mlevel & models == model)),
                             bsub=bres[runID])
          print(paste(Sys.time(), ": Evaluation for", dataForRun$i, "started"))
          mpi.send.Robj(dataForRun, slave_id, 1)
          runs.completed <- runs.completed + 1
        } else {
          mpi.send.Robj(junk, slave_id, 2)
        }
      } else if (tag == 2) { # The message contains results
        idone[[complete$i]] <- complete$r
        print(paste(Sys.time(), ": Evaluation for", complete$i, "ended"))
      } else if (tag == 3) { # A slave has closed down.
        closed_slaves <- closed_slaves + 1
      }		
    }								
    
    mpi.bcast.cmd(rm(list=ls()))
    mpi.bcast.cmd(gc())
  }
  if(identical(parallel_backend, "parallel")){
    clusterExport(cl, c("bres", "runRequests", "runRequestIDs", "runEvals", "runEvalIDs", "runIDs", list.export)) #TODO: exporting large bres will fail
    clusterEvalQ(cl, lapply(libraries, FUN=require, character.only=TRUE))
    
    list.noexport <- (temp <- ls())[!(temp %in% list.export)]
    idone <- foreach(i=1:nrow(runEvals), .errorhandling="pass", .noexport=list.noexport) %dopar% {
      type <- runEvals[i, "types"]
      error <- runEvals[i,"errors"]
      mlevel <- runEvals[i, "mlevels"]
      model <- runEvals[i, "models"]
      runID <- with(runRequests[runIDs, ], which(types == type & errors == error & mlevels == mlevel & models == model))
      eval2.SDMs(i, runEval=runEvals[i, ], type=type, error = error, mlevel=mlevel, model=model, runID=runID, bsub=bres[runID])
    }
    
    clusterEvalQ(cl, rm(list=ls()))
    clusterEvalQ(cl, gc())
  }

  # Get data from disk
  beval <- list()
  ib <- 1
  for (i in xt) {
  	temp <- readRDS(file = get_temp_fname(runEvals[i, ], runEvalIDs[i]))
  	if (!inherits(temp, "try-error")) {
  		beval[[ib]] <- temp
  		ib <- ib + 1
  	}
  }

  
  #Identify runs
  goodRuns <- sapply(beval, FUN=function(l) !(inherits(l, c("try-error", "simpleError"))))
  if(any(!goodRuns)) print(beval[[which(!goodRuns)[1]]])
  if(all(!goodRuns)) stop(paste(Sys.time(), ": No evaluation successful"))
  print(paste(Sys.time(), ":", sum(goodRuns), "out of", length(goodRuns), "evaluations successful"))
  beval <- beval[goodRuns]
  temp <- t(sapply(beval, FUN=function(l) l$evalID))
  #evalIDs = index for which row of runEvals corresponds to elements of beval
  evalIDs <- na.exclude(match(apply(temp, 1, paste, collapse="_"), table=apply(runEvals, 1, paste, collapse="_")))
  
  save(beval, evalIDs, file=file.path(dir.sdm, filename.saveEvals))
  
  print(paste(Sys.time(), ": Evaluation done"))
}


## Summarize model evaluations
if(do.EvaluationSummary){
  print(paste(Sys.time(), ": Evaluation summary started"))
  
  if(!exists("beval") || !exists("evalIDs")){
    load(file.path(dir.sdm, filename.saveEvals))
  }
  
  #Fill evaluation arrays with values from bres
  evalA_SDMs <- array(NA, dim=c(length(types), length(errors), length(mlevels), length(sdm.models), length(eval.methods)+1, 2), dimnames=list(types,errors, names(mlevels), sdm.models, c(eval.methods, "Deviance"), c("mean", "sd")))
  evalA_Proj <- array(NA, dim=c(length(types), length(errors), length(mlevels), length(regions), length(sdm.models), length(eval.methods)+2, 2), dimnames=list(types,errors, names(mlevels), paste0("region", regions), sdm.models, c(eval.methods, "RMSE", "MAE"), c("mean", "sd")))
  evalA_ProjDiffs <- array(NA, dim=c(length(types), length(errors), length(mlevels), length(regions), length(sdm.models), length(eval.methods)+2, 2), dimnames=list(types,errors, names(mlevels), paste0("region", regions), sdm.models, c(eval.methods, "RMSE", "MAE"), c("mean", "sd")))
  
  for(i in evalIDs){
    evalID <- evalIDs[i]
    it <- which(runEvals[i, "types"] == types)
    ie <- which(runEvals[i, "errors"] == errors)
    il <- which(runEvals[i, "mlevels"] == names(mlevels))
    im <- which(runEvals[i, "models"] == sdm.models)
    
    evalA_SDMs[it, ie, il, im, , "mean"] <- as.numeric(c(t(beval[[evalID]]$Eval$mean[, "Testing.data"]), beval[[evalID]]$Deviance["mean"]))
    evalA_SDMs[it, ie, il, im, , "sd"] <- as.numeric(c(t(beval[[evalID]]$Eval$sd[, "Testing.data"]), beval[[evalID]]$Deviance["sd"]))
    for(ir in seq_along(regions)){
      evalA_Proj[it, ie, il, ir, im, , "mean"] <- as.numeric(c(t(beval[[evalID]]$Proj[[ir]]$Eval$mean[, "Testing.data"]), beval[[evalID]]$Proj[[ir]]$EvalProb$mean))
      evalA_Proj[it, ie, il, ir, im, , "sd"] <- as.numeric(c(t(beval[[evalID]]$Proj[[ir]]$Eval$sd[, "Testing.data"]), beval[[evalID]]$Proj[[ir]]$EvalProb$sd))
      evalA_ProjDiffs[it, ie, il, ir, im, , "mean"] <- as.numeric(c(t(beval[[evalID]]$Proj[[ir]]$EvalDiffToBase$mean[, "Testing.data"]), beval[[evalID]]$Proj[[ir]]$EvalProbDiffToBase$mean))
      evalA_ProjDiffs[it, ie, il, ir, im, , "sd"] <- as.numeric(c(t(beval[[evalID]]$Proj[[ir]]$EvalDiffToBase$sd[, "Testing.data"]), beval[[evalID]]$Proj[[ir]]$EvalProbDiffToBase$sd))
    }
    
  }
  
  #Reshape arrays into matrices
  evalT_SDMs <- acast(melt(evalA_SDMs), formula=Var1+Var2+Var3+Var4~Var5+Var6)
  evalT_Proj <- acast(melt(evalA_Proj), formula=Var1+Var2+Var3+Var4+Var5~Var6+Var7)
  evalT_ProjDiffs <- acast(melt(evalA_ProjDiffs), formula=Var1+Var2+Var3+Var4+Var5~Var6+Var7)
  
  write.csv(evalT_SDMs, file=file.path(dir.tables, "Table_EvaluationModels.csv"))
  write.csv(evalT_Proj, file=file.path(dir.tables, "Table_EvaluationProjections.csv"))
  write.csv(evalT_ProjDiffs, file=file.path(dir.tables, "Table_EvaluationDifferencesProjections.csv"))
  
  #Dave's figures
	try(source(file.path(path.functions, "ResponseExtrapolation_PlotPerformance_Bell_and_Schlaepfer.R")))
 
  print(paste(Sys.time(), ": Evaluation summary done"))
}


#does a serialized version work
if(FALSE){

	print(paste(Sys.time(), ": Figures started"))

	if (!exists("bres")) bres <- readRDS(file = file.path(dir.sdm, filename.saveSDMs))
	if (!exists("runIDs")) runIDs <- readRDS(file = file.path(dir.sdm, filename.saveRunIDs))

	for(i in 1:nrow(runEvals)){	
		type <- runEvals[i, "types"]
		error <- runEvals[i, "errors"]
		mlevel <- runEvals[i, "mlevels"]
		model <- runEvals[i, "models"]
		runID <- with(runRequests[runIDs, ], which(types == type & errors == error & mlevels == mlevel & models == model))
		bsub <- bres[runID]
	
		temp <- make2.figures(i,type,error,mlevel,model,runID,bsub)
	}

}

## Create figures
if(do.Figures){
	print(paste(Sys.time(), ": Figures started"))

	if (!exists("bres")) bres <- readRDS(file = file.path(dir.sdm, filename.saveSDMs))
	if (!exists("runIDs")) runIDs <- readRDS(file = file.path(dir.sdm, filename.saveRunIDs))


	## Read rasters
	grids <- list()
#	for(ir in regions) grids[[ir]] <- get_GeographicRaster_FromFile(ir)

	list.export <- c("libraries", "grids", "regions", "baseRegion","climData","probData","obsData","centerMeansData", 
		"make2.figures", "map_distributions", "plot_scatterPredvsTrueProbs", "plot_responseCurves2",
		"dir.sdm", "dir.in", "dir.maps", "dir.figs")
	if(identical(parallel_backend, "mpi")){
		exportObjects(c("work", list.export))
		mpi.bcast.cmd(lapply(libraries, FUN=require, character.only=TRUE))
		mpi.bcast.cmd(work())

		jobs <- 1:nrow(runEvals)
		res <- 0 #Result container
		workersN <- (mpi.comm.size() - 1)
		junk <- 0
		closed_slaves <- 0
		runs.completed <- 1

		while(closed_slaves < workersN) {
			complete <- mpi.recv.Robj(mpi.any.source(), mpi.any.tag())
			complete_info <- mpi.get.sourcetag()
			slave_id <- complete_info[1]
			tag <- complete_info[2]
			
			if(tag == 1) {# slave is ready for a task. Give it the next task, or tell it tasks are done if there are none.
				if(runs.completed <= length(jobs)){# Send a task, and then remove it from the task list
					dataForRun <- list(fun="figures",
										i=i <- jobs[runs.completed],
										type=type <- runEvals[i, "types"],
										error=error <- runEvals[i, "errors"],
										mlevel=mlevel <- runEvals[i, "mlevels"],
										model=model <- runEvals[i, "models"],
										runID=runID <- with(runRequests[runIDs, ], which(types == type & errors == error& mlevels == mlevel & models == model)),
										bsub=bres[runID])
										
					print(paste(Sys.time(), ": Figures for", dataForRun$i, "started"))
					mpi.send.Robj(dataForRun, slave_id, 1)
					runs.completed <- runs.completed + 1
				} else {
					mpi.send.Robj(junk, slave_id, 2)
				}
			} else if (tag == 2) { # The message contains results
				res <- res + complete$r
				print(paste(Sys.time(), ": Figures for", complete$i, "ended"))
			} else if (tag == 3) { # A slave has closed down.
				closed_slaves <- closed_slaves + 1
			}		
		}								

		mpi.bcast.cmd(rm(list=ls()))
		mpi.bcast.cmd(gc())
	}
	if(identical(parallel_backend, "parallel")){
		clusterExport(cl, c("bres", "runRequests", "runEvals", "runIDs", list.export)) #TODO: exporting large bres will fail
		clusterEvalQ(cl, lapply(libraries, FUN=require, character.only=TRUE))
	
		list.noexport <- (temp <- ls())[!(temp %in% list.export)]
		res <- foreach(i=1:nrow(runEvals), .combine="+", .errorhandling="remove", .noexport=list.noexport) %dopar% {
					type <- runEvals[i, "types"]
					error <- runEvals[i, "errors"],
					mlevel <- runEvals[i, "mlevels"]
					model <- runEvals[i, "models"]
					runID <- with(runRequests[runIDs, ], which(types == type & errors == error& mlevels == mlevel & models == model))
					make2.figures(i, type=type, error = error, mlevel=mlevel, model=model, runID=runID, bsub=bres[runID])
				}
	
		clusterEvalQ(cl, rm(list=ls()))
		clusterEvalQ(cl, gc())
	}

	#Identify runs
	if(res == nrow(runEvals)) warning(paste(Sys.time(), ": No figure set completely successful"))
	print(paste(Sys.time(), ":", nrow(runEvals) - res, "out of", nrow(runEvals), "figure sets successful"))

  #plot full map
    #make data, (Obs, Fit, XY, model, fun = mean, maxPred = 1000, figname)
   


plot.type<- 'SCT'
plot.error <- "binom"
plot.var <- 'LnP'

save.resp <- TRUE

if(save.resp){
  jpeg(file.path(dir.figs, paste0('ResponseCurves_',plot.type,'_',plot.error,'_',plot.var,'_','.jpeg')),width=6,height=6,units='in',res=300)
  par(mfrow=c(4,4),mar=c(1,1,0,0))
}

for(re in which(runEvals[,'types'] == plot.type & runEvals[,'errors'] == plot.error)[order(runEvals[runEvals[,'types'] == plot.type & runEvals[,'errors'] == plot.error,'models'],decreasing=TRUE)]){    

    ids <- which(apply(runRequests[,1:3],1,paste,collapse='.') == paste(runEvals[re,1:3],collapse='.'))
    ids <- cbind(ids,runRequests[ids,])
    for(ir in seq_along(regions)){

      Ftmp <- Otmp <- matrix(NA, nrow=length(bres[[ids[1,1]]]$Proj[[ir]]$Proj$pred), ncol=nrow(ids))
      
      Rtmp <- list()
      
      if(ir == 1){
        coordsXY <- climData[(climData[, "region"] == regions[ir]), c("x", "y")]
  
        for(rr in 1:nrow(ids)){
          Ftmp[,rr] <- bres[[ids[rr,1]]]$Proj[[ir]]$Proj$pred
          Otmp[,rr] <- obsData[[paste(ids[rr,'types'],ids[rr,'errors'],sep="_")]][(climData[, "region"] == regions[ir]), ids[rr, 'realizations']]
          Rtmp[[rr]] <- bres[[ids[rr,1]]]$ResponseCurvePreds[[ir]]
          
        }

        Fit <- Ftmp
        Obs <- Otmp
        
      }
        if(ir > 1){
          coordsXY <- rbind(coordsXY,
                            climData[(climData[, "region"] == regions[ir]), c("x", "y")])
          
          for(rr in 1:nrow(ids)){
            Ftmp[, rr] <- bres[[ids[rr,1]]]$Proj[[ir]]$Proj$pred
            Otmp[, rr] <- obsData[[paste(ids[rr,'types'],ids[rr,'errors'],sep="_")]][(climData[, "region"] == regions[ir]), ids[rr, 'realizations']]
            Rtmp[[rr]] <- bres[[ids[rr,1]]]$ResponseCurvePreds[[ir]]
          }
        
          Fit <- rbind(Fit,Ftmp)
          Obs <- rbind(Obs, Otmp)
          
        }
            
      plot_responseCurves_CurvesOnly(newdata = climData[climData[,'region'] == regions[ir],c('LnP','LnPScaled','MinT','MinTScaled')],
                           Prob = probData[[paste(runEvals[re,'types'],runEvals[re,'errors'],sep="_")]][climData[,'region'] == regions[ir]],
                           Obs = Otmp, Fit = Ftmp,respCurvePreds = Rtmp,
                           model = runEvals[re,'models'],
                           maxPred=1000,
                           centerMeans = centerMeansData,
                           env = plot.var)
      print(paste(ir,ids[1,'models'],ids[1,'mlevels']))

      }



map_distributions(Obs = Obs,Fit = Fit,XY = coordsXY,model = "",
                  figname = file.path(dir.maps, paste0("Map_TrueVsPredicted_", 
                                                       ids[1,'types'], "_", 
                                                       ids[1,'errors'], "_", 
                                                       ids[1,'mlevels'], "_", 
                                                       ids[1,'models'],'FULL', ".jpg")))

map_distributions(Obs = Obs,Fit = Fit,XY = coordsXY,model = "",fun=sd,
                  figname = file.path(dir.maps, paste0("Map_UncertaintyVsPredicted_", 
                                                       ids[1,'types'], "_", 
                                                       ids[1,'errors'], "_", 
                                                       ids[1,'mlevels'], "_", 
                                                       ids[1,'models'],'FULL', ".jpg")))


}
dev.off()


  jpeg(file.path(dir.maps, paste0('maps_',plot.type,'_',plot.error,'_','.jpeg')),width=5,height=6.5,units='in',res=600)

plt.mat <- rbind(c(0.06,0.36,0.75,0.95), #a
                 c(0.38,0.68,0.75,0.95), #e
                 c(0.70,1.00,0.75,0.95), #i
                 c(0.06,0.36,0.53,0.73), #b
                 c(0.38,0.68,0.53,0.73), #f
                 c(0.70,1.00,0.53,0.73), #j
                 c(0.06,0.36,0.31,0.51), #c
                 c(0.38,0.68,0.31,0.51), #g
                 c(0.70,1.00,0.31,0.51), #k
                 c(0.06,0.36,0.09,0.29), #d
                 c(0.38,0.68,0.09,0.29), #h
                 c(0.70,1.00,0.09,0.29), #l
                 c(0.09,0.33,0.01,0.02), #pred legend
                 c(0.41,0.65,0.01,0.02), #change legend
                 c(0.73,0.97,0.01,0.02))#sd legend

plot.new()
             
let <- c("(a)","(e)","(i)",
         "(b)","(f)","(j)",
         "(c)","(g)","(k)",
         "(d)","(h)","(l)")
mod <- c("GLM w/o Inter.", "GLM w Inter.", "GAM w/o Inter.", "GAM w Inter.")
runmods <- which(runEvals[,'types'] == plot.type & runEvals[,'errors'] == plot.error)[order(runEvals[runEvals[,'types'] == plot.type & runEvals[,'errors'] == plot.error,'models'],decreasing=TRUE)]

for(re in runmods){    
  
  ids <- which(apply(runRequests[,1:3],1,paste,collapse='.') == paste(runEvals[re,1:3],collapse='.'))
  ids <- cbind(ids,runRequests[ids,])
  for(ir in seq_along(regions)){
    
    Ftmp <- Otmp <- matrix(NA, nrow=length(bres[[ids[1,1]]]$Proj[[ir]]$Proj$pred), ncol=nrow(ids))
    
    Rtmp <- list()
    
    if(ir == 1){
      coordsXY <- climData[(climData[, "region"] == regions[ir]), c("x", "y")]
      
      for(rr in 1:nrow(ids)){
        Ftmp[,rr] <- bres[[ids[rr,1]]]$Proj[[ir]]$Proj$pred
        Otmp[,rr] <- obsData[[paste(ids[rr,'types'],id[rr,"errors"],sep="_")]][(climData[, "region"] == regions[ir]), ids[rr, 'realizations']]
        Rtmp[[rr]] <- bres[[ids[rr,1]]]$ResponseCurvePreds[[ir]]
        
      }
      
      Fit <- Ftmp
      Obs <- Otmp
      
    }
    if(ir > 1){
      coordsXY <- rbind(coordsXY,
                        climData[(climData[, "region"] == regions[ir]), c("x", "y")])
      
      for(rr in 1:nrow(ids)){
        Ftmp[, rr] <- bres[[ids[rr,1]]]$Proj[[ir]]$Proj$pred
        Otmp[, rr] <- obsData[[paste(ids[rr,'types'],id[rr,"errors"],sep="_")]][(climData[, "region"] == regions[ir]), ids[rr, 'realizations']]
        Rtmp[[rr]] <- bres[[ids[rr,1]]]$ResponseCurvePreds[[ir]]
      }
      
      Fit <- rbind(Fit,Ftmp)
      Obs <- rbind(Obs, Otmp)
      
    }
    
      print(paste(ir,ids[1,'models'],ids[1,'mlevels']))
    
  }
  
  #Obs, Fit, XY, model, fun = mean, maxPred = 1000, figname,save=TRUE)
  maxPred <- 1000
  
  zmax <- max(Fit, maxPred)/1000
  
  pred <- apply(Fit/maxPred, MARGIN=1, FUN=function(x) mean(x))
  pred.diff <- apply(Fit/maxPred - Obs, MARGIN=1, FUN=mean)
  pred.sd <- apply(Fit/maxPred, MARGIN=1, FUN=sd)
    
  
  cols.pred  <- rev(terrain.colors(n=maxPred))
  cols.diff <- c(cm.colors(n=maxPred))
  cols.sd   <- rev(c(heat.colors(n=maxPred),"#F2F2F2FF"))
  
  pred.seq  <- seq(0,1,length=maxPred)
  diff.seq <- c(-10,seq(-1,1,length=maxPred))
  sd.seq   <- c(0,seq(0.01,.5,length=maxPred))
  
  par(plt=plt.mat[(which(runmods == re) - 1) * 3 + 1,],new = TRUE)

  plot(coordsXY, pch=15, cex=0.7, col="black", asp=1, xlab="", ylab="", axes=FALSE)
  points(coordsXY, pch=15, cex=0.43, col=cols.pred[findInterval(pred,pred.seq)])
    legend("bottomleft",legend=let[(which(runmods == re) - 1) * 3 + 1],bty="n")

    mtext(mod[which(runmods == re)],side=2,las=3,line=-.25)
    if(re == runmods[1]) mtext("Prediction",side=3,line=.2)

  par(plt=plt.mat[(which(runmods == re) - 1) * 3 + 2,],new = TRUE)

  plot(coordsXY, pch=15, cex=0.7, col="black", asp=1, xlab="", ylab="", axes=FALSE)
  points(coordsXY, pch=15, cex=0.43, col=cols.diff[findInterval(pred.diff,diff.seq)], asp=1, xlab="", ylab="", axes=FALSE)
    legend("bottomleft",legend=let[(which(runmods == re) - 1) * 3 + 2],bty="n")
    if(re == runmods[1]) mtext("Bias",side=3,line=.2)

  par(plt=plt.mat[(which(runmods == re) - 1) * 3 + 3,],new = TRUE)

  plot(coordsXY, pch=15, cex=0.7, col="black", asp=1, xlab="", ylab="", axes=FALSE)
  points(coordsXY, pch=15, cex=0.43, col=cols.sd[findInterval(pred.sd,sd.seq)], asp=1, xlab="", ylab="", axes=FALSE)
    legend("bottomleft",legend=let[(which(runmods == re) - 1) * 3 + 3],bty="n")
    if(re == runmods[1]) mtext("SD",side=3,line=.2)

  
 # axis(side=1)
#  axis(side=2)
#  points(XY, pch=16, cex=ifelse(obs > 0, 0.06 + 0.25*obs, 0), col = adjustcolor("black", alpha.f = 0.4))
#  mtext(text=model)
  
  #Legend
#  if(save) par(mar=c(2,0,0.1,0))
#  plot(x=seq(from=0, to=zmax, length=maxPred), y=rep(0, maxPred), xlim=c(0, zmax), col=cols, xlab="Suitability score", axes=FALSE)
#  axis(side=1, pos=-0.1, cex=0.6)
  
 
  
}

  par(plt=plt.mat[13,],new = TRUE)
  image(pred.seq,1:2,matrix(pred.seq[-1],ncol=1),col=cols.pred[-length(cols.pred)],
        ylab="",xlab="",axes=FALSE,frame.plot=TRUE)
  box()
  mtext(c("0.0","0.5","1.0"), side=3, at = c(0,.5,1),line=.2,cex=.8)

  par(plt=plt.mat[14,],new = TRUE)
  image(diff.seq[-1],1:2,matrix(diff.seq[-1],ncol=1),col=cols.diff[-length(cols.diff)],
        ylab="",xlab="",axes=FALSE,frame.plot=TRUE)
  box()
  mtext(c("-1.0","0.0","1.0"), side=3, at = c(-1,0,1),line=.2,cex=.8)

  par(plt=plt.mat[15,],new = TRUE)
  image(sd.seq,1:2,matrix(sd.seq[-1],ncol=1),col=cols.sd[-length(cols.sd)],
        ylab="",xlab="",axes=FALSE,frame.plot=TRUE)
  box()
  mtext(c("0.0","0.5","1.0"), side=3, at = c(0,.5,1),line=.2,cex=.8)


dev.off()

rm(Ftmp,Otmp)

  
#finish
  print(paste(Sys.time(), ": Figures done"))
}

## Clean up parallel backends
print(paste(Sys.time(), ": Clean-up"))

if(identical(parallel_backend, "mpi")){
	#mpi.close.Rslaves(dellog=FALSE)
	mpi.exit()
}
