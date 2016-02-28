our.response.plot2 <- function(modelObj, modelName, Data, orig.variables, scaled.variables=NULL, centerMeans, data_species, fixed.var.metric = 'mean'){
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
		data_species <- rep(1,nrow(Data))
	} else {
		data_species[data_species!=1 | is.na(data_species)] <- 0
	}


	# 2. build function outputs
	factor_id <- which(sapply(Data,is.factor))
	list.out <- list()

	# Create a ranged data table
	ref_table <- Data[1,,drop=F]
	rownames(ref_table) <- NULL

	for(i in 1:ncol(Data)){
		temp <- Data[data_species==1,i]
		if(is.numeric(Data[,i])){
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

	for(vari in show.variables){
		# creating Tmp data
		if(is.factor(Data[,vari])){
			pts.tmp <- as.factor(levels(Data[,vari]))
		} else {
			pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)
			temp <- grep(iorig <- sub("Scaled", "", vari), names(centerMeans))
			pts.tmp.orig <- if(isScaled) pts.tmp + centerMeans[temp] else NULL
		}

		Data.r.tmp <- eval(parse(text=paste("cbind(",vari,"=pts.tmp,ref_table[,-which(colnames(ref_table)==vari),drop=F])",sep="")))
		Data.r.tmp <- Data.r.tmp[,colnames(ref_table),drop=F]
		if(length(factor_id)){
			for(f in factor_id){
				Data.r.tmp[,f] <- factor(as.character(Data.r.tmp[,f]), levels=levels(Data[,f]))
			}
		}
		if(!is.null(pts.tmp.orig)) Data.r.tmp[, iorig] <- pts.tmp.orig

		# 2. make projections 
		proj.tmp <- predict(object=modelObj, newdata=Data.r.tmp, type='response', se.fit=FALSE)

		# 5. Storing results
		if(length(list.out[[vari]]) == 0){ #init
			eval(parse(text=paste("list.out[['",vari,"']] <- data.frame(",vari,"=pts.tmp, ",modelName,"=proj.tmp)",sep="")))
		} else {
			eval(parse(text=paste("list.out[['",vari,"']] <- cbind(list.out[['",vari,"']],",modelName,"=proj.tmp)",sep="")))
		}

	}

	invisible(list.out)

}
