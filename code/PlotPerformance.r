#clear workspace

 # rm(list=ls())

#load table

  setwd(dir.tables)

  model.perf <- read.csv('Table_EvaluationModels.csv',header=TRUE)
  proj.perf  <- read.csv('Table_EvaluationProjections.csv',header=TRUE)
  projDiff.perf  <- read.csv('Table_EvaluationDifferencesProjections.csv',header=TRUE)
  
  
  model.perf <- cbind(t(matrix(unlist(strsplit(as.character(model.perf[,1]),'_')),nrow=4)),model.perf[,-1])
    colnames(model.perf) <- c('Source','Error','Int','Model',colnames(model.perf)[-(1:4)])
  
  proj.perf <- cbind(t(matrix(unlist(strsplit(as.character(proj.perf[,1]),'_')),nrow=5)),proj.perf[,-1])
    colnames(proj.perf) <- c('Source','Error','Int','Region','Model',colnames(proj.perf)[-(1:5)])
  
  projDiff.perf <- cbind(t(matrix(unlist(strsplit(as.character(projDiff.perf[,1]),'_')),nrow=5)),projDiff.perf[,-1])
    colnames(projDiff.perf) <- c('Source','Error','Int','Region','Model',colnames(projDiff.perf)[-(1:5)])
  
  pchs <- rep(1,nrow(model.perf))
    pchs[as.character(model.perf[,"Model"]) == 'GAM' & as.character(model.perf[,"Int"]) == 'wInt'] <- 16
    pchs[as.character(model.perf[,"Model"]) == 'GAM' & as.character(model.perf[,"Int"]) == 'woInt'] <- 15
    pchs[as.character(model.perf[,"Model"]) == 'GLM' & as.character(model.perf[,"Int"]) == 'woInt'] <- 0
  cols <- rep("black",nrow(model.perf))
    cols[as.character(model.perf[,"Error"]) == 'binom+res'] <- "darkred"
    cols[as.character(model.perf[,"Error"]) == 'spatial'] <- "gray40"
    cols[as.character(model.perf[,"Error"]) == 'spatial+res'] <- "red"
  
  jpeg('Model.Performance.jpg',width=6,height=8,units = "in", res = 1000)
  par(mfrow=c(4,1),mar=c(2,5,2,1))
  
  plot(model.perf[,'TSS_mean'],axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(min(model.perf[,'TSS_mean'] - model.perf[,'TSS_sd']),
              max(model.perf[,'TSS_mean'] + model.perf[,'TSS_sd'])),
       pch = pchs,
       col = cols)
  arrows(1:nrow(model.perf),model.perf[,'TSS_mean'] - model.perf[,'TSS_sd'],
         1:nrow(model.perf),model.perf[,'TSS_mean'] + model.perf[,'TSS_sd'],
         angle = 90, code = 3, length = .05,
         col = cols)
  abline(v=c(16.5,32.5,48.5),lty=2)
  
  mtext(unique(model.perf[,1]),side=3,at=c(8.5,24.5,40.5,56.5))
  legend('bottomright',
         legend=c('GLM; No Interaction',
                  'GAM; No Interaction',
                  'GLM; Interaction',
                  'GAM; Interaction'),
         pch = c(0,15,1,16),
         bg='white',
         cex=.8)
  mtext(seq(.85,.95,by=.05),side=2,at=seq(.85,.95,by=.05),las=2)
  mtext('Mean TSS (+/- sd)',side=2,line=3.0)
  
  plot(model.perf[,'KAPPA_mean'],axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(min(model.perf[,'KAPPA_mean'] - model.perf[,'KAPPA_sd']),
              max(model.perf[,'KAPPA_mean'] + model.perf[,'KAPPA_sd'])),
       pch = pchs,
       col = cols)
  arrows(1:nrow(model.perf),model.perf[,'KAPPA_mean'] - model.perf[,'KAPPA_sd'],
         1:nrow(model.perf),model.perf[,'KAPPA_mean'] + model.perf[,'KAPPA_sd'],
         angle = 90, code = 3, length = .05,
         col = cols)
  abline(v=c(16.5,32.5,48.5),lty=2)
  
  mtext(seq(.75,.8,by=.05),side=2,at=seq(.75,.8,by=.05),las=2)
  mtext('Mean KAPPA (+/- sd)',side=2,line=3.0)

  plot(model.perf[,'ROC_mean'],axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(min(model.perf[,'ROC_mean'] - model.perf[,'ROC_sd']),
              max(model.perf[,'ROC_mean'] + model.perf[,'ROC_sd'])),
       pch = pchs,
       col = cols)
  arrows(1:nrow(model.perf),model.perf[,'ROC_mean'] - model.perf[,'ROC_sd'],
         1:nrow(model.perf),model.perf[,'ROC_mean'] + model.perf[,'ROC_sd'],
         angle = 90, code = 3, length = .05,
         col = cols)
  abline(v=c(16.5,32.5,48.5),lty=2)
  
  mtext(seq(.95,1.0,by=.01),side=2,at=seq(.95,1.0,by=.01),las=2)
  mtext('Mean ROC (+/- sd)',side=2,line=3.0)

  legend('bottomright',
         legend=c("no error","error","autocovariation","error + autocovariation"),
         fill = c("black","darkred","gray40","red"),
         bg='white',
         cex=.8)
  
  plot(model.perf[,'Deviance_mean'],axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(min(model.perf[,'Deviance_mean'] - model.perf[,'Deviance_sd']),
              max(model.perf[,'Deviance_mean'] + model.perf[,'Deviance_sd'])),
       pch = pchs,
       col = cols)
  arrows(1:nrow(model.perf),model.perf[,'Deviance_mean'] - model.perf[,'Deviance_sd'],
         1:nrow(model.perf),model.perf[,'Deviance_mean'] + model.perf[,'Deviance_sd'],
         angle = 90, code = 3, length = .05,
         col = cols)
  abline(v=c(16.5,32.5,48.5),lty=2)
  
  #mtext(seq(.95,1.0,by=.01),side=2,at=seq(.95,1.0,by=.01),las=2)
  axis(side=2)
  mtext('Mean Deviance (+/- sd)',side=2,line=3.0)
  
  dev.off()
  



  
draw.panel <- function(dat, meanCol, sdCol, subset){

	plot(dat[subset, meanCol],axes=FALSE,frame.plot=TRUE,xlab="",ylab="",
       ylim=c(min(dat[, meanCol] - dat[, sdCol]),
              max(dat[, meanCol] + dat[, sdCol])),
       pch = pchs,
       col = cols)
	arrows(1:nrow(dat),dat[subset, meanCol] - dat[subset, sdCol],
         1:nrow(dat),dat[subset, meanCol] + dat[subset, sdCol],
         angle = 90, code = 3, length = .05,
         col = cols)
	abline(v=c(16.5,32.5,48.5),lty=2)
}

meanCols <- c("TSS_mean", "KAPPA_mean", "ROC_mean", "RMSE_mean", "MAE_mean")
sdCols <- gsub("_mean", "_sd", meanCols)

postscript('Projection.Performance.eps',width=18,height=9,onefile=FALSE,horizontal=FALSE)
par(mfrow=c(length(meanCols),length(regions)),mar=c(2,5,2,1))

for(ic in seq_along(meanCols)){
	for(ir in seq_along(regions)){
		draw.panel(dat=proj.perf, meanCol=meanCols[ic], sdCol=sdCols[ic],
					subset=proj.perf[, "Region"] == paste0("region", ir))
		axis(side=2, at=axTicks(2), labels=if(ir == 1) axTicks(2) else FALSE, cex=2, las=2)
		
		if(ir == 1) mtext(paste0("Mean ", sub("_mean", "", meanCols[ic]), " (+/- sd)"),side=2,line=3.0)
		if(ic == 1) mtext(unique(proj.perf[,"Source"]),side=3,at=c(8.5,24.5,40.5,56.5))
		if(ic == length(meanCols)) mtext(paste0("Region", ir), side=1, line=0.5)
		if(ic == 1 && ir == 1){
			legend('bottomright',
				legend=c('GLM; No Interaction',
						'GAM; No Interaction',
						'GLM; Interaction',
						'GAM; Interaction'),
				pch = c(0,15,1,16),
				bg='white',
				cex=.8)
		}
	}
}	

dev.off()
  


jpeg('ProjectionDifferences.Performance.jpg',width=18,height=9,units="in",res=1000)
par(mfrow=c(length(meanCols),length(regions)-1),mar=c(2,5,2,1))

for(ic in seq_along(meanCols)){
	for(ir in seq_along(regions)[-baseRegion]){
		draw.panel(dat=projDiff.perf, meanCol=meanCols[ic], sdCol=sdCols[ic],
					subset=projDiff.perf[, "Region"] == paste0("region", ir))
		abline(a=0, b=0)
		axis(side=2, at=axTicks(2), labels=if(ir == 1) axTicks(2) else FALSE, cex=2, las=2)
		
		if(ir == 1) mtext(paste0("Mean ", sub("_mean", "", meanCols[ic]), " (+/- sd)"),side=2,line=3.0)
		if(ic == 1) mtext(unique(projDiff.perf[,"Source"]),side=3,at=c(8.5,24.5,40.5,56.5))
		if(ic == length(meanCols)) mtext(paste0("Region", ir), side=1, line=0.5)
		if(ic == 1 && ir == 1){
			legend('bottomright',
				legend=c('GLM; No Interaction',
						'GAM; No Interaction',
						'GLM; Interaction',
						'GAM; Interaction'),
				pch = c(0,15,1,16),
				bg='white',
				cex=.8)
		}
		if(ic == 1 && ir == 4){
		  legend('bottomright',
		         legend=c("no error","error","autocovariation","error + autocovariation"),
		         fill = c("black","darkred","gray40","red"),
		         bg='white',
		         cex=.8)
		}
	}
}	

dev.off()
  
