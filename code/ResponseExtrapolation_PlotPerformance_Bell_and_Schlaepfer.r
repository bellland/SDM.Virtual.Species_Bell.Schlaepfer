###############################################################################
#
# Bell, D. M., and D. R. Schlaepfer. Impacts of the data-generating processes and species distribution model complexity on ecological fidelity and global change predictions.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
###############################################################################

#---Helper functions
	nice_eval_labels <- function(x) sub("KAPPA", "expression(kappa)", sub("ROC", "AUC", x))

	add_yaxis <- function(limited_to01 = TRUE) {
		temp <- par("usr")[3:4]
		if (limited_to01) {
			temp[1] <- max(0, temp[1])
			temp[2] <- min(1, temp[2])
		}
		fact <- 10 ^ (1 + ceiling(log(1 + ceiling(1 / diff(temp)), base = 10)))
		at <- temp
		at[1] <- if (temp[1] >= 0) max(0, ceiling(fact * temp[1] * 1.04) / fact) else floor(fact * temp[1] * 0.96) / fact
		at[2] <- if (temp[2] >= 0) min(1, floor(fact * temp[2]) / fact) else ceiling(fact * temp[2] * 0.96) / fact
		at <- c(at[1], round(mean(at), floor(log(fact, base = 10))), at[2])
		if (!limited_to01 && at[1] < 0 && at[3] > 0) at <- sort(c(0, at))
		
		axis(side = 2, at = at, las = 2, cex = cex)
	}

	draw.panel <- function(dat, meanCol, sdCol, pchs, cols, subset = 1:nrow(dat), divCol = colnames(dat)[1], ylab = "", cex = 1){
		ynegs <- dat[, meanCol] - dat[, sdCol]
		ypos <- dat[, meanCol] + dat[, sdCol]

		N <- if (is.logical(subset)) sum(subset) else length(subset)
		exx <- 0.5
		dividers <- seq(from = exx + 0, to = exx + N, length.out = length(unique(dat[subset, divCol])) + 1)
	
		ylab <- if (nchar(ylab) > 1 && grepl("expression", ylab)) eval(parse(text = ylab)) else ylab
	
		plot(dat[subset, meanCol], axes=FALSE, frame.plot=TRUE, xlab = "", ylab = "",
				xlim = c(dividers[1] - exx, tail(dividers, n = 1) + exx), ylim = c(min(ynegs), max(ypos)),
				pch = pchs, col = cols, cex = cex)
		arrows(1:nrow(dat), ynegs[subset], 1:nrow(dat), ypos[subset],
				angle = 90, code = 3, length = .05, col = cols)
	
		mtext(side = 2, line = 2, text = ylab, xpd = NA)
	
		abline(v = dividers[-c(1, length(dividers))], lty = 2)
	
		invisible(dividers)
	}

	get_pchs <- function(data, col1 = "Model", col2 = "Int") {
		temp <- expand.grid(levels(data[, col2]), levels(data[, col1]))[, 2:1]
		ptemp <- list(c(0, 15), c(1, 16), c(2, 17), c(5, 18), c(8, 11))
		pch_ids <- data.frame(temp, pch = unlist(ptemp)[1:nrow(temp)])
		pchs <- rep(NA, nrow(data))
		for (i in 1:nrow(pch_ids)) {
			pchs[as.character(data[, col1]) == pch_ids[i, 1] & as.character(data[, col2]) == pch_ids[i, 2]] <- pch_ids[i, "pch"]
		}
		
		list(pch_ids = pch_ids, pchs = pchs)
	}
  
	get_cols <- function(data, col1 = "Error") {
		temp <- expand.grid(levels(data[, col1]))
		ptemp <- c("black", "darkred", "gray40", "red")
		col_ids <- data.frame(temp, col = unlist(ptemp)[1:nrow(temp)], stringsAsFactors = FALSE)
		cols <- rep(NA, nrow(data))
		for (i in 1:nrow(col_ids)) {
			cols[as.character(data[, col1]) == col_ids[i, 1]] <- col_ids[i, "col"]
		}
		
		list(col_ids = col_ids, cols = cols)
	}

	convert_data <- function(data, include_region = FALSE, include_error = TRUE) {
		n_header <- (if (include_region) 5 else 4) - (if (include_error) 0 else 1)
		col_header <- c('Source', if (include_error) 'Error' else NULL, 'Int', if (include_region) 'Region' else NULL, 'Model')
		temp <- c('Source', 'Model', 'Int')
		sort_header <- c(temp, col_header[!(col_header %in% temp)])
		
		data <- cbind(t(matrix(unlist(strsplit(as.character(data[, 1]), '_')), nrow = n_header)), data[, -1])
		colnames(data) <- c(col_header, colnames(data)[-(1:n_header)])
		data[do.call(order, data[, sort_header]), ]
	}


#---Loop through aggregation levels and draw figures
	collapse.levels <- c("full", "collerr")
	
	for (colllev in collapse.levels) {
		include_error <- switch(EXPR = colllev, full = TRUE, collerr = FALSE)
	
		#---Load data  
		model.perf <- read.csv(file.path(dir.tables, paste0('Table_EvaluationModels_', colllev, '.csv')), header = TRUE)
		model.perf <- convert_data(model.perf, include_region = FALSE, include_error = include_error)
		
		proj.perf <- read.csv(file.path(dir.tables, paste0('Table_EvaluationProjections_', colllev, '.csv')), header = TRUE)
		proj.perf <- convert_data(proj.perf, include_region = TRUE, include_error = include_error)
  
		projDiff.perf <- read.csv(file.path(dir.tables, paste0('Table_EvaluationDifferencesProjections_', colllev, '.csv')), header = TRUE)
		projDiff.perf <- convert_data(projDiff.perf, include_region = TRUE, include_error = include_error)


		#---Prepare plotting information
		n_units <- nrow(unique(model.perf[, seq_len(if (include_error) 4 else 3)]))
		pchs <- get_pchs(model.perf, col1 = "Model", col2 = "Int")
		cols <- get_cols(model.perf, col1 = if (include_error) "Error" else "Int")
		
		cex <- 1
		h.panel <- 2.6
		w.panel <- 0.075 * n_units


		#---FIGURE. MODEL PERFORMANCE IN TRAINING REGION
		temp <- eval_disc.methods[c(1, 3, 2)]
		meanCols <- paste(temp, "mean", sep = "_")
		sdCols <- paste(temp, "sd", sep = "_")
		ylabels <- nice_eval_labels(temp)
	
		npanelsX <- 1; npanelsY <- length(ylabels)
		h.edge_lo <- 0.0; h.edge_up <- 0.2
		w.edge_left <- 0.55; w.edge_right <- 0.0
		panels <- matrix(0, nrow = 1 + npanelsY + 1, ncol = 1 + npanelsX + 1, byrow=FALSE)
		panels[-c(1, nrow(panels)), -c(1, ncol(panels))] <- 1:(npanelsX * npanelsY)

		png(height = h.edge_lo + h.panel * npanelsY + h.edge_up,
			width = w.edge_left + w.panel * npanelsX + w.edge_right, units = "in",
			res = 600, file = file.path(dir.figs, paste0('Model.Performance_', colllev, '.png')))

		layout(panels,
				heights = c(h.edge_up, rep(h.panel, times = npanelsY), h.edge_lo),
				widths = c(w.edge_left, rep(w.panel, times = npanelsX), w.edge_right))
		op_old <- par(mgp = c(2, 0.25, 0), mar = c(0.2, 0.2, 0.2, 0.2), tcl = 0.3, cex = cex, xaxs = "i")

		for (i in seq_along(meanCols)) {
			divs <- draw.panel(dat = model.perf, meanCol = meanCols[i], sdCol = sdCols[i],
								pchs = pchs$pchs, cols = cols$cols,
								subset = 1:nrow(model.perf), divCol = "Source", ylab = ylabels[i], cex = cex)
		
			#at <- axisTicks(usr = par("usr")[3:4], log = FALSE, nint = 2)
			add_yaxis()
		
			if (i == 1) {
				mtext(unique(model.perf[, "Source"]), side = 3, at = zoo::rollmean(divs, 2), cex = cex)
		
				legend("bottomright", legend = apply(pchs$pch_ids[, 1:2], 1, paste, collapse = "; "),
								ncol = 3, pch = pchs$pch_ids[, "pch"], bg = "white", cex = 0.75 * cex)
		
				legend("bottomright", inset = c(0.45, 0), legend = cols$col_ids[, 1],
								ncol = 2, fill = cols$col_ids[, "col"], bg = "white", cex = 0.75 * cex)
			}
		}
	
		par(op_old)
		dev.off()
  

  
		#---FIGURE 3. REGIONAL MODEL PERFORMANCE; ABSOLUTE FOR TRAINING REGION; RELATIVE FOR PROJECTION REGIONS
		temp <- c(eval_disc.methods[c(1, 3, 2)], eval_cont.methods)
		meanCols <- paste(temp, "mean", sep = "_")
		sdCols <- paste(temp, "sd", sep = "_")
		ylabels <- nice_eval_labels(temp)
		ylabelsDiff <- sapply(ylabels, function(x) {
							if (grepl("expression", x)) {
								temp <- strsplit(strsplit(x, ")", fixed = TRUE)[[1]], "(", fixed = TRUE)[[1]]
								paste0(temp[1], "(paste(", temp[2], ", ' Difference'))")
							} else {
								paste(x, "Difference")
							}
						})
	
		npanelsX <- length(regions); npanelsY <- length(ylabels)
		h.edge_lo <- 0.5; h.edge_up <- 0.2
		w.edge_left <- 0.55; w.edge_right <- 0.0
		panels <- matrix(0, nrow = 1 + npanelsY + 1, ncol = 2 + npanelsX + 1, byrow=FALSE)
		panels[-c(1, nrow(panels)), -c(1, 3, ncol(panels))] <- 1:(npanelsX * npanelsY)

		png(height = h.edge_lo + h.panel * npanelsY + h.edge_up,
			width = 2 * w.edge_left + w.panel * npanelsX + w.edge_right, units = "in",
			res = 600, file = file.path(dir.figs, paste0('Fig3_Performance_Training_Transferability_', colllev, '.png')))

		layout(panels,
				heights = c(h.edge_up, rep(h.panel, times = npanelsY), h.edge_lo),
				widths = c(w.edge_left, w.panel, w.edge_left, rep(w.panel, times = npanelsX - 1), w.edge_right))
		op_old <- par(mgp = c(2, 0.25, 0), mar = c(0.2, 0.2, 0.2, 0.2),
						tcl = 0.3, cex = cex, xaxs = "i")
	
		for (j in seq_along(reg.sort)) {
			rsubset <- proj.perf[, "Region"] == paste0("region", reg.sort[j])
		
			for (i in seq_along(meanCols)) {
				dats <- if (reg.sort[j] == baseRegion) proj.perf else projDiff.perf
				ylab <- if (j == 1) ylabels[i] else if (j == 2) ylabelsDiff[i] else ""
	
				divs <- draw.panel(dat = dats, meanCol = meanCols[i], sdCol = sdCols[i],
									pchs = pchs$pchs, cols = cols$cols,
									subset = rsubset, divCol = "Source", ylab = ylab, cex = cex)
			
				if (j %in% 1:2) add_yaxis(j == 1)
				
				if (j > 1) abline(h = 0, lty = 2, col = "gray")
		
				if (i == 1) {
					mtext(unique(proj.perf[, "Source"]), side = 3, at = zoo::rollmean(divs, 2), cex = cex)
				
					if (j == 1) {
						legend("bottomright", legend = apply(pchs$pch_ids[, 1:2], 1, paste, collapse = "; "),
										ncol = 3, pch = pchs$pch_ids[, "pch"], bg = "white", cex = 0.75 * cex)
		
						legend("bottomright", inset = c(0.45, 0), legend = cols$col_ids[, 1],
										ncol = 2, fill = cols$col_ids[, "col"], bg = "white", cex = 0.75 * cex)
					}
				}
			
				if (i == length(meanCols)) {
					mtext(names(reg.sort)[j], side = 1, at = mean(divs), cex = cex)
				
					mtemp <- if (reg.sort[j] == baseRegion) {
									"Training Region"
								} else if (j == 1 + round((length(reg.sort) - 1) / 2)) {
									"Transferability Regions"
								} else ""
					mtext(mtemp, side = 1, at = mean(divs), line = 1, cex = cex)
				}
			}
		}
	
		par(op_old)
		dev.off()
	} 

warning("Not all figures of 'ResponseExtrapolation_PlotPerformance_Bell_and_Schlaepfer' re-implemented")
if (FALSE) {
	#####################
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
}  
