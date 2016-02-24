library(raster)

dir_grid <- "~/Dropbox/Work_Stuff/2_Research/200907_UofWyoming_PostDoc/Product21_SDM_AsymmetryResponseCurve/3_Simulations/inst/extdata/rasters"

regions <- 1:4

# Read in mask of each region
r <- vector(mode = "list", length = length(regions))
for (i in regions) r[[i]] <- raster(file.path(dir_grid, "old", paste0("extent.raster.r", i, ".grd")))


# Set region value
for (i in regions) r[[i]] <- calc(r[[i]], function(x) ifelse(is.na(x), NA, i))

#---Save new rasters with region values (and datatype = "INT1S" for 8-times smaller files)
for (i in regions) writeRaster(r[[i]], filename = file.path(dir_grid, "res_1to32_degree", paste0("extent.raster.r", i, ".grd")), datatype = "INT1S")


#---Merge/mosaic
compare_raster <- sapply(r[-1], function(x) all.equal(r[[1]], x, extent = FALSE, rowcol = FALSE, orig = TRUE, res = TRUE))

if (all(compare_raster)) {
	rall <- merge(r[[1]], r[[2]], r[[3]], r[[4]])
} else {
	print(sapply(r, origin))
#             [,1]          [,2]         [,3]         [,4]
#[1,] -0.002309340 -1.341035e-02 0.0081491800 1.046753e-02
#[2,]  0.004203793 -2.205594e-05 0.0003128027 1.273406e-05
	
	ro <- lapply(r, function(x) {origin(x) <- origin(r[[1]]); x})
	
	rall <- merge(ro[[1]], ro[[2]], ro[[3]], ro[[4]])
}

writeRaster(rall, filename = file.path(dir_grid, "res_1to32_degree", "extent.raster.r.all.grd"), datatype = "INT1S")


#---Aggregate all rasters by 8 to match simulated gridded points
for (i in regions) writeRaster(aggregate(r[[i]], fact = 8), filename = file.path(dir_grid, "res_1to4_degree", paste0("extent.raster.r", i, ".grd")), datatype = "INT1S")
writeRaster(aggregate(rall, fact = 8), filename = file.path(dir_grid, "res_1to4_degree", "extent.raster.r.all.grd"), datatype = "INT1S")


