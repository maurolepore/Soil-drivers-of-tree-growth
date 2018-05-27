# Script to krige the soil variables used in the growth model for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest


# Check if the kriged soils file does not exist
if ( !file.exists( paste0( dataPath, "kriged soils-10m.csv" ) ) ) {

  # The kriging functionality resides in a package, managed by ForestGEO.
  # Install the package if need be
  #devtools::install_github("forestgeo/fgeo.habitat")
  
  # Get the soil data
  df.soil <- read.csv( paste0( rawDataPath, "soils.csv" ) )
  # Get the kriging semivariogram parameters
  df.prms <- read.csv( paste0( dataPath, "krige prms.csv" ) )
  
  df.res <- data.frame()
  for ( i in 1:nrow(df.prms) ) {
    prm <- df.prms[i,]
    prms <- list( model=as.character(prm$model),
                  kappa=prm$kappa,
                  nugget=prm$nugget, sill=prm$sill,
                  range=prm$range )
    # Do the kriging. Note that the argument useKsLine=T or F should give the same result
    ks <- fgeo.habitat::GetKrigedSoil( df.soil, var=as.character(prm$var),
                                       krigeParams=prms, gridSize=10, useKsLine=F )
    
    # Store the results in a data frame
    if ( i == 1 ) {
      df.res <- ks$df 
    } else {
      df.res <- cbind( df.res, ks$df[,"z"] )
    }
  }
  
  # Set the column names correctly
  names( df.res )[3:ncol(df.res)] <- as.character( df.prms$var )
  
  # Save the kriged values
  write.csv( df.res, file=paste0( dataPath, "kriged soils-10m.csv" ), row.names=F )
}