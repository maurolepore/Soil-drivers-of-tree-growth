# Script to calculate (and save) the DBH growth and canopy light data for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest
# Various paths, e.g. mainSourcePath, should already have been set up

source( paste0( mainSourcePath, "Growth functions.R" ) )
source( paste0( mainSourcePath, "Shade estimation.R" ) )

censusPath1 <- paste0( growthDataPath, "BCI tree census 6.rds" )
censusPath2 <- paste0( growthDataPath, "BCI tree census 7.rds" )
sppPath <- paste0( growthDataPath, "bci_spptable.rds" )
topoPath <- paste0( growthDataPath, "BCI topography.rds" )

if ( !file.exists( paste0( growthDataPath, "dbhgrowth-tree censuses 6-7.rds" ) ) ) {
  # Calculate the DBH growth. Note that for the data from conditdatacenter the
  # tables do not have the pom column, so pomCut refers to homCut.
  # Because buttressed trees might have substantial hom differences, and the
  # buttressing is being corrected for, the homCut is relaxed somewhat, 0.5
  # vs 0.05
  df <- GetDBHGrowthEx( censusPath1, censusPath2, speciesPath=sppPath, 
                        topoPath=topoPath,
                        imposeMinGrowth=0.1, minimumNumber=10,
                        pomCut=0.5, gridSize=10 )

  # Set the dbh per growth interval to be the midpoint
  df$dbh <- (df$dbh1 + df$dbh2) / 2

  # Clean up the table - remove unnecessary columns
  df <- df[ , c("binX", "binY", "gx", "gy", "tag", "sp", "Latin",
                "incgr", "convex", "slope", "dbh") ]

  saveRDS( df, paste0( growthDataPath, "dbhgrowth-tree censuses 6-7.rds" ) )
}

if ( !file.exists( paste0( growthDataPath, "DBH growth data with light-2005-2010.rds" ) ) ) {
  # In order to filter out a basic subset of individuals for which light estimation
  # is problematic, e.g. palms, I use the DBH growth function even though I'm
  # not actually interested in the growth per se, just the tree identities.
  # homCut and the minimum number per species are thus relaxed
  df <- GetDBHGrowthEx( censusPath1, censusPath2, speciesPath=sppPath, 
                        topoPath=topoPath,
                        imposeMinGrowth=0.1, minimumNumber=0,
                        pomCut=2, gridSize=10 )
  
  # Set the dbh per growth interval to be the midpoint
  df$dbh <- (df$dbh1 + df$dbh2) / 2
  
  # Calculate the tree heights
  df.ssa <- read.csv( paste0( rawDataPath, "BCI species-specific allometries.csv" ) )
  df$height <- taper.H( df, df.ssa=df.ssa )
  
  listCanopy <- GetCanopyCensuses()
  start.0 <- proc.time()
  # For manual parallelisation, a subset of censuses can be used and later merged
  censusSubset <- 4:6  
  for ( i in censusSubset ) {
    start <- proc.time()
    s <- paste0( "shade", names(listCanopy)[i] )
    l <- paste0( "light", names(listCanopy)[i] )
    print( s )
    df[,s] <- GetShadingIndex( df, listCanopy[[i]] )
    # Convert the shade into a light value using the formula from Appendix S1 of
    # Ruger et al. 2009
    df[,l] <- exp(-0.01351 - 0.08043*df[,s] - 0.00315*df[,s]^2)
    print( proc.time() - start )
    saveRDS( df[,c("tag",s,l)], paste0( growthDataPath, "Shade-light-", names(listCanopy)[i], ".rds" ) )
  }
  print( proc.time() - start.0 )

  # Save the full df, after merging the shade values into the main
  # growth table
  df.growth <- readRDS( paste0( growthDataPath, "dbhgrowth-tree censuses 6-7.rds" ) )
  files <- dir( growthDataPath )[ grepl( "Shade\x2Dlight", dir(growthDataPath) ) ]
  for ( i in 1:length(files) ) {
    df.tmp <- readRDS( paste0( growthDataPath, files[i] ) )
    cols <- names(df.tmp)[ grepl("shade|light", names(df.tmp)) ]
    df.growth <- merge( df.growth, df.tmp[,c("tag",cols)], by="tag" )
  }
  saveRDS( df.growth, paste0( growthDataPath, "DBH growth data with light-2005-2010.rds" ) )
}



