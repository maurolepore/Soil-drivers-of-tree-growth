# Utility functions

GetGrowthModelData <- function( noSwamp=F )
{
  # Kriged soil data
  df.soil <- read.csv( paste0( dataPath, "kriged soils-10m.csv" ) )
  
  # Get the presaved growth data with light
  df <- readRDS( paste( growthDataPath, "DBH growth data with light-2005-2010.rds", sep="" ) )
  
  # Remove the 20m edges of the plot
  df <- subset( df, binX >= 20 & binX <= 980 & binY >= 20 & binY <= 480 )
  df <- merge( df, df.soil, by.x=c("binX","binY"), by.y=c("x","y") )
  
  if ( noSwamp ) {
    # Remove the swamp
    df.swamp <- read.csv( paste0( dataPath, "BCI swamp.csv" ) )
    library(spatstat)
    w <- owin( poly=df.swamp )
    sel <- inside.owin( df$gx, df$gy, w )
    df.tmp <- df[ sel, ]
    df <- df[ !( df$tag %in% df.tmp$tag ), ]
  }
  
  # Filter out those species with too few (<=10) individuals
  tab <- table( df$sp )
  species <- names( tab[ tab > 10 ] )
  df <- subset( df, sp %in% species )
  
  # Standardise all predictor variables
  df$M3P <- scale( df$M3P )
  df$Mn3 <- scale( df$M3Mn )
  df$Al3 <- scale( df$M3Al )
  df$K3 <- scale( df$M3K )
  df$pH <- scale( df$pH_water )
  df$l.slope <- scale( log(df$slope ) ) # log-transformed slope due to log-normal distribution
  df$light <- rowMeans( df[ , c("light2005","light2006","light2007","light2008","light2009","light2010")] )
  df$l.light <- scale( log(df$light) )
  df$l.dbh <- scale( log(df$dbh) )
  
  # Reduce the size of the df so that the model object isn't bloated unnecessarily
  df <- df[, c("tag", "sp", "incgr", "l.slope", "l.dbh", "dbh", "l.light", "pH", "M3P", "Mn3", "Al3","K3","convex" ) ]
  df
}

# GetDistributionAffinities
# Returns a df with the distributional affinities of tree species across the Panama isthmus,
# as identified in the Condit et al 2013 PNAS paper, for the given tree species.
# The affinities correspond to the first order affinities, not the second order ones.
GetDistributionAffinities <- function( species )
{
  require(XLConnect) 
  
  wb <- loadWorkbook( paste( STRIRawDataPath, "TreeCommunityDrySeasonSpeciesResponse.xlsx", sep="") )
  df <- readWorksheet( wb, sheet = "TreeCommunityDrySeasonSpeciesRe" )
  
  df.result <- data.frame()
  for ( i in 1:length(species) ) {
    df.tmp <- df[ grepl( species[i], df$Latin), 3:11]
    if ( nrow(df.tmp) > 0 ) {
      df.tmp$sp <- species[i]
      df.tmp <- df.tmp[, c(ncol(df.tmp), 1:(ncol(df.tmp) - 1))]
      df.result <- rbind( df.result, df.tmp )
    }
  }
  df.result
}

# GetIsthmusGrowthCoefs
# Returns a df with the growth coefs from Brenes-Arguedas et al of tree species across 
# the Panama isthmus, for the given tree species.
# The coefs correspond to the first order affinities, not the second order ones.
GetIsthmusGrowthCoefs <- function( species )
{
  require(XLConnect) 
  
  wb <- loadWorkbook( paste0( STRIRawDataPath, "Panama Plot Growth Coefs Table S5.xlsx") )
  df <- readWorksheet( wb, sheet = "Table S4" )
  # Clean up the bottom of the worksheet and remove columns not wanted
  df <- df[ -((nrow(df)-19):nrow(df)),-(11:17)]
  
  df.result <- data.frame()
  for ( i in 1:length(species) ) {
    df.tmp <- df[ grepl( species[i], df$Species), 4:10]
    if ( nrow(df.tmp) > 0 ) {
      df.tmp$sp <- species[i]
      df.tmp <- df.tmp[, c(ncol(df.tmp), 1:(ncol(df.tmp) - 1))]
      df.result <- rbind( df.result, df.tmp )
    }
  }
  df.result
}

GetQuadratIndices <- function( x, y, gridSize=20 )
{
  xCentre <- floor( x / gridSize )*gridSize + gridSize*0.5
  yCentre <- floor( y / gridSize )*gridSize + gridSize*0.5
  return( data.frame( binX=xCentre, binY=yCentre ) )
}

GetModelSoilTopoCorrelations <- function()
{
  df.slope <- read.csv( paste0( dataPath, "BCI topography.csv" ) )
  df.soil <- read.csv( paste0( dataPath, "kriged soils-10m.csv" ) )

  df <- merge( df.soil, df.slope, by.x=c("x","y"), by.y=c("binX", "binY") )
  
  vars <- c( "M3P", "M3Mn", "M3Al", "M3K", "convex", "slope", "pH_water", "M3Ca" )
  logs <- c( F, F, F, F, F, T, F, F )
  df <- df[ , vars ]
  for ( i in 1:length(vars) ) {
    if ( !logs[i] ) df[,vars[i]] <- scale( df[,vars[i]] )
    else df[,vars[i]] <- scale( log( df[,vars[i]] ) )
  }
  cor( df )
}

r2.Xu <- function( model )
{
  1 - var(residuals(model)) /
    var(getResponse(model))
}



