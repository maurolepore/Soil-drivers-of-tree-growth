# Growth functions
# Functionality pertaining to calculating the DBH growth from CTFS data

source( paste0( mainSourcePath, "BCI palms.R" ) )

# CalculateDBHGrowth
# Calculates the DBH growth for all living species in the censuses
# contained in the rds files path1 and path2.
# doButtressCorrection: boolean indicating whether buttress corrections 
#                       should be done
# excludedSpecies are species (as binomials) to be excluded from the analysis, e.g. palms
# pomCut and excludeStemChange are passed to growth.indiv
# imposeMinGrowth sets the minimum growth to that given
CalculateDBHGrowth <- function( path1, path2, pathSpecies,
        doButtressCorrection=T, excludedSpecies="", pomCut=1, excludeStemChange=F,
        imposeMinGrowth=NULL )
{
  if ( doButtressCorrection ) {
    # Do the buttress correction
    df.first <- ApplyButtressCorrection( readRDS(path1) )
    df.second <- ApplyButtressCorrection( readRDS(path2) )
  } else {
    df.first <- readRDS(path1)
    df.second <- readRDS(path2)
  }
  
  # Do the calculation
  df <- growth.indiv.ex( df.first, df.second, pomcut=pomCut,
                         excludeStemChange=excludeStemChange )
  
  # Filter out the results based on status/NAs
  selection <- df$status1 == "A" & df$status2 == "A"
  df <- df[ selection, ]
  
  # NA's in the growth column are now of no use
  df <- df[ !is.na(df$incgr), ]
  
  if ( !is.null(imposeMinGrowth) ) {
    df[ df$incgr <= 0, "incgr" ] <- imposeMinGrowth
  }
  
  # Add the Latin name
  df <- merge( df, readRDS(pathSpecies)[ , c("sp","Latin") ] )
  
  # Remove species to be excluded, based on the species binomial
  df <- df[ !(df$Latin %in% excludedSpecies), ]
  
  df
}

# growth.indiv.ex
# Modified form of growth.indiv from the CTFS package to allow
# passing pomcut and excludeStemChange parameters to the 
# trim.growth function
# pomcut gets passed on to trim.growth.ex as the homcut parameter; pom and hom
# represent the same thing but have different precision (hom is more precise)
growth.indiv.ex <- function(census1, census2, rounddown=FALSE, 
                            mindbh=10, dbhunit="mm", err.limit= 4,
                            maxgrow=75, pomcut=0.05, excludeStemChange=T ) 
{
  if (is.null(census2$codes)) 
    census2$codes = rep(".", length(census2$dbh))
  time = (census2$date - census1$date)/365.25
  if (rounddown) {
    sm = ((census1$dbh < 55 | census2$dbh < 55) & !is.na(census1$dbh) & 
            !is.na(census2$dbh))
    census1$dbh[sm] = rndown5(census1$dbh[sm])
    census2$dbh[sm] = rndown5(census2$dbh[sm])
  }
  incgrowth = (census2$dbh - census1$dbh)/time
  expgrowth = (log(census2$dbh) - log(census1$dbh))/time
  good = trim.growth.ex(census1, census2, time, err.limit=err.limit, 
                     maxgrow=maxgrow, mindbh=mindbh, dbhunit=dbhunit,
                     homcut=pomcut, exclude.stem.change=excludeStemChange)
  good[census1$dbh < mindbh] = FALSE
  incgrowth[!good] = expgrowth[!good] = NA
  growthdata = data.frame(tag = I(census1$tag), #treeID = census1$treeID, Yasuni/conditdatacenter doesn't have treeID
                          sp = I(census1$sp), gx = census1$gx, gy = census1$gy, 
                          dbh1 = census1$dbh, dbh2 = census2$dbh, time = time, 
                          incgr = incgrowth, expgr = expgrowth, status1 = census1$status, 
                          status2 = census2$status)
  return(growthdata)
}

# trim.growth.ex
# Modified version of trim.growth from the CTFS package, modified to exclude trees
# based on changes in hom, rather than pom. This is necessary for data from
# conditdatacenter.org, which does not have the pom column
trim.growth.ex <- function( cens1, cens2, time, slope=0.006214, intercept=.9036,
                            err.limit=4, maxgrow=75, homcut=0.05, mindbh=10,
                            dbhunit='mm', exclude.stem.change=TRUE)
{
  if (dbhunit=='cm') intercept = intercept/10
  stdev.dbh1 = slope*cens1$dbh + intercept
  
  growth = (cens2$dbh-cens1$dbh)/time
  
  bad.neggrow = which(cens2$dbh <= (cens1$dbh - err.limit*stdev.dbh1)) 
  bad.posgrow = which(growth>maxgrow)
  
  homdiff = abs( as.numeric(cens2$hom) - as.numeric(cens1$hom)) / as.numeric(cens1$hom)
  
  accept = rep(TRUE,length(cens1$dbh))
  accept[homdiff>homcut] = FALSE
  accept[bad.neggrow] = FALSE
  accept[bad.posgrow] = FALSE
  accept[is.na(growth)] = FALSE
  
  if(exclude.stem.change) accept[cens1$stemID != cens2$stemID] = FALSE
  
  accept[cens1$dbh<mindbh] = FALSE
  accept[is.na(cens1$dbh) | is.na(cens2$dbh) | cens2$dbh<=0] = FALSE
  
  return(accept)
}


# GetDBHGrowthEx
# Uses all the CTFS functions (or slight mods thereof) to calculate
# the DBH growth for the given censuses
GetDBHGrowthEx <- function( censusPath1, censusPath2, speciesPath, 
                            topoPath=NULL,
                            imposeMinGrowth=NULL,
                            minimumNumber=10,
                            pomCut=1,
                            gridSize=10,
                            plotSize=c(1000,500) )
{
  palmSpecies <- GetBCIPalms()
  df <- CalculateDBHGrowth( censusPath1, censusPath2, pomCut=pomCut,
                            pathSpecies=speciesPath,
                            doButtressCorrection=F, excludedSpecies=palmSpecies,
                            imposeMinGrowth=imposeMinGrowth )   
  
  # Filter out those with too few individuals
  tab <- table( df$sp )
  species <- names( tab[ tab > minimumNumber ] )
  df <- subset( df, sp %in% species )
  
  # Add in the grid coords which correspond to the quadrats from the kriged soil
  df.grid <- GetQuadratIndices( df$gx, df$gy, gridSize )
  df <- cbind( df, df.grid )
  
  # Do a correction for individuals on the outer border - they'll be placed in the last quadrat
  df[ df$gx == plotSize[1], "binX" ] <- plotSize[1] - gridSize*0.5
  df[ df$gy == plotSize[2], "binY" ] <- plotSize[2] - gridSize*0.5
  
  # Add topographic data?
  if ( !is.null(topoPath) ) {
    df.topo <- readRDS( topoPath )
    df.slope <- allquadratslopes( df.topo, gridsize=gridSize, plotdim=plotSize, edgecorrect=TRUE)
    # df.slope doesn't have the quadrat x,y coords so I need to add them. The function used
    # rowcol.to.index, so I'm following its format which takes the y coord first
    df.slope <- cbind( df.slope, expand.grid( binY=seq(gridSize/2,500-gridSize/2,by=gridSize), 
                                              binX=seq(gridSize/2,1000-gridSize/2,by=gridSize) ) )
    df <- merge( df, df.slope, by=c("binX", "binY") )
  }
  
  df  
}

GetQuadratIndices <- function( x, y, gridSize=20 )
{
  xCentre <- floor( x / gridSize )*gridSize + gridSize*0.5
  yCentre <- floor( y / gridSize )*gridSize + gridSize*0.5
  return( data.frame( binX=xCentre, binY=yCentre ) )
}

# ApplyButtressCorrection
# Returns the df with the DBHs corrected for buttress trees, i.e. those measured above
# 1.3 m. The corrections are from Cushman et al. 2014, using equation 1
ApplyButtressCorrection <- function( df, dbh="dbh", hom="hom", id="treeID" )
{
  # There's already a table of corrections for some species:
  df.params <- read.csv( paste0( dataPath, "taper.parameters.csv" ) )
  
  # Only use equation 1 for the correction, as it was the most accurate
  df.params <- df.params[ df.params$eqn == 1, c( "sp", "DBH", "b1" ) ]
  
  # For the species covered by the buttress correction study, the mean (modelled) values
  # for the species should be used. For all others, the overall mean correction is used
  require(plyr)
  df.means <- ddply( df.params, .(sp), numcolwise(mean) )
  
  selection <- ( df$sp %in% levels(df.means$sp) ) & is.finite(df[,hom]) & df[,hom] > 1.3
  
  # Calculate and merge these values
  df.tmp <- df[ selection, ]
  df.tmp <- merge( df.tmp, df.means[, c("sp", "b1")], by="sp" )
  
  # This is the actual correction calculation: From equation 1 of Cushman et al 2014,
  # d = DBH * exp(-b1 * (h - 1.3))
  # => DBH = d * exp(b1 * (h - 1.3))
  # I round the values because the additional precision is beyond the original
  df.tmp$DBH <- round( df.tmp[,dbh] * exp( df.tmp$b1 * (df.tmp[,hom] - 1.3) ) )
  
  # Merge those corrected values back into the original data
  # The following selection, if used with a stem table rather than a tree table, assumes
  # that buttress trees only have one buttressed stem (because the selection is done by
  # tree ID and doesn't include the stem ID)
  selection <- match( df.tmp[,id], df[,id] )
  df[ selection, dbh ] <- df.tmp$DBH
  
  # The mean values need to be used for all other species with dbh > 1.3
  selection <- !( df$sp %in% levels(df.means$sp) ) & is.finite(df[,hom]) & df[,hom] > 1.3
  b1 <- mean( df.means$b1 )
  df[ selection, dbh ] <- round( df[selection, dbh] * exp( b1 * (df[selection, hom] - 1.3) ) )
  
  df
}

# taper.H
# Height allometry calculation from the Cushman 2014 paper, modified to do the
# calculation from dbh in mm, rather than cm
# df.ssa allows species-specific allometries to be used to improve the estimations
# for certain species
# dbh: the name of the column with the DBHs in it.
taper.H <- function( df, dbh="dbh", df.ssa=NULL )
{
  if ( is.null(df.ssa) ) {
    # Original formula for dbh in cm: 43.4375 * ( 1- exp(-0.04469 * dbh^0.78339) ) 
    return( 43.4375 * ( 1- exp(-0.007359027 * df[,dbh]^0.78339) ) )
  } else {
    df.tmp <- df[ , c("sp", dbh) ]
    df.tmp$h <- 43.4375 * ( 1 - exp(-0.007359027 * df[,dbh]^0.78339) )
    for ( i in 1:nrow(df.ssa) ) {
      selection <- df.tmp$sp == df.ssa$sp[i]
      df.tmp$h[ selection ] <- df.ssa$a[i] * 
        ( 1 - exp(-df.ssa$b[i] * (df[selection,dbh]*0.1) ^ df.ssa$c[i] ) )
    }
    return( df.tmp$h )
  }
}
