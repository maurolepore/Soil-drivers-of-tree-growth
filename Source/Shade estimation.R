# Routines for estimating light and shade at BCI

# GetCanopyCensuses
# Get the canopy data in a list, with the name of each element being the year
# censusYears controls which censuses to retrieve
GetCanopyCensuses <- function( censusYears=c( "2005", "2006", "2007", "2008", "2009", "2010" ),
                               canopyPath=paste0(dataPath,"Raw Data/") )
{
  listCanopy <- vector( "list", l=length(censusYears) )
  for ( i in 1:length(censusYears) ) {
    listCanopy[[i]] <- read.delim( paste( canopyPath, "Canopy", censusYears[i], ".csv", sep="" ) )
    head( listCanopy[[i]] )
    # Select only the columns we want - the 0_2 column, used until 2003, is not included
    listCanopy[[i]] <- listCanopy[[i]][ , c("x", "y", "ht0_1", "ht1_2", "ht2_5", "ht5_10",
                                        "ht10_20", "ht20_30", "ht30_") ]
  }
  names( listCanopy ) <- censusYears
  listCanopy
}

# I'm using the same model and approximations as the Ruger 2009 paper, but I
# modify the calculations so that the shade/light at the top of each tree is
# estimated, rather than at a fixed 2 m.

# GetObscuredAngle
# Calculates the angle cast by a horizontal disc of radius r to a point lte
# 20 m away.
# d: horizontal distance between the tree in question and the centre of the cell
# h: vertical distance between the top of the tree and the centre of the cell
# r: radius of the shading disc
GetObscuredAngle <- function( d, h, r=2.5)
{
  if ( d > 20 ) return( 0 )
  
  if ( length(h) > 1 ) {
    h <- h[ !is.na(h) ]   # There might be NAs in the canopy data
    h <- h[ h > 0 ]
    if ( length(h) == 0 ) return( 0 )
  }

  atan( (d + r)/h ) - atan( (d - r)/h ) 
}

# GetCanopyCellRange
# Returns a list of ranges that are extracted from cellNames.
# An open range is specified as having an equal upper and lower bounds
GetCanopyCellRange <- function( cellNames )
{
  require( stringr )
  cellNames <- cellNames[ grepl("[0-9]+", cellNames) ]
  ranges <- vector( "list", l=length(cellNames) )
  for ( i in 1:length(cellNames) ) {
    r <- str_extract_all( cellNames[i], "[0-9]+" )
    r1 <- r[[1]][1]
    r2 <- ifelse( length(r[[1]]) > 1, r[[1]][2], r[[1]][1] )
    ranges[[i]] <- as.numeric( c( r1, r2 ) )
  }
  # Return a list of ranges and the names they apply to
  return( list( ranges=ranges, names=cellNames ) )
}

GetQuadrat <- function( x, y, gridSize=20 )
{
  xCentre <- gridSize * (floor( x / gridSize ) + 0.5)
  yCentre <- gridSize * (floor( y / gridSize ) + 0.5)
  c( xCentre, yCentre )
}

# GetShadingIndex
# Returns the sum of all obscured angles, relating to the shade cast by discs 
# approximating layers of leaves in the canopy
# There will be edge effects (reductions) for points near the edge of the plot
# Vegetation data: gx, gy, height of the trees
# Canopy data: x, y, (level occupancy)
GetShadingIndex <- function( df.veg, df.canopy, plotWidth=1000, plotHeight=500 )
{
  # For the canopy data, I use a function to parse the column names and then 
  # return the max height info
  ranges <- GetCanopyCellRange( names(df.canopy) )
  
  # The main code is a loop on every tree in the plot. A more efficient
  # way might be to somehow lump trees in an area, but that could be tricky
  indices <- vector( m="numeric", l=nrow(df.veg) )
  discHeights <- vector( m="numeric", l=length(ranges$ranges) )
  heights <- vector( m="numeric", l=length(ranges$ranges) )
  for ( i in 1:nrow(df.veg) ) {
    indices[i] <- 0
    
    gx <- df.veg$gx[i]
    gy <- df.veg$gy[i]
    quadrat <- GetQuadrat( gx, gy, gridSize=5 )
    
    for ( r in 1:length(ranges$ranges) ) {
      # sum * 0.5 = mean, because n = 2
      discHeights[ r ] <- (sum(ranges$ranges[[r]])*0.5 - df.veg$height[i])
    }
    
    # Loop on all the cells that can give a centre within 20 m of the tree
    # Assuming a cell size of 5 m, the 3 cells around each corner of the grid
    # are guaranteed to be >20 m away
    for ( x in -4:4 ) {
      for ( y in -4:4 ) {
        if ( (abs(x) == 4 & abs(y) >=3) |
             (abs(y) == 4 & abs(x) >=3) ) next # corner cells
        
        # Calculate the actual distance to the centre of the cell and 
        # jump the loop if it's >20 m
        cx <- quadrat[1] + x*5
        cy <- quadrat[2] + y*5
        distance <- sqrt( (gx-cx)^2 + (gy-cy)^2 )
        if ( distance > 20 | cx < 0 | cy < 0 |
             cx > plotWidth | cy > plotHeight ) next

        # Select the canopy census cell for the current location
        m <- df.canopy$x==cx & df.canopy$y==cy
        df.cell <- df.canopy[ m, ranges$name ]
        
        for ( r in 1:length(ranges$ranges) ) {
          # Get all the heights. The height to each "disc" is the same
          # for tree i. If the disc is shaded then its value is 100%. The
          # 0.01 factor converts this to 1
          heights[ r ] <- discHeights[ r ] * df.cell[ 1, r ] * 0.01
        }
        
        # Calculate the index as the sum of the obscured angles
        indices[i] <- indices[i] + sum( GetObscuredAngle( distance, heights ) )
      }
    }
  }

  # Return the indices
  indices
}















