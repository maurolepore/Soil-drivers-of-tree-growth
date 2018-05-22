# Script to install and include the required packages associated with the analyses for the
# paper: Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest
#
# This script handles installing packages that are used and that may not have been
# installed.

# I use a modified version of a function from stackoverflow

InstallPackages <- function( requiredPackages )
{
  remaining <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
  
  if ( length( remaining ) ) {
    install.packages( remaining, dependencies=T )
  }
}

# These are the packages that are needed
requiredPackages = c( "car", "effects", "lazyeval", "lsmeans", "MuMIn", "nlme",
                      "plyr", "reshape2", "spatstat", "stringr" )

# Install the packages if need be
InstallPackages( requiredPackages )
