# Standalone script to run all the analyses and create the various source and 
# result model data frames

if (Sys.getenv("R_USER") == "") {
  basePath <- "~"   # The typical R base path on Unix/Mac systems
} else {
  basePath <- Sys.getenv("R_USER")  # Windows
}

# This path needs to be set to the correct path
baseRSourcePath <- paste( basePath, "/R/Github/Soil-drivers-of-tree-growth/", sep="" )

source( paste( baseRSourcePath, "Source/Env.R", sep="" ) )
source( paste( mainSourcePath, "packages.R", sep="" ) )

# Do the analyses
source( paste( mainSourcePath, "All analyses.R", sep="" ) )
