# Standalone script to run all the analyses and create the various source and 
# result model data frames

if (Sys.getenv("R_USER") == "") {
  basePath <- "~"   # The typical R base path on Unix/Mac systems
} else {
  basePath <- Sys.getenv("R_USER")  # Windows
}

# This path needs to be set to the correct path
baseRSourcePath <- paste0( basePath, "/R/Github/Soil-drivers-of-tree-growth/" )

source( paste0( baseRSourcePath, "Source/Env.R" ) )
source( paste0( mainSourcePath, "packages.R" ) )

# Do the analyses
source( paste0( mainSourcePath, "All analyses.R" ) )
