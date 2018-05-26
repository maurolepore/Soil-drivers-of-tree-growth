# Script to run all the analyses for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest

# All the analyses for the BCI growth model are run, including
# creating all the growth data

# Set up the paths and base environment if need be
if ( !exists("PathsAndEnvironmentInitialised") ) {
  if ( !exists("baseRSourcePath") ) {
    # Set up the variable baseRSourcePath to the correct path
    # e.g.:
    # baseRSourcePath <- Sys.getenv("R_USER") 
  }
  source( paste0( baseRSourcePath, "Source/Env.R" ) )
}

# Set up directories to reference the raw and the reformatted and created growth data
growthDataPath <- paste0( dataPath, "Growth Data/" )
rawDataPath <- paste0( dataPath, "Raw Data/" )

# 1. Load in the CTFS package
attach( paste0( rawDataPath, "CTFSRPackage.rdata" ) )

# 2. Ensure all data is there and in the right format
source( paste0( mainSourcePath, "Reformat raw data.R" ) )

# 3. Create the growth data files if they don't exist
source( paste0( mainSourcePath, "Create BCI growth and light data.R" ) )

# 4. Krige the soils, if not already done so
source( paste0( mainSourcePath, "Krige soils.R" ) )

# 5. Do the analysis to find the best/optimum random effects
source( paste0( mainSourcePath, "Find best model random effects.R" ) )

# 6. Do the backward selection models
source( paste0( mainSourcePath, "Stepwise model selection.R" ) )

# 7. Do the model averaging
source( paste0( mainSourcePath, "Model averaging.R" ) )

##### Additional analyses ####
# 8. Response of leguminous species to P and Mn
source( paste0( mainSourcePath, "Legume Growth Analysis.R" ) )

# 9. Correlation tables and analyses presented in the paper
source( paste0( mainSourcePath, "Correlations Analyses.R" ) )

