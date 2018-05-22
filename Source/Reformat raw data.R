# Script to create rds files from the raw data for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest


if ( !exists("rawDataPath") ) stop( "Please define the raw data path in the variable rawDataPath" )

CreateRDSObject <- function( inPath, outPath )
{
  attach_if_needed( list( inPath ) )
  name <- ls(2)
  saveRDS( get( name ), outPath )
  rm( name )
}

if ( !file.exists( paste0( growthDataPath, "BCI tree census 6.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "bci.tree6.rdata" ),
                   paste0( growthDataPath, "BCI tree census 6.rds" ) )
}

if ( !file.exists( paste0( growthDataPath, "BCI tree census 7.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "bci.tree7.rdata" ),
                   paste0( growthDataPath, "BCI tree census 7.rds" ) )
}

# Legacy files - i.e. Data obtained from the CTFS servers, not from 
# conditdatacenter.org
if ( !file.exists( paste0( growthDataPath, "BCI full census 6.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "bci.full6.rdata" ),
                   paste0( growthDataPath, "BCI full census 6.rds" ) )
}

if ( !file.exists( paste0( growthDataPath, "BCI full census 7.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "bci.full7.rdata" ),
                   paste0( growthDataPath, "BCI full census 7.rds" ) )
}

if ( !file.exists( paste0( growthDataPath, "bci_spptable.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "bci.spptable.rdata" ),
                   paste0( growthDataPath, "bci_spptable.rds" ) )
}

if ( !file.exists( paste0( growthDataPath, "BCI topography.rds" ) ) ) {
  CreateRDSObject( paste0( rawDataPath, "CTFSElev_bci.rdata" ),
                   paste0( growthDataPath, "BCI topography.rds" ) )
}
