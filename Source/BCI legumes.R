# Script to get the legume species in the BCI plot

GetBCILegumes <- function( abbreviation=T, verbose=F )
{
  df.legumes <- read.csv( paste0( dataPath, "Leguminosae.csv" ) ) # All known Fabaceae species
  df.latin <- readRDS( paste0( dataPath, "bci_spptable.rds" ) )
  if ( verbose ) {
    str(df.latin)
    str(df.legumes)
  }
  
  allLegumes <- levels( factor(df.legumes$Genus) )
  df.BCIlegumes <- df.latin[ df.latin$Genus %in% df.legumes$Genus, ]
  
  if( abbreviation ) {
    return( df.BCIlegumes$sp )
  } else {
    return( df.BCIlegumes$Latin )
  }
}