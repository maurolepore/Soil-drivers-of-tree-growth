# Script to get the palm species in the BCI plot

GetBCIPalms <- function( verbose=F )
{
  df.palms <- read.csv( paste0( dataPath, "Arecaceae.csv" ) ) # All known palm species
  df.latin <- readRDS( paste0( dataPath, "bci_spptable.rds" ) )
  if ( verbose ) {
    str(df.latin)
    str(df.palms)
  }
  
  allPalms <- levels( factor(df.palms$Genus) )
  df.BCIpalms <- df.latin[ df.latin$Genus %in% df.palms$Genus, ]
  
  df.BCIpalms$Latin
}