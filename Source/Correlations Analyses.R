# Script to compute correlations of soil variables and correlations of P 
# affinities (Condit et al 2013) against P responses, for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest

source( paste0( mainSourcePath, "Utilities.R" ) )


############# Soil variables correlations
# Presented as Table S1 in Appendix S1
print( GetModelSoilTopoCorrelations() )



############# P affinities comparison
mod <- readRDS( paste0( dataPath, "Selected growth model.rds" ) )

# The random effects coefficients in this model are for four REs: light, DBH, Mn and P
df.coef <- coef( mod )

# The affinity table only has Latin names, so to be able to merge it with the responses
# we need to get the Latin names from the species table
df.latin <- readRDS( paste0( growthDataPath, "bci_spptable.rds" ) )

# Merge the affinities with the species responses
df.aff <- data.frame( sp=row.names( df.coef ) )
df.aff <- cbind( df.aff, df.coef[df.aff$sp,] )
df.aff <- merge( df.aff, df.latin[, 1:2], by="sp")
df.aff <- merge( df.aff, GetDistributionAffinities(df.aff$Latin), by.x="Latin", by.y="sp" )

# Inspect correlations
# Graphically
# pairs(df.aff[, c("M3P", "P", "Mn3", "Al") ]) # no correlation in affinites
# Matrix
affCor <- cor(df.aff[, c("M3P", "P", "Mn3", "Al") ])

# R2 for P affinity
print( paste( "P affinity correlation R2 with P response: ", affCor["M3P", "P"]^2) )  # 0.0126

# Assess the relationship with a model
mod.P <- lm( M3P ~ P, data=df.aff)
# plot( mod.P )
summary( mod.P )  # No relationship for P, p = 0.11, adjusted R2 0.0076, multiple R2 as above


