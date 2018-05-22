# Script to compare the responses of legume species versus non-legumes for the paper:
# Zemunik et al. (2018) Soil drivers of local-scale tree growth in a lowland tropical forest


library( nlme )
source( paste0( mainSourcePath, "BCI legumes.R" ) )
mod.c.reml <- readRDS( paste0( dataPath, "Selected growth model.rds" ) )

# The random effects coefficients in this model are for four REs: light, DBH, Mn and P
df.coef <- coef( mod.c.reml )
df.legumes <- df.coef[ row.names(df.coef) %in% GetBCILegumes(), ]
df.nonlegumes <- df.coef[ !(row.names(df.coef) %in% GetBCILegumes()), ]
df.coef$legume <- F
df.coef[ row.names(df.coef) %in% GetBCILegumes(), "legume" ] <- T

mod.Mn <- gls( Mn3 ~ legume, data=df.coef )
plot( mod.Mn )
summary( mod.Mn )
mod.Mn.gls <- gls( Mn3 ~ legume, data=df.coef, weights = varExp() )
mod.Mn.gls.2 <- gls( Mn3 ~ legume, data=df.coef, weights = varPower(2) )
AIC( mod.Mn, mod.Mn.gls, mod.Mn.gls.2 ) # varExp or varPower is better
plot(mod.Mn.gls)
anova(mod.Mn.gls)

mod.P <- gls( M3P ~ legume, data=df.coef )
mod.P.gls <- gls( M3P ~ legume, data=df.coef, weights = varExp() )
AIC( mod.P, mod.P.gls ) # varExp is better
plot(mod.P.gls)
anova(mod.P.gls)

# The main result reported in the paper is the different responses to P and Mn
# in the legume species
print( summary(mod.P.gls) ) # p = 0.24
print( summary(mod.Mn.gls) ) # Legumes more negative by -0.011 (p=0.004)

# Light and DBH analyses
mod.light <- gls( l.light ~ legume, data=df.coef )
mod.light.gls <- gls( l.light ~ legume, data=df.coef, weights = varExp() )
AIC( mod.light, mod.light.gls )
summary( mod.light.gls ) # p = 0.7854

mod.dbh <- gls( l.dbh ~ legume, data=df.coef )
mod.dbh.gls <- gls( l.dbh ~ legume, data=df.coef, weights = varExp() )
AIC( mod.dbh, mod.dbh.gls )
summary( mod.dbh.gls ) # p = 0.1082
