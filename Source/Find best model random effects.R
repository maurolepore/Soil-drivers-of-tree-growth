# Script to estimate the best/optimal random effects for the DBH growth model
# using the 2010 soil sampling by Wolf with light estimates using Ruger 2009
# and edge effects (bordering 20m) removed

library(nlme)

source( paste( mainSourcePath, "Utilities.R", sep="") )
df <- GetGrowthModelData()  # From Utilities.R

# Find the best random effects using REML. Zuur (2009) advocates using
# maximal random effects to begin with and then backward select the REs.
# I start by using the full model with all predictors for the fixed effects
# variables. For the REs, the full set is those with at least some support
# from data exploration: DBH, light, Mn, P, K
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + pH + l.dbh + l.light + l.slope + convex
modCtrl <- lmeControl( maxIter=150, msMaxIter=150, niterEM=50, msMaxEval=250 )

listRE <- list( "~ l.dbh + l.light + Mn3 + M3P + K3 | sp",
                "~ l.dbh + l.light + Mn3 + K3 | sp",
                "~ l.dbh + l.light + Mn3 + M3P | sp",
                "~ l.dbh + l.light + M3P + K3 | sp",
                "~ l.dbh + l.light + Mn3 | sp",
                "~ l.dbh + l.light + K3 | sp",
                "~ l.dbh + l.light + M3P | sp",
                "~ l.dbh + l.light | sp")

# Test each RE structure. REML must be used
listMods <- vector( mode="list", length=length(listRE) )
for ( i in 1:length(listRE) ) {
  startTime <- proc.time()
  randomEffect <- as.formula( listRE[[i]] )
  listMods[[i]] <- lme( fixed=fixedEffect, random=randomEffect,
                        data=df, method="REML", control=modCtrl,
                        weights=varExp() )
  print( paste( "Time for", listRE[[i]], ":" ) )
  print( proc.time() - startTime )
}

# Compare them all
anova( listMods[[1]], listMods[[2]], listMods[[3]], listMods[[4]],
       listMods[[5]], listMods[[6]], listMods[[7]], listMods[[8]])

# Summary:
# RE 3 has slightly worse AIC than 2 but a lot better BIC
# So, although a bit equivocal RE 2 could be considered the best option.
# In the model averaging, however, several different RE structures are considered
# to test whether altering these makes substantial differences
