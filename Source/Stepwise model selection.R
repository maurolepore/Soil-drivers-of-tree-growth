# Stepwise selection of the "best" model using the 2010 soil sampling.
# All variables are standardised and 20 m around the plot edge is excluded to
# eliminate edge effects with the light proxy

library( nlme )
library( car )

source( paste( mainSourcePath, "Utilities.R", sep="") )
df <- GetGrowthModelData()  # From Utilities.R

# The models need to be fitted with ML since the fixed effects
# are being changed in each step of the process

# The "optimal" random effects included the soil vars Mn and P.

fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + pH + l.dbh + l.light + l.slope + convex
randomEffect <- ~ l.dbh + l.light + Mn3 + M3P | sp
modCtrl <- lmeControl( maxIter=150, msMaxIter=150, niterEM=50, msMaxEval=250 )
startTime <- proc.time()
mod.full <- lme( fixed=fixedEffect, random=randomEffect,
                data=df, method="ML", control=modCtrl,
                weights=varExp() )
print( proc.time() - startTime )
summary( mod.full ) # Slope least significant

# Remove slope
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + pH + l.dbh + l.light + convex
randomEffect <- ~ l.dbh + l.light + Mn3 + M3P | sp
startTime <- proc.time()
mod.a <- lme( fixed=fixedEffect, random=randomEffect,
                      data=df, method="ML", control=modCtrl,
                      weights=varExp() )
print( proc.time() - startTime )
summary(mod.a) # pH least significant

# Remove pH
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + l.dbh + l.light + convex
startTime <- proc.time()
mod.b <- lme( fixed=fixedEffect, random=randomEffect,
                data=df, method="ML", control=modCtrl,
                weights=varExp() )
print( proc.time() - startTime )
summary(mod.b) # Convexity least significant

# Remove convex
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + l.dbh + l.light
startTime <- proc.time()
mod.c <- lme( fixed=fixedEffect, random=randomEffect,
                data=df, method="ML", control=modCtrl,
                weights=varExp() )
print( proc.time() - startTime )
summary(mod.c) # All signif
anova( mod.b, mod.c ) # 2.c better, as expected

vif( mod.c ) # K and Al have highest vif's - Al the highest

# Test removing Al due to collinearity.
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + l.dbh + l.light
startTime <- proc.time()
mod.d <- lme( fixed=fixedEffect, random=randomEffect,
                data=df, method="ML", control=modCtrl,
                weights=varExp() )
print( proc.time() - startTime )
summary(mod.d) # P not signif
anova( mod.b, mod.c, mod.d ) # mod.d much worse by AIC and BIC
# suggesting that maintaining the collinear vars is better.

# REML model
fixedEffect <- log(incgr) ~ Mn3 + M3P + K3 + Al3 + l.dbh + l.light
startTime <- proc.time()
mod.c.reml <- lme( fixed=fixedEffect, random=randomEffect,
                data=df, method="REML", control=modCtrl,
                weights=varExp() )
print( proc.time() - startTime )
summary(mod.c.reml)

# Save the model? If not, comment out
saveRDS( mod.c.reml, paste0( dataPath, "Selected growth model.rds" ) )
df.coef <- cbind( data.frame( sp=row.names(coef( mod.c.reml )) ), coef( mod.c.reml ) )
write.csv( df.coef, paste0( dataPath, "BCI growth model Coefs.csv" ), row.names = F )


