# Model-averaging of the growth model using the 2010 soil sampling by Wolf
# with light estimates using Ruger 2009 and edge effects (bordering 20m)
# removed

library(nlme)
library(MuMIn)
options(na.action = "na.fail")

source( paste0( mainSourcePath, "Utilities.R" ) )
df <- GetGrowthModelData()  # From Utilities.R

# The best random effect structure is found by "Find best model random effects.R"
# Since the "best" RE chosen would not have been the same depending on whether
# the AIC or BIC was used to decide, several different RE structures are tested
# to ascertain any differences based on differing RE structures.

# Model for the model averaging with the "best" RE, and the full fixed effects
fixedEffect <- log(incgr) ~ l.dbh + l.light + Mn3 + M3P + Al3 + K3 + pH + l.slope + convex
randomEffect <- ~ l.dbh + l.light + Mn3 + M3P | sp
modCtrl <- lmeControl( maxIter=150, msMaxIter=150, niterEM=50, msMaxEval=250 )
startTime <- proc.time()
mod.full <- lme( fixed=fixedEffect, random=randomEffect,
                 data=df, method="ML", control=modCtrl,
                 weights=varExp() )
print( proc.time() - startTime )

# Model averaging
startTime <- proc.time()
mods.ldmp.AIC <- dredge( mod.full, subset=l.dbh & l.light & Mn3 & M3P, 
                           beta="none", extra=c(r2.Xu), trace=T )
print( proc.time() - startTime )

ma <- model.avg( mods.ldmp.AIC, cumsum(weight) <= .95, fit = F )
df.ci <- as.data.frame( confint( ma ) ) 
df.ci <- cbind( fit=coef(ma), df.ci)
df.ci <- df.ci[ order(row.names(df.ci)), ]
write.csv( df.ci, file=paste0( dataPath, "BCI dredged growth model 95CI light,dbh,Mn,P-RE 2010 soil.csv") )

# Test the effect of different REs
# RE with light, dbh and Mn (i.e. no P)
fixedEffect <- log(incgr) ~ l.dbh + l.light + Mn3 + M3P + Al3 + K3 + pH + l.slope + convex
randomEffect <- ~ l.dbh + l.light + Mn3 | sp
startTime <- proc.time()
mod.full.2 <- lme( fixed=fixedEffect, random=randomEffect,
                 data=df, method="ML", control=modCtrl,
                 weights=varExp() )
print( proc.time() - startTime )
startTime <- proc.time()
mods.ldm.AIC <- dredge( mod.full.2, subset=l.dbh & l.light & Mn3, 
                         beta="none", extra=c(r2.Xu), trace=T )
print( proc.time() - startTime )
ma.ldm <- model.avg( mods.ldm.AIC, cumsum(weight) <= .95, fit = F )
df.ci.2 <- as.data.frame( confint( ma.ldm ) ) 
df.ci.2 <- cbind( fit=coef(ma.ldm), df.ci.2)
df.ci.2 <- df.ci.2[ order(row.names(df.ci.2)), ]
write.csv( df.ci.2, file=paste0( dataPath, "BCI dredged growth model 95CI light,dbh,Mn-RE 2010 soil.csv") )

# RE with light, dbh
fixedEffect <- log(incgr) ~ l.dbh + l.light + Mn3 + M3P + Al3 + K3 + pH + l.slope + convex
randomEffect <- ~ l.dbh + l.light | sp
startTime <- proc.time()
mod.full.3 <- lme( fixed=fixedEffect, random=randomEffect,
                   data=df, method="ML", control=modCtrl,
                   weights=varExp() )
print( proc.time() - startTime )
startTime <- proc.time()
mods.ld.AIC <- dredge( mod.full.3, subset=l.dbh & l.light, 
                        beta="none", extra=c(r2.Xu), trace=T )
print( proc.time() - startTime )
ma.ld <- model.avg( mods.ld.AIC, cumsum(weight) <= .95, fit = F )
df.ci.3 <- as.data.frame( confint( ma.ld ) ) 
df.ci.3 <- cbind( fit=coef(ma.ld), df.ci.3)
df.ci.3 <- df.ci.3[ order(row.names(df.ci.3)), ]
write.csv( df.ci.3, file=paste0( dataPath, "BCI dredged growth model 95CI light,dbh-RE 2010 soil.csv") )


# Model averaging with the size subsets 
df <- readRDS( paste0( growthDataPath, "DBH growth data with light-2005-2010.rds" ) )
# Remove the edge of the plot
df <- subset( df, binX >= 20 & binX <= 980 & binY >= 20 & binY <= 480 )
# Merge with the soil
df.soil <- read.csv( paste( dataPath, "kriged soils-10m.csv", sep="" ) )
df <- merge( df, df.soil, by.x=c("binX","binY"), by.y=c("x","y") )
# Add in the light estimate
df$light <- rowMeans( df[ , c("light2005","light2006","light2007","light2008","light2009","light2010")] )

listSizes <- list( list( df$dbh < 25, "dbh < 25", "dbh25" ),
                   list( df$dbh >= 25 & df$dbh < 250, "25 >= dbh < 250", "dbh25-250" ),
                   list( df$dbh >= 250, "dbh >= 250", "dbh250" ) )

fixedEffect <- log(incgr) ~ l.dbh + l.light + Mn3 + M3P + Al3 + K3 + pH + l.slope + convex
# RE of light, DBH, Mn, P
randomEffect <- ~ l.dbh + l.light + Mn3 + M3P | sp
mods <- vector(mode="list", length=length(listSizes))
modAvg <- vector(mode="list", length=length(listSizes))

findBestRE <- F # Set to T to find the optimal random effects per size class

for ( i in 1:length(listSizes) ) {
  print( paste( "Iteration:", listSizes[[i]][[2]] ) )
  df.data <- df[ listSizes[[i]][[1]], ]
  # Filter out those with too few individuals
  tab <- table( df.data$sp )
  species <- names( tab[ tab > 10 ] )
  df.data <- subset( df.data, sp %in% species )
  
  # Standardise
  df.data$M3P <- scale( df.data$M3P )
  df.data$Mn3 <- scale( df.data$M3Mn )
  df.data$Al3 <- scale( df.data$M3Al )
  df.data$K3 <- scale( df.data$M3K )
  df.data$pH <- scale( df.data$pH_water )
  df.data$l.slope <- scale( log(df.data$slope ) )
  df.data$l.light <- scale( log(df.data$light) )
  df.data$l.dbh <- scale( log(df.data$dbh) )
  # Reduce the number of columns in the data set, as it just bloats the size of the model object if saved
  df.data <- df.data[, c("tag", "sp", "incgr", "l.slope", "l.dbh", "dbh", "l.light", "pH", "M3P", "Mn3", "Al3", "K3", "convex") ]

  if ( findBestRE ) {
    # For finding the optimal random effect
    listRE <- list( "~ l.dbh + l.light + Mn3 + M3P + K3 | sp",
                    "~ l.dbh + l.light + Mn3 + K3 | sp",
                    "~ l.dbh + l.light + Mn3 + M3P | sp",
                    "~ l.dbh + l.light + M3P + K3 | sp",
                    "~ l.dbh + l.light + Mn3 | sp",
                    "~ l.dbh + l.light + K3 | sp",
                    "~ l.dbh + l.light + M3P | sp",
                    "~ l.dbh + l.light | sp",
                    "~ l.light | sp",
                    "~ l.dbh | sp",
                    "~ 1 | sp")
    
    # Test each RE structure. REML must be used
    listMods <- vector( mode="list", length=length(listRE) )
    for ( j in 1:length(listRE) ) {
      startTime <- proc.time()
      randomEffect <- as.formula( listRE[[j]] )
      listMods[[j]] <- lme( fixed=fixedEffect, random=randomEffect,
                            data=df.data, method="REML", control=modCtrl,
                            weights=varExp() )
      print( paste( "Time for", listRE[[j]], ":" ) )
      print( proc.time() - startTime )
    }
    print( anova( listMods[[1]], listMods[[2]], listMods[[3]], listMods[[4]],
                  listMods[[5]], listMods[[6]], listMods[[7]], listMods[[8]],
                  listMods[[9]], listMods[[10]], listMods[[11]] ) )
    # For < 25 mm, #5 - dbh, light, Mn - is best
    # For 25 >= dbh < 250 mm, it's equivocal: #1 is the best by AIC, but #7
    # is the best by BIC. #7 - dbh, light, P - is much simpler
    # For >= 250 mm, #5 - dbh, light, Mn - is best
  } else {
    # These are the best RE structures
    listRE <- list( "~ l.dbh + l.light + Mn3 | sp",
                    "~ l.dbh + l.light + M3P | sp",
                    "~ l.dbh + l.light + Mn3 | sp")
    listSubsets <- list( "~ l.dbh & l.light & Mn3",
                         "~ l.dbh & l.light & M3P",
                         "~ l.dbh & l.light & Mn3" )
    randomEffect <- as.formula( listRE[[i]] )
    print( paste( "Random effect:", listRE[[i]] ) )
    
    startTime <- proc.time()
    # Run the (full) model
    mods[[i]] <- lme( fixed=fixedEffect, random=randomEffect, control=modCtrl, weights=varExp(),
                      data=df.data, method="ML" )
    print( summary( mods[[i]] ) )
    print( proc.time() - startTime )
    
    startTime <- proc.time()
    modAvg[[i]] <- dredge( mods[[i]], subset=as.formula(listSubsets[[i]]), 
                           beta="none", extra=c(r2.Xu), trace=T )
    print( proc.time() - startTime )

    ma <- model.avg(modAvg[[i]], cumsum(weight) <= .95, fit = F)    
    df.ci <- as.data.frame( confint( ma ) ) 
    df.ci <- cbind( fit=coef(ma), df.ci)
    df.ci <- df.ci[ order(row.names(df.ci)), ]
    write.csv( df.ci, paste0( dataPath, "BCI dredged growth model 95CI ", listSizes[[i]][[3]], ".csv" ) )
  }
}


####### Model averaging of the data without the BCI swamp
df <- GetGrowthModelData( noSwamp=T )  # From Utilities.R

# Model for the model averaging with the same RE structure as for the model with the swamp
fixedEffect <- log(incgr) ~ l.dbh + l.light + Mn3 + M3P + Al3 + K3 + pH + l.slope + convex
randomEffect <- ~ l.dbh + l.light + Mn3 + M3P | sp
modCtrl <- lmeControl( maxIter=150, msMaxIter=150, niterEM=50, msMaxEval=250 )
startTime <- proc.time()
mod.full.noswamp <- lme( fixed=fixedEffect, random=randomEffect,
                 data=df, method="ML", control=modCtrl,
                 weights=varExp() )
print( proc.time() - startTime )

# Model averaging
startTime <- proc.time()
noswamp.AIC <- dredge( mod.full.noswamp, subset=l.dbh & l.light & Mn3 & M3P, 
                       beta="none", extra=c(r2.Xu), trace=T )
print( proc.time() - startTime )

ma <- model.avg(noswamp.AIC, cumsum(weight) <= .95, fit = F)
df.ci <- as.data.frame( confint( ma ) ) 
df.ci <- cbind( fit=coef(ma), df.ci)
df.ci <- df.ci[ order(row.names(df.ci)), ]
write.csv( df.ci, file=paste0( dataPath, "BCI dredged growth model 95CI light,dbh,Mn,P-RE-no swamp-AIC.csv" ) )
x <- colMeans(df[,4:13])
write.csv( x, file=paste0( dataPath, "BCI growth model col means-no swamp.csv" ) )








