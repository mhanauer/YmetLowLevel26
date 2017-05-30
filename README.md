# Jags-Ymet-XmetSsubj-MrobustHier.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="X" , x2Name = "X2", x3Name = "X3", x4Name = "X4", x5Name = "X5", x6Name = "X6", x7Name = "X7", x8Name = "X8", x9Name = "X9", x10Name = "X10", x11Name = "X11", x12Name = "X12", x13Name="X13" , x14Name = "X14", x15Name = "X15", x16Name = "X16", x17Name = "X17", x18Name = "X18", x19Name = "X19", x20Name = "X20", x21Name = "X21", x22Name = "X22", x23Name = "X23", x24Name = "X24", x25Name = "25", x26Name = "X26", yName="y" , sName="s" ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault) { 
  
  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  x2 = data[,x2Name]
  x3 = data[,x3Name]
  x4 = data[,x4Name]
  x5 = data[,x5Name]
  x6 = data[,x6Name]
  x7 = data[,x7Name]
  x8 = data[,x8Name]
  x9 = data[,x9Name]
  x10 = data[,x10Name]
  x11 = data[,x11Name]
  x12 = data[,x12Name]
  x13 = data[,x13Name]
  x14 = data[,x14Name]
  x15 = data[,x15Name]
  x16 = data[,x16Name]
  x17 = data[,x17Name]
  x18 = data[,x18Name]
  x19 = data[,x19Name]
  x20 = data[,x20Name]
  x21 = data[,x21Name]
  x22 = data[,x22Name]
  x23 = data[,x23Name]
  x24 = data[,x24Name]
  x25 = data[,x25Name]
  x26 = data[,x26Name]
  # Convert sName to consecutive integers:
  s = as.numeric(factor(data[,sName]))
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    x2 = x2,
    x3 = x3,
    x4 = x4,
    x5 = x5,
    x6 = x6,
    x7 = x7,
    x8 = x8,
    x9 = x9,
    x10 = x10,
    x11 = x11,
    x12 = x12,
    x13 = x13,
    x14 = x14,
    x15 = x15,
    x16 = x16,
    x17 = x17,
    x18 = x18,
    x19 = x19,
    x20 = x20,
    x21 = x21,
    x22 = x22,
    x23 = x23,
    x24 = x24,
    x25 = x25,
    x26 = x26,
    y = y ,
    s = s ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
  Ntotal <- length(y)
  xm <- mean(x)
  x2m <- mean(x2)
  x3m <- mean(x3)
  x4m <- mean(x4)
  x5m <- mean(x5)
  x6m <- mean(x6)
  x7m <- mean(x7)
  x8m <- mean(x8)
  x9m <- mean(x9)
  x10m <- mean(x10)
  x11m <- mean(x11)
  x12m <- mean(x12)
  x13m <- mean(x13)
  x14m <- mean(x14)
  x15m <- mean(x15)
  x16m <- mean(x16)
  x17m <- mean(x17)
  x18m <- mean(x18)
  x19m <- mean(x19)
  x20m <- mean(x20)
  x21m <- mean(x21)
  x22m <- mean(x22)
  x23m <- mean(x23)
  x24m <- mean(x24)
  x25m <- mean(x25)
  x26m <- mean(x26)
  ym <- mean(y)
  xsd <- sd(x)
  x2sd <- sd(x2)
  x3sd <- sd(x3)
  x4sd <- sd(x4)
  x5sd <- sd(x5)
  x6sd <- sd(x6)
  x7sd <- sd(x7)
  x8sd <- sd(x8)
  x9sd <- sd(x9)
  x10sd <- sd(x10)
  x11sd <- sd(x11)
  x12sd <- sd(x12)
  x13sd <- sd(x13)
  x14sd <- sd(x14)
  x15sd <- sd(x15)
  x16sd <- sd(x16)
  x17sd <- sd(x17)
  x18sd <- sd(x18)
  x19sd <- sd(x19)
  x20sd <- sd(x20)
  x21sd <- sd(x21)
  x22sd <- sd(x22)
  x23sd <- sd(x23)
  x24sd <- sd(x24)
  x25sd <- sd(x25)
  x26sd <- sd(x26)
  
  ysd <- sd(y)
  for ( i in 1:length(y) ) {
  zx[i] <- ( x[i] - xm ) / xsd
  zx2[i] <- ( x2[i] - x2m ) / x2sd
  zx3[i] <- ( x3[i] - x3m ) / x3sd
  zx4[i] <- ( x4[i] - x4m ) / x4sd
  zx5[i] <- ( x5[i] - x5m ) / x5sd
  zx6[i] <- ( x6[i] - x6m ) / x6sd
  zx7[i] <- ( x7[i] - x7m ) / x7sd
  zx8[i] <- ( x8[i] - x8m ) / x8sd
  zx9[i] <- ( x9[i] - x9m ) / x9sd
  zx10[i] <- ( x10[i] - x10m ) / x10sd
  zx11[i] <- ( x11[i] - x11m ) / x11sd
  zx12[i] <- ( x12[i] - x12m ) / x12sd
  zx13[i] <- ( x13[i] - x13m ) / x13sd
  zx14[i] <- ( x14[i] - x14m ) / x14sd
  zx15[i] <- ( x15[i] - x15m ) / x15sd
  zx16[i] <- ( x16[i] - x16m ) / x16sd
  zx17[i] <- ( x17[i] - x17m ) / x17sd
  zx18[i] <- ( x18[i] - x18m ) / x18sd
  zx19[i] <- ( x19[i] - x19m ) / x19sd
  zx20[i] <- ( x20[i] - x20m ) / x20sd
  zx21[i] <- ( x21[i] - x21m ) / x21sd
  zx22[i] <- ( x22[i] - x22m ) / x22sd
  zx23[i] <- ( x23[i] - x23m ) / x23sd
  zx24[i] <- ( x24[i] - x24m ) / x24sd
  zx25[i] <- ( x25[i] - x25m ) / x25sd
  zx26[i] <- ( x26[i] - x26m ) / x26sd
  zy[i] <- ( y[i] - ym ) / ysd
  }
  }
  # Specify the model for standardized data:
  model {
  for ( i in 1:Ntotal ) {
  zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx2[i] + zbeta3[s[i]] * zx3[i] + zbeta4[s[i]] * zx4[i]  + zbeta5[s[i]] * zx5[i] + zbeta6[s[i]] * zx6[i] + zbeta7[s[i]] * zx7[i] + zbeta8[s[i]] * zx8[i] + zbeta9[s[i]] * zx9[i] + zbeta10[s[i]] * zx10[i] + zbeta11[s[i]] * zx11[i] + zbeta12[s[i]] * zx12[i] + zbeta13[s[i]] * zx13[i] + zbeta14[s[i]] * zx14[i] + zbeta15[s[i]] * zx15[i] + zbeta16[s[i]] * zx16[i]  + zbeta17[s[i]] * zx17[i] + zbeta18[s[i]] * zx18[i] + zbeta19[s[i]] * zx19[i] + zbeta20[s[i]] * zx20[i] + zbeta21[s[i]] * zx21[i] + zbeta22[s[i]] * zx22[i] + zbeta23[s[i]] * zx23[i] + zbeta24[s[i]] * zx24[i] + zbeta25[s[i]] * zx25[i] + zbeta26[s[i]] * zx26[i], 1/zsigma^2 , nu )
  }
  for ( j in 1:Nsubj ) {
  zbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )  
  zbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )
  zbeta2[j] ~ dnorm( zbeta2mu , 1/(zbeta2sigma)^2 )
  zbeta3[j] ~ dnorm( zbeta3mu , 1/(zbeta3sigma)^2 )
  zbeta4[j] ~ dnorm( zbeta4mu , 1/(zbeta4sigma)^2 )
  zbeta5[j] ~ dnorm( zbeta5mu , 1/(zbeta5sigma)^2 )
  zbeta6[j] ~ dnorm( zbeta6mu , 1/(zbeta6sigma)^2 )
  zbeta7[j] ~ dnorm( zbeta7mu , 1/(zbeta7sigma)^2 )
  zbeta8[j] ~ dnorm( zbeta8mu , 1/(zbeta8sigma)^2 )
  zbeta9[j] ~ dnorm( zbeta9mu , 1/(zbeta9sigma)^2 )
  zbeta10[j] ~ dnorm( zbeta10mu , 1/(zbeta10sigma)^2 )
  zbeta11[j] ~ dnorm( zbeta11mu , 1/(zbeta11sigma)^2 )
  zbeta12[j] ~ dnorm( zbeta12mu , 1/(zbeta12sigma)^2 )
  zbeta13[j] ~ dnorm( zbeta13mu , 1/(zbeta13sigma)^2 )
  zbeta14[j] ~ dnorm( zbeta14mu , 1/(zbeta14sigma)^2 )
  zbeta15[j] ~ dnorm( zbeta15mu , 1/(zbeta15sigma)^2 )
  zbeta16[j] ~ dnorm( zbeta16mu , 1/(zbeta16sigma)^2 )
  zbeta17[j] ~ dnorm( zbeta17mu , 1/(zbeta17sigma)^2 )
  zbeta18[j] ~ dnorm( zbeta18mu , 1/(zbeta18sigma)^2 )
  zbeta19[j] ~ dnorm( zbeta19mu , 1/(zbeta19sigma)^2 )
  zbeta20[j] ~ dnorm( zbeta20mu , 1/(zbeta20sigma)^2 )
  zbeta21[j] ~ dnorm( zbeta21mu , 1/(zbeta21sigma)^2 )
  zbeta22[j] ~ dnorm( zbeta22mu , 1/(zbeta22sigma)^2 )
  zbeta23[j] ~ dnorm( zbeta23mu , 1/(zbeta23sigma)^2 )
  zbeta24[j] ~ dnorm( zbeta24mu , 1/(zbeta24sigma)^2 )
  zbeta25[j] ~ dnorm( zbeta25mu , 1/(zbeta25sigma)^2 )
  zbeta26[j] ~ dnorm( zbeta26mu , 1/(zbeta26sigma)^2 )
  }
  # Priors vague on standardized scale:
  zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta2mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta3mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta4mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta5mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta6mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta7mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta8mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta9mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta10mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta11mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta12mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta13mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta14mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta15mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta16mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta17mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta18mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta19mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta20mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta21mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta22mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta23mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta24mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta25mu ~ dnorm( 0 , 1/(10)^2 )
  zbeta26mu ~ dnorm( 0 , 1/(10)^2 )


  zsigma ~ dnorm( 1.0E-3 , 1.0E+3 )
  zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta2sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta3sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta4sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta5sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta6sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta7sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta8sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta9sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta10sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta11sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta12sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta13sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta14sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta15sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta16sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta17sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta18sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta19sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta20sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta21sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta22sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta23sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta24sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta25sigma ~ dunif( 1.0E-3 , 1.0E+3 )
  zbeta26sigma ~ dunif( 1.0E-3 , 1.0E+3 )


  nu ~ dexp(1/30.0)
  # Transform to original scale:
  for ( j in 1:Nsubj ) {
  beta1[j] <- zbeta1[j] * ysd / xsd 
  beta2[j] <- zbeta2[j] * ysd / x2sd
  beta3[j] <- zbeta3[j] * ysd / x3sd
  beta4[j] <- zbeta4[j] * ysd / x4sd
  beta5[j] <- zbeta5[j] * ysd / x5sd
  beta6[j] <- zbeta6[j] * ysd / x6sd
  beta7[j] <- zbeta7[j] * ysd / x7sd
  beta8[j] <- zbeta8[j] * ysd / x8sd
  beta9[j] <- zbeta9[j] * ysd / x9sd
  beta10[j] <- zbeta10[j] * ysd / x10sd
  beta11[j] <- zbeta11[j] * ysd / x11sd
  beta12[j] <- zbeta12[j] * ysd / x12sd
  beta13[j] <- zbeta13[j] * ysd / x13sd 
  beta14[j] <- zbeta14[j] * ysd / x14sd
  beta15[j] <- zbeta15[j] * ysd / x15sd
  beta16[j] <- zbeta16[j] * ysd / x16sd
  beta17[j] <- zbeta17[j] * ysd / x17sd
  beta18[j] <- zbeta18[j] * ysd / x18sd
  beta19[j] <- zbeta19[j] * ysd / x19sd
  beta20[j] <- zbeta20[j] * ysd / x20sd
  beta21[j] <- zbeta22[j] * ysd / x22sd
  beta22[j] <- zbeta22[j] * ysd / x22sd
  beta23[j] <- zbeta23[j] * ysd / x23sd
  beta24[j] <- zbeta24[j] * ysd / x24sd
  beta25[j] <- zbeta25[j] * ysd / x25sd
  beta26[j] <- zbeta26[j] * ysd / x26sd
  beta0[j] <- zbeta0[j] * ysd  + ym - zbeta1[j] * xm * ysd / xsd + zbeta2[j] * x2m * ysd / x2sd + zbeta3[j] * x3m * ysd / x3sd + zbeta4[j] * x4m * ysd / x4sd + zbeta5[j] * x5m * ysd / x5sd + zbeta6[j] * x6m * ysd / x6sd + zbeta7[j] * x7m * ysd / x7sd + zbeta8[j] * x8m * ysd / x8sd + zbeta9[j] * x9m * ysd / x9sd + zbeta10[j] * x10m * ysd / x10sd + zbeta11[j] * x11m * ysd / x11sd + zbeta12[j] * x12m * ysd / x12sd + ym - zbeta13[j] * x13m * ysd / x13sd + zbeta14[j] * x14m * ysd / x14sd + zbeta15[j] * x15m * ysd / x15sd + zbeta16[j] * x16m * ysd / x16sd + zbeta17[j] * x17m * ysd / x17sd + zbeta18[j] * x18m * ysd / x18sd + zbeta19[j] * x19m * ysd / x19sd + zbeta20[j] * x20m * ysd / x20sd + zbeta21[j] * x21m * ysd / x21sd + zbeta22[j] * x22m * ysd / x22sd + zbeta23[j] * x23m * ysd / x23sd + zbeta24[j] * x24m * ysd / x24sd + zbeta25[j] * x25m * ysd / x25sd + zbeta26[j] * x26m * ysd / x24sd
  }
  beta1mu <- zbeta1mu * ysd / xsd
  beta2mu <- zbeta2mu * ysd / x2sd
  beta3mu <- zbeta3mu * ysd / x3sd
  beta4mu <- zbeta4mu * ysd / x4sd
  beta5mu <- zbeta5mu * ysd / x5sd
  beta6mu <- zbeta6mu * ysd / x6sd
  beta7mu <- zbeta7mu * ysd / x7sd
  beta8mu <- zbeta8mu * ysd / x8sd
  beta9mu <- zbeta9mu * ysd / x9sd
  beta10mu <- zbeta10mu * ysd / x10sd
  beta11mu <- zbeta11mu * ysd / x11sd
  beta12mu <- zbeta12mu * ysd / x12sd
  beta13mu <- zbeta13mu * ysd / x13sd
  beta14mu <- zbeta14mu * ysd / x14sd
  beta15mu <- zbeta15mu * ysd / x15sd
  beta16mu <- zbeta16mu * ysd / x16sd
  beta17mu <- zbeta17mu * ysd / x17sd
  beta18mu <- zbeta18mu * ysd / x18sd
  beta19mu <- zbeta19mu * ysd / x19sd
  beta20mu <- zbeta20mu * ysd / x20sd
  beta21mu <- zbeta21mu * ysd / x21sd
  beta22mu <- zbeta22mu * ysd / x22sd
  beta23mu <- zbeta23mu * ysd / x23sd
  beta24mu <- zbeta24mu * ysd / x24sd
  beta25mu <- zbeta25mu * ysd / x25sd
  beta26mu <- zbeta26mu * ysd / x26sd


  beta0mu <- zbeta0mu * ysd  + ym - zbeta1mu * xm * ysd / xsd + zbeta2mu * x2m * ysd / x2sd + zbeta3mu * x3m * ysd / x3sd + zbeta4mu * x4m * ysd / x4sd + zbeta5mu * x5m * ysd / x5sd + zbeta6mu * x6m * ysd / x6sd + zbeta7mu * x7m * ysd / x7sd + zbeta8mu * x8m * ysd / x8sd + zbeta9mu * x9m * ysd / x9sd + zbeta10mu * x10m * ysd / x10sd + zbeta11mu * x11m * ysd / x11sd + zbeta12mu * x12m * ysd / x12sd + zbeta12mu * x12m * ysd / x12sd + zbeta13mu * x13m * ysd / x13sd + zbeta14mu * x14m * ysd / x14sd + zbeta15mu * x15m * ysd / x15sd + zbeta16mu * x16m * ysd / x16sd + zbeta17mu * x17m * ysd / x17sd + zbeta18mu * x18m * ysd / x18sd + zbeta19mu * x19m * ysd / x19sd + zbeta20mu * x20m * ysd / x20sd + zbeta21mu * x21m * ysd / x21sd + zbeta22mu * x22m * ysd / x22sd + zbeta23mu * x23m * ysd / x23sd + zbeta24mu * x24m * ysd / x24sd+ zbeta25mu * x25m * ysd / x25sd+ zbeta26mu * x26m * ysd / x26sd
  sigma <- zsigma * ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1","beta2" , "beta3",  "beta4", "beta5", "beta6",  "beta7",  "beta8",  "beta9",  "beta10",  "beta11",  "beta12", "beta13" ,  "beta14","beta15" , "beta16",  "beta17", "beta18", "beta19",  "beta20",  "beta21",  "beta22",  "beta23",  "beta24",  "beta25", "beta26", "beta0mu" , "beta1mu" , "beta2mu", "beta3mu", "beta4mu", "beta5mu", "beta6mu", "beta7mu", "beta8mu", "beta9mu", "beta10mu", "beta11mu", "beta12mu","beta13mu" , "beta14mu" , "beta15mu", "beta16mu", "beta17mu", "beta18mu", "beta19mu", "beta20mu", "beta21mu", "beta22mu", "beta23mu", "beta24mu", "beta25mu", "beta26mu", 
                  "zbeta0" , "zbeta1" , "zbeta2", "zbeta3", "zbeta4", "zbeta5", "zbeta6", "zbeta7", "zbeta8", "zbeta9", "zbeta10", "zbeta11", "zbeta12", "zbeta13" , "zbeta14" , "zbeta15", "zbeta16", "zbeta17", "zbeta18", "zbeta19", "zbeta20", "zbeta21", "zbeta22", "zbeta23", "zbeta24", "zbeta25", "zbeta26",  "zbeta0mu" , "zbeta1mu" ,"zbeta2mu", "zbeta3mu", "zbeta4mu", "zbeta5mu", "zbeta6mu", "zbeta7mu", "zbeta8mu", "zbeta9mu", "zbeta10mu", "zbeta11mu", "zbeta12mu", "zbeta13mu" , "zbeta14mu" ,"zbeta15mu", "zbeta16mu", "zbeta17mu", "zbeta18mu", "zbeta19mu", "zbeta20mu", "zbeta21mu", "zbeta22mu", "zbeta23mu", "zbeta24mu", "zbeta25mu","zbeta26mu",
                  "zsigma", "sigma", "nu" , "zbeta0sigma" , "zbeta1sigma","zbeta2sigma", "zbeta3sigma", "zbeta4sigma" , "zbeta5sigma" , "zbeta6sigma" , "zbeta7sigma" , "zbeta8sigma" , "zbeta9sigma" , "zbeta10sigma" , "zbeta11sigma" , "zbeta12sigma", "zbeta13sigma" , "zbeta14sigma","zbeta15sigma", "zbeta16sigma", "zbeta17sigma" , "zbeta18sigma" , "zbeta19sigma" , "zbeta20sigma" , "zbeta21sigma" , "zbeta22sigma" , "zbeta23sigma" , "zbeta24sigma" , "zbeta25sigma", "zbeta26sigma")
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 2000
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}
