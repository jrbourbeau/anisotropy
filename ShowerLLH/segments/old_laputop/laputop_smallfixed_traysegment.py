#!/usr/bin/python

## With small showers, one cannot use timing/curvature as part of the likelihood. 
## If ONLY performing a charge LDF reconstruction, then only ONE step of Laputop is necessary.
## Note: the "Trigger" parameter must be set to 3 (since its default is 5).
## The parameters are borrowed from 1st step of the March 2012 parameters in "laputop_standard_traysegment.py".

from I3Tray import *

load("libgulliver")
load("liblilliput")
load("libtoprec")

@icetray.traysegment
def LaputopSmallFixed(tray, name, 
                       pulses='IceTopVEMPulses_0',
                       snowfactor=1.5,
                       If = lambda frame: True):
    ## Some more defaults
    fixcore = True     # do NOT keep core fixed
    fitsnow = False     # do NOT let the snow factor float


    ########## SERVICES FOR GULLIVER ##########

    ## The "simple lambda" snowservice
    tray.AddService("I3SimpleSnowCorrectionServiceFactory",name+"SimpleSnow")(
        ("Lambda", snowfactor)
        )


    ## This one is the standard one.
    tray.AddService("I3GulliverMinuitFactory",name+"Minuit")(
        ("MinuitPrintLevel",-2),  
        ("FlatnessCheck",True),  
        ("Algorithm","SIMPLEX"),  
        ("MaxIterations",2500),
        ("MinuitStrategy",2),
        ("Tolerance",0.01),
        )

    ## The Seed service
    tray.AddService("I3LaputopSeedServiceFactory",name+"ToprecSeed")(
        ("InCore", "ShowerLLH_proton"),
        ("InPlane", "ShowerPlane"),
        #("SnowCorrectionFactor", snowfactor),   # <-- snow correction factor here
        ("Beta", 2.6),                    # first guess for Beta
        ("InputPulses", pulses) # this'll let it first-guess at S125 automatically
        )

    ## Step 1:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam2")(
        ("FixCore", fixcore),
        ("FixTrackDir", True),   # Yes, fix the direction!
        ("FitSnowCorrection", fitsnow),
        ("IsBeta", True),
        ("MinBeta", 1.5),   ## From toprec... 2nd iteration (DLP, using beta)
        ("MaxBeta", 5.0),
        ("LimitCoreBoxSize", 200.0) 
    )

    tray.AddService("I3LaputopLikelihoodServiceFactory",name+"ToprecLike2")(
        ("DataReadout", pulses),
        ("BadStations", "IceTopExcludedStations"),
        ("DynamicCoreTreatment", 11.0),     # do the 11-meter core cut
        ("Trigger", 3),       ## Reduce min number of stations (the default is 5)
        ("Curvature","")    # NO timing likelihood (at first; will be overridden)
        )


    ################# GULLIVERIZED FITTER MODULE #######################

    ## This module performs the three steps
    tray.AddModule("I3LaputopFitter",name,
        SeedService = name+"ToprecSeed",
        NSteps = 1,     # <--- tells it how many services to look for and perform
        Parametrization1 = name+"ToprecParam2",   # the one parametrization
        StoragePolicy = "OnlyBestFit",
        Minimizer = name+"Minuit",
        LogLikelihoodService = name+"ToprecLike2",     # the one likelihood
        LDFFunctions = ["dlp"],
        CurvFunctions = [""],     # NO curvature or timing likelihood
        If = If,
        )


