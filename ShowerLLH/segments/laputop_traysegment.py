#!/usr/bin/python

## This tray segment will perform a three-step Laputop fit, using Tom Feusels' favorite parameters,
## presented at the Berkeley Collaboration Meeting in March 2012.

from I3Tray import *


load("libgulliver")
load("liblilliput")
load("libtoprec")
load("libparaboloid")

@icetray.traysegment
def Laputop(tray, name, 
                     pulses='CleanedHLCTankPulses',
                     excluded='ClusterCleaningExcludedStations',
                     fixcore  = False,
                     core_pos = 'ShowerCOG',
                     snowfactor=2.1):

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
        ("InCore", core_pos),
        ("InPlane", "ShowerPlane"),
#        ("SnowCorrectionFactor", snowfactor),   # <--- now obsolete
        ("Beta",2.6),                    # first guess for Beta
        ("InputPulses",pulses)  # this'll let it first-guess at S125 automatically
        )
    
    ## Step 1:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam2")(
        ("FixCore", fixcore),        
        ("FixTrackDir", True),
        ("IsBeta", True),
        ("MinBeta", 2.9),   ## From toprec... 2nd iteration (DLP, using beta)
        ("MaxBeta", 3.1),
        ("maxLogS125",8.0),        # Default is 6., be a bit safer, although should never happen to be this large
        ("VertexStepsize",10.0),   # The COG is very good for contained events, don't step too far
        ("LimitCoreBoxSize", 200.0) 
    )

    ## Step 2:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam3")(
        ("FixCore", fixcore),        
        ("FixTrackDir", False),      # FREE THE DIRECTION!
        ("IsBeta", True),
        ("MinBeta", 2.0),   ## From toprec... 3rd iteration (DLP, using beta)
        ("MaxBeta", 4.0),
        ("LimitCoreBoxSize", 15.0),
        ("maxLogS125",8.0),                   
        ## Use these smaller stepsizes instead of the defaults:
        ("VertexStepsize",5.0),      # default is 20
        ("SStepsize", 0.045),        # default is 1
        ("BetaStepsize",0.15)        # default is 0.6    
        )
        
    ## Step 3:
    tray.AddService("I3LaputopParametrizationServiceFactory",name+"ToprecParam4")(
        ("FixCore", fixcore),        
        ("FixTrackDir", True),
        ("IsBeta", True),
        ("MinBeta", 0.0),   
        ("MaxBeta", 10.0),
        ("LimitCoreBoxSize", 45.0),
        ("maxLogS125",8.0),        
        ## Use these smaller stepsizes instead of the defaults:
        ("VertexStepsize", 4.0),     # default is 20
        ("SStepsize", 0.045),        # default is 1
        ("BetaStepsize",0.15)        # default is 0.6 
        )
    
    tray.AddService("I3LaputopLikelihoodServiceFactory",name+"ToprecLike2")(
        ("DataReadout", pulses),
        ("BadStations", excluded),
        ("DynamicCoreTreatment", 11.0),     # do the 11-meter core cut (TF : using 5m, but didn't see a big difference in the plots)
        ("SaturationLikelihood", True),
        ("MaxIntraStationTimeDiff",80.0),    # Don't use time fluctuating tanks for timing fits, could really mess up the hard work
        ("Curvature",""),      # NO timing likelihood (at first; this will be overridden)
        ("SnowServiceName",name+"SimpleSnow")
        )
    
    
    ################# GULLIVERIZED FITTER MODULE #######################
    
    ## This module performs the three steps
    tray.AddModule("I3LaputopFitter",name)(
        ("SeedService",name+"ToprecSeed"),
        ("NSteps",3),            # <--- tells it how many services to look for and perform
        ("Parametrization1",name+"ToprecParam2"),   # the three parametrizations
        ("Parametrization2",name+"ToprecParam3"),
        ("Parametrization3",name+"ToprecParam4"),
        ("StoragePolicy","OnlyBestFit"),
        ("Minimizer",name+"Minuit"),
        ("LogLikelihoodService",name+"ToprecLike2"),     # the three likelihoods
        ("LDFFunctions",["dlp","dlp","dlp"]),
        ("CurvFunctions",["","gausspar","gausspar"])   # VERY IMPORTANT : use time Llh for step 3, but fix direction!
        )
    
    '''
    # THIS IS RECOMMENDED
    tray.AddModule("I3ParaboloidFitter","paraboloid")(
        ("SeedService",name+"ToprecSeed"),
        ("LogLikelihood",name+"ToprecLike2"),
        ("VertexStepSize",5.0*I3Units.m),
        ("MaxMissingGridPoints",1),
        ("GridpointVertexCorrection","seedprep"),
        ("Minimizer",name+"minuit"),
    )
    '''
