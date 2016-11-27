import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.inputFiles = 'csc_forsync.root'
options.outputFile = 'evttree.root'
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( SkipEvent =
cms.untracked.vstring('ProductNotFound') )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles )        

)
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet( limit = cms.untracked.int32(100) ),
        FwkJob = cms.untracked.PSet( limit = cms.untracked.int32(0) )
    ),
    categories = cms.untracked.vstring('FwkJob'),
    destinations = cms.untracked.vstring('cout')
)

process.load('Configuration.Geometry.GeometryExtended2023D1Reco_cff')
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

#call to customisation function cust_2023tilted imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023tilted(process)


process.newrh = cms.EDProducer("ME0StripRecHitProducer",
    me0DigiLabel = cms.InputTag("simMuonME0Digis"),
    recAlgo = cms.string('ME0RecHitStandardAlgo'),
    recAlgoConfig = cms.PSet(
        nEtaParts         = cms.int32(10),        
        nStrips           = cms.int32(768)    
    )
)


process.analyze = cms.EDAnalyzer("ME0SimHitAnalyzer",
        outFileName       = cms.untracked.string(options.outputFile)    
        )



process.p = cms.Path(process.analyze)

