import FWCore.ParameterSet.Config as cms
import re

process = cms.Process("TEST")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( SkipEvent =
cms.untracked.vstring('ProductNotFound') )

process.source = cms.Source ("PoolSource",
        fileNames = cms.untracked.vstring(options.inputFiles )   ,
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")      

)
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet( limit = cms.untracked.int32(100) ),
        FwkJob = cms.untracked.PSet( limit = cms.untracked.int32(0) )
    ),
    categories = cms.untracked.vstring('FwkJob'),
    destinations = cms.untracked.vstring('cout')
)

process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')


process.analyze = cms.EDAnalyzer("ME0TrackDensity",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_std.root',options.outputFile)),
        tracks            = cms.InputTag("generalTracks")
        )
process.analyzeHP = cms.EDAnalyzer("ME0TrackDensity",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_hp.root',options.outputFile)),
        tracks            = cms.InputTag("probeTracks")
        )


import PhysicsTools.RecoAlgos.recoTrackSelector_cfi
process.probeTracks = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.probeTracks.quality = cms.vstring('highPurity')
process.probeTracks.tip = cms.double(3.5)
process.probeTracks.lip = cms.double(30.)
process.probeTracks.ptMin = cms.double(0.5)
process.probeTracks.minRapidity = cms.double(-3.0)
process.probeTracks.maxRapidity = cms.double(3.0)

process.p = cms.Path(process.probeTracks*process.analyze*process.analyzeHP)
