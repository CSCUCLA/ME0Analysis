import FWCore.ParameterSet.Config as cms
import re

process = cms.Process("TEST")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( SkipEvent =
cms.untracked.vstring('ProductNotFound') )

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(options.inputFiles )        ,
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet( limit = cms.untracked.int32(-1) ),
        FwkJob = cms.untracked.PSet( limit = cms.untracked.int32(0) )
    ),
    categories = cms.untracked.vstring('FwkJob'),
    destinations = cms.untracked.vstring('cout')
)

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

#call to customisation function cust_2023tilted imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# process = cust_2023tilted(process)
# process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
def ME0ReReco(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1), minSegmentLayers = 4, usePads = False) :
    from SimMuon.GEMDigitizer.muonME0ReDigis_cfi import simMuonME0ReDigis
    
    digi_name          = name + "Digis"
    setattr(process.RandomNumberGeneratorService, digi_name, cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    ))
    
    setattr( process, digi_name, simMuonME0ReDigis.clone(useBuiltinGeo = cms.bool(False),numberOfStrips=cms.uint32(nStrips), 
                                                            numberOfPartitions =cms.uint32(nPartitions),
                                                            neutronAcceptance  =cms.double(neutBKGAcc),
                                                            layerReadout = layerRO,usePadsForDefaultGeo =cms.bool(usePads)))
    seq += getattr(process, digi_name)    
    
    
    rh_name          = name + "RecHits"
    setattr( process, rh_name, process.me0RecHits.clone(me0DigiLabel = cms.InputTag(digi_name)))
    seq += getattr(process, rh_name)
    
    seg_name          = name + "Segments"
    setattr( process, seg_name, cms.EDProducer("ME0SegmentProducer",
    algo_psets = cms.VPSet(cms.PSet(
        algo_name = cms.string('ME0SegmentAlgorithm'),
        algo_pset = cms.PSet(
            ME0Debug = cms.untracked.bool(True),            
            dEtaChainBoxMax = cms.double(0.05),
            dPhiChainBoxMax = cms.double(0.02),
            dTimeChainBoxMax = cms.double(1.5),
            dXclusBoxMax = cms.double(1.0),
            dYclusBoxMax = cms.double(5.0),
            maxRecHitsInCluster = cms.int32(6),
            minHitsPerSegment = cms.uint32(3),
            preClustering = cms.bool(True),
            preClusteringUseChaining = cms.bool(True)
        )
    ), 
        cms.PSet(
            algo_name = cms.string('ME0SegAlgoRU'),
            algo_pset = cms.PSet(
                allowWideSegments = cms.bool(True),
                doCollisions = cms.bool(True),
                requireCentralBX = cms.bool(True),
                maxChi2Additional = cms.double(100.0),
                maxChi2Prune = cms.double(50),
                maxChi2GoodSeg = cms.double(50),
                maxPhiSeeds = cms.double(1.2*0.35/nStrips),
                maxPhiAdditional = cms.double(1.2*0.35/nStrips),
                maxETASeeds = cms.double(0.8/nPartitions),
                maxTOFDiff = cms.double(25),
                minNumberOfHits = cms.uint32(minSegmentLayers),
            )
        )),
    algo_type = cms.int32(2),
    me0RecHitLabel = cms.InputTag(rh_name)
    ))
    seq += getattr(process, seg_name)
    
    
def doAnalysis(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1), minSegmentLayers = 4, usePads = False) :
    ME0ReReco(process,seq,name,nStrips,nPartitions,neutBKGAcc,layerRO,minSegmentLayers,usePads)
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0SegmentTreeMaker",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string(name+"Digis"),
        segmentCollection = cms.string(name+"Segments"),
        recHitCollection = cms.string(name+"RecHits"),
        runName           = cms.untracked.string(name+"_")
        )        
    )
    seq += getattr(process, anName)
    print re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)


def doAnalysisOnly(process, seq, name) :
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0SegmentTreeMaker",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string("simMuonME0ReDigis"),
        segmentCollection = cms.string("me0Segments"),
        recHitCollection = cms.string("me0RecHits"),
        runName           = cms.untracked.string(name+"_"),
        )
    )
    seq += getattr(process, anName)

process.newdigiseq  = cms.Sequence()

# doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),4)
# doAnalysis(process,process.newdigiseq,"p8s192" ,192,8 ,2.0,cms.vint32(1,1,1,1,1,1),4,True)
doAnalysisOnly(process,process.newdigiseq,"p8s384")


process.p = cms.Path(process.newdigiseq)