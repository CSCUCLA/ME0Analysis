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
        fileNames = cms.untracked.vstring(options.inputFiles ),
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

process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

#call to customisation function cust_2023tilted imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# process = cust_2023tilted(process)
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
def ME0ReReco(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1),  doRUSegmentAlgo = True, onlyDigis = False, useDefault = False, doMerge = True) :

    from SimMuon.GEMDigitizer.muonME0ReDigis_cfi import simMuonME0ReDigis
    
    digi_name          = name + "Digis"
    setattr(process.RandomNumberGeneratorService, digi_name, cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    ))
    
    setattr( process, digi_name, simMuonME0ReDigis.clone(useBuiltinGeo = cms.bool(useDefault),mergeDigis = cms.bool(doMerge),numberOfSrips=cms.uint32(nStrips), 
                                                             numberOfPartitions =cms.uint32(nPartitions),
                                                             neutronAcceptance  =cms.double(neutBKGAcc),
                                                             layerReadout = layerRO,  minBXReadout       =cms.int32(-2),  maxBXReadout       =cms.int32(2)))
    seq += getattr(process, digi_name)    
    
    if onlyDigis : return
    
    rh_name          = name + "RecHits"
    setattr( process, rh_name, process.me0RecHits.clone(me0DigiLabel = cms.InputTag(digi_name)))
    seq += getattr(process, rh_name)
    
    seg_name          = name + "Segments"
    setattr( process, seg_name, process.me0Segments.clone(algo_type = (2 if doRUSegmentAlgo else 1),me0RecHitLabel = cms.InputTag(rh_name)))
    seq += getattr(process, seg_name)
    
    segM_name          = name + "SegmentMatching"
    setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag(seg_name)))
    seq += getattr(process, segM_name)
    
    me0Muon_name          = name + "Me0Muon"
    setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
    seq += getattr(process, me0Muon_name)
    
def doAnalysis(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1), doRUSegmentAlgo = True, onlyDigis = False, useDefault = False, doMerge = True) :
    anName = name + "Analysis"

    if name == "def" :
        setattr( process, anName, cms.EDAnalyzer("ME0DigiAnalyzer",
                                                 outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
                                                 newDigiCollection = cms.string("simMuonME0ReDigis"),
                                                 runName           = cms.untracked.string(name+"_")
                                                 )
                )
    else :
        ME0ReReco(process,seq,name,nStrips,nPartitions,neutBKGAcc,layerRO,doRUSegmentAlgo,onlyDigis,useDefault,doMerge)
        setattr( process, anName, cms.EDAnalyzer("ME0DigiAnalyzer",
                                                 outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
                                                 newDigiCollection = cms.string(name+"Digis"),
                                                 runName           = cms.untracked.string(name+"_")
                                                 )
                )
    seq += getattr(process, anName)


process.newdigiseq  = cms.Sequence()

doAnalysis(process,process.newdigiseq,"def",384,8,2.0,cms.vint32(1,1,1,1,1,1), True,True,False)
# doAnalysis(process,process.newdigiseq,"p8s384Def",384,8,2.0,cms.vint32(1,1,1,1,1,1), True,True,True)
# doAnalysis(process,process.newdigiseq,"p8s384DefNoMerge",384,8,2.0,cms.vint32(1,1,1,1,1,1), True,True,True,False)

process.p = cms.Path(process.newdigiseq)

