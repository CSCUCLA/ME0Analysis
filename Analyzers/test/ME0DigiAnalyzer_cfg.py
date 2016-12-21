import FWCore.ParameterSet.Config as cms
import re

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
# no longer needed in CMSSW_9
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

#call to customisation function cust_2023tilted imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# process = cust_2023tilted(process)


def ME0ReReco(process, seq, name, nStrips, nPartitions, doMerge = True, onlyDigis = False) :
    process.load('Configuration.StandardSequences.Services_cff')
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    from SimMuon.GEMDigitizer.muonME0NewGeoDigis_cfi import simMuonME0NewGeoDigis
    
    digi_name          = name + "Digis"
    setattr(process.RandomNumberGeneratorService, digi_name, cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    ))
    
    setattr( process, digi_name, simMuonME0NewGeoDigis.clone(numberOfSrips=cms.uint32(nStrips), numberOfPartitions =cms.uint32(nPartitions), mergeDigis=cms.bool(doMerge)))
    seq += getattr(process, digi_name)    
    
    if onlyDigis : return
    
    rh_name          = name + "RecHits"
    setattr( process, rh_name, process.me0RecHits.clone(me0DigiLabel = cms.InputTag(digi_name)))
    seq += getattr(process, rh_name)
    
    seg_name          = name + "Segments"
    setattr( process, seg_name, process.me0Segments.clone(algo_type = 1,me0RecHitLabel = cms.InputTag(rh_name)))
    seq += getattr(process, seg_name)
    
    segM_name          = name + "SegmentMatching"
    setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag(seg_name)))
    seq += getattr(process, segM_name)
    
    me0Muon_name          = name + "Me0Muon"
    setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
    seq += getattr(process, me0Muon_name)
    
def doAnalysis(process, seq, name, nStrips, nPartitions, doMerge = True, onlyDigis = False) :
    ME0ReReco(process,seq,name,nStrips,nPartitions,doMerge,onlyDigis)
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0DigiAnalyzer",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string(name+"Digis"),
        runName           = cms.untracked.string(name+"_")
        )
    )
    seq += getattr(process, anName)


process.newdigiseq  = cms.Sequence()
doAnalysis(process,process.newdigiseq,"p8s768Merge",768,8,True,True)
doAnalysis(process,process.newdigiseq,"p8s640Merge",640,8,True,True)
doAnalysis(process,process.newdigiseq,"p8s512Merge",512,8,True,True)
doAnalysis(process,process.newdigiseq,"p8s384Merge",384,8,True,True)
doAnalysis(process,process.newdigiseq,"p8s256Merge",256,8,True,True)
doAnalysis(process,process.newdigiseq,"p8s128Merge",128,8,True,True)
# doAnalysis(process,process.newdigiseq,"p4s64Merge",64,4,True,True)
# doAnalysis(process,process.newdigiseq,"p4s4Merge",4,4,True,True)
# doAnalysis(process,process.newdigiseq,"p4s4NoMerge",4,4,False,True)
process.p = cms.Path(process.newdigiseq)

# process.newdigiseq1 = cms.Sequence()
# ME0ReReco(process,process.newdigiseq1,"p8s768Merge",768,8,True,True)
# process.newdigiseq2 = cms.Sequence()
# ME0ReReco(process,process.newdigiseq2,"p4s64Merge",64,4,True,True)
# process.newdigiseq3 = cms.Sequence()
# ME0ReReco(process,process.newdigiseq3,"p4s64NoMerge",64,4,False,True)
# 
# process.analyze = cms.EDAnalyzer("ME0DigiAnalyzer",
#         outFileName       = cms.untracked.string("p8s768Merge_" + options.outputFile),   
#         newDigiCollection = cms.untracked.string("p8s768MergeDigis"),
#         runName           = cms.untracked.string("p8s768Merge")
#         )

# process.newdigiseq = cms.Sequence()
# ME0ReReco(process,process.newdigiseq,"LessStrips",64,4)


# process.analyze = cms.EDAnalyzer("ME0DigiAnalyzer",
#         outFileName       = cms.untracked.string(options.outputFile),   
#         newDigiCollection = cms.untracked.string("simMuonME0NewGeoDigis"),
#         runName           = cms.untracked.string("")
#         )



# process.p = cms.Path(process.newdigiseq*process.analyze)

