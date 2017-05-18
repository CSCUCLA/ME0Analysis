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

process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023tilted

#call to customisation function cust_2023tilted imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
# no longer needed in CMSSW_9
# process = cust_2023tilted(process)
# process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
def ME0ReReco(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1),  doRUSegmentAlgo = True, minSegmentLayers = 4,timeWindow = 25, onlyDigis = False) :
    from SimMuon.GEMDigitizer.muonME0NewGeoDigis_cfi import simMuonME0NewGeoDigis
    
    digi_name          = name + "Digis"
    setattr(process.RandomNumberGeneratorService, digi_name, cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    ))
    
    setattr( process, digi_name, simMuonME0NewGeoDigis.clone(numberOfSrips=cms.uint32(nStrips), 
                                                             numberOfPartitions =cms.uint32(nPartitions),
                                                             neutronAcceptance  =cms.double(neutBKGAcc),
                                                             layerReadout = layerRO))
    seq += getattr(process, digi_name)    
    
    if onlyDigis : return
    
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
                maxChi2Additional = cms.double(100.0),
                maxChi2Prune = cms.double(50),
                maxChi2GoodSeg = cms.double(50),
                maxPhiSeeds = cms.double(1.2*0.35/nStrips),
                maxPhiAdditional = cms.double(1.2*0.35/nStrips),
                maxETASeeds = cms.double(0.8/nPartitions),
                maxTOFDiff = cms.double(timeWindow),
                minNumberOfHits = cms.uint32(minSegmentLayers),
            )
        )),
    algo_type = cms.int32(2 if doRUSegmentAlgo else 1),
    me0RecHitLabel = cms.InputTag(rh_name)
    ))
    seq += getattr(process, seg_name)
    
#     segM_name          = name + "SegmentMatching"
#     setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag(seg_name)))
#     seq += getattr(process, segM_name)
#      
#     me0Muon_name          = name + "Me0Muon"
#     setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
#     seq += getattr(process, me0Muon_name)
    
def doAnalysis(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1), doRUSegmentAlgo = True, minSegmentLayers = 4,timeWindow = 25, onlyDigis = False) :
    ME0ReReco(process,seq,name,nStrips,nPartitions,neutBKGAcc,layerRO,doRUSegmentAlgo,minSegmentLayers,timeWindow,onlyDigis)
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0SegmentAnalyzer",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string(name+"Digis"),
        segmentCollection = cms.string(name+"Segments"),
        recHitCollection = cms.string(name+"RecHits"),
        runName           = cms.untracked.string(name+"_")
        )
    )
    seq += getattr(process, anName)


process.newdigiseq  = cms.Sequence()
# doAnalysis(process,process.newdigiseq,"p128s768Merge",768,128,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p64s768Merge",768,64,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p64s4Merge",4,64,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p4s4Merge",4,4,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p4s768Merge",768,4,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p8s128Merge",128,8,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p8s256Merge",256,8,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p8s512Merge",512,8,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p4s512Merge",512,4,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p8s1024Merge",1024,8,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p16s768Merge",768,16,0.0,cms.vint32(1,1,1,1,1,1), True,True)
# doAnalysis(process,process.newdigiseq,"p8s640Merge",640,8,True,True)
# doAnalysis(process,process.newdigiseq,"p8s512Merge",512,8,True,True)
# doAnalysis(process,process.newdigiseq,"p8s384Merge",384,8,True,True)

# doAnalysis(process,process.newdigiseq,"p12s768",768,12,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p12s512",512,12,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p12s384",384,12,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p12s256",128,12,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p12s128",128,12,2.0,cms.vint32(1,1,1,1,1,1),True)

# doAnalysis(process,process.newdigiseq,"p8s768",768,8,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p8s512",512,8,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p8s384",384,8,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p8s128",128,8,2.0,cms.vint32(1,1,1,1,1,1),True)
# 
# doAnalysis(process,process.newdigiseq,"p6s768",768,6,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512",512,6,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s384",384,6,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s128",128,6,2.0,cms.vint32(1,1,1,1,1,1),True)
# 
# doAnalysis(process,process.newdigiseq,"p4s768",768,4,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p4s512",512,4,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p4s384",384,4,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p4s128",128,4,2.0,cms.vint32(1,1,1,1,1,1),True)

# doAnalysis(process,process.newdigiseq,"p8s256Merge",256,8,True,True)
# doAnalysis(process,process.newdigiseq,"p8s128Merge",128,8,True,True)
# doAnalysis(process,process.newdigiseq,"p4s64Merge",64,4,True,True)
# doAnalysis(process,process.newdigiseq,"p4s4Merge",4,4,True,True)
# doAnalysis(process,process.newdigiseq,"p4s4NoMerge",4,4,False,True)

# doAnalysis(process,process.newdigiseq,"p12s256",128,12,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p8s768" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512" ,512,6 ,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p4s768" ,768,4 ,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),True)
# 
# 
# doAnalysis(process,process.newdigiseq,"p6s512M1" ,512,6 ,2.0,cms.vint32(0,1,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512M2" ,512,6 ,2.0,cms.vint32(1,0,1,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512M3" ,512,6 ,2.0,cms.vint32(1,1,0,1,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512M4" ,512,6 ,2.0,cms.vint32(1,1,1,0,1,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512M5" ,512,6 ,2.0,cms.vint32(1,1,1,1,0,1),True)
# doAnalysis(process,process.newdigiseq,"p6s512M6" ,512,6 ,2.0,cms.vint32(1,1,1,1,1,0),True)


#std
# doAnalysis(process,process.newdigiseq,"p12s256",128,12,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s768" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s512" ,512,6 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p4s768" ,768,4 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# 
# #test removeing layer, looser layer reqs
# doAnalysis(process,process.newdigiseq,"p8s768M6" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p8s768M6L3" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p8s768L3" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,3)
# 
# #do fine grain...comment out already processed
# doAnalysis(process,process.newdigiseq,"p2s768"  ,768,2 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# #doAnalysis(process,process.newdigiseq,"p4s768"  ,768,4 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s768"  ,768,6 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# #doAnalysis(process,process.newdigiseq,"p8s768"  ,768,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p10s768" ,768,10 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# 
# doAnalysis(process,process.newdigiseq,"p8s256" ,256,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s512" ,512,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s640" ,640,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)
# #doAnalysis(process,process.newdigiseq,"p8s768" ,768,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4)

#neutron stress test
#for full 112.5 sample
##FOR TDR

# doAnalysis(process,process.newdigiseq,"p8s38475NL4"   ,384,8 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s38437p5NL4" ,384,8 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s38422p5NL4" ,384,8 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s38410NL4"   ,384,8 ,0.1333,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s3847p5NL4"  ,384,8 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s38475NL5"   ,384,8 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s38437p5NL5" ,384,8 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s38422p5NL5" ,384,8 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s38410NL5"   ,384,8 ,0.1333,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s3847p5NL5"  ,384,8 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,5)




# doAnalysis(process,process.newdigiseq,"p6s51275NL3"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL3" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL3" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NL3"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL3"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51275NL4"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL4" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL4" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NL4"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL4"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51275NL5"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL5" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL5" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51210NL5"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL5"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51275NL6"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL6" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL6" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s51210NL6"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL6"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s51275NM6L3"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM6L3" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM6L3" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L3"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L3"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51275NM6L4"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM6L4" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM6L4" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L4"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L4"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51275NM6L5"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM6L5" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM6L5" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L5"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L5"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51275NM56L3"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM56L3" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM56L3" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NM56L3"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM56L3"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51275NM56L4"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM56L4" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM56L4" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NM56L4"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM56L4"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51275NM456L3"   ,512,6 ,2.0   ,cms.vint32(1,1,1,0,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM456L3" ,512,6 ,0.5   ,cms.vint32(1,1,1,0,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM456L3" ,512,6 ,0.3   ,cms.vint32(1,1,1,0,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NM456L3"   ,512,6 ,0.1333,cms.vint32(1,1,1,0,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM456L3"  ,512,6 ,0.1   ,cms.vint32(1,1,1,0,0,0),True,3)


# doAnalysis(process,process.newdigiseq,"p6s51275NL4T0"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL4T0" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL4T0" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51210NL4T0"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL4T0"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51275NL5T0"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,5,0)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NL5T0" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,5,0)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NL5T0" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,5,0)
# doAnalysis(process,process.newdigiseq,"p6s51210NL5T0"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,1),True,5,0)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL5T0"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,5,0)
# doAnalysis(process,process.newdigiseq,"p6s51275NM6L3T0"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,0),True,3,0)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM6L3T0" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,0),True,3,0)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM6L3T0" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,0),True,3,0)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L3T0"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,0),True,3,0)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L3T0"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,0),True,3,0)
# doAnalysis(process,process.newdigiseq,"p6s51275NM6L4T0"   ,512,6 ,2.0   ,cms.vint32(1,1,1,1,1,0),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51237p5NM6L4T0" ,512,6 ,0.5   ,cms.vint32(1,1,1,1,1,0),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51222p5NM6L4T0" ,512,6 ,0.3   ,cms.vint32(1,1,1,1,1,0),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L4T0"   ,512,6 ,0.1333,cms.vint32(1,1,1,1,1,0),True,4,0)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L4T0"  ,512,6 ,0.1   ,cms.vint32(1,1,1,1,1,0),True,4,0)

# doAnalysis(process,process.newdigiseq,"p12s51275NL4"   ,512,12 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p12s51237p5NL4" ,512,12 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p12s51222p5NL4" ,512,12 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p12s51210NL4"   ,512,12 ,0.1333,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p12s5127p5NL4"  ,512,12 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p12s51275NL5"   ,512,12 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p12s51237p5NL5" ,512,12 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p12s51222p5NL5" ,512,12 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p12s51210NL5"   ,512,12 ,0.1333,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p12s5127p5NL5"  ,512,12 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,5)

# doAnalysis(process,process.newdigiseq,"p8s76875NL4"   ,768,8 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s76837p5NL4" ,768,8 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s76822p5NL4" ,768,8 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s76810NL4"   ,768,8 ,0.1333,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p8s7687p5NL4"  ,768,8 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,4)

# doAnalysis(process,process.newdigiseq,"p8s76875NL5"   ,768,8 ,2.0   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s76837p5NL5" ,768,8 ,0.5   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s76822p5NL5" ,768,8 ,0.3   ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s76810NL5"   ,768,8 ,0.1333,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p8s7687p5NL5"  ,768,8 ,0.1   ,cms.vint32(1,1,1,1,1,1),True,5)





#for the 15 sample
# doAnalysis(process,process.newdigiseq,"p6s51210NL3"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL3"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,1),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NL4"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL4"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,1),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NL5"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL5"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,1),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51210NL6"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NL6"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,1),True,6)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L3"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L3"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L4"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L4"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NM6L5"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM6L5"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,1,0),True,5)
# doAnalysis(process,process.newdigiseq,"p6s51210NM56L3"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM56L3"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s51210NM56L4"     ,512,6 ,2.0000,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM56L4"    ,512,6 ,0.75  ,cms.vint32(1,1,1,1,0,0),True,4)
# doAnalysis(process,process.newdigiseq,"p6s51210NM456L3"     ,512,6 ,2.0000,cms.vint32(1,1,1,0,0,0),True,3)
# doAnalysis(process,process.newdigiseq,"p6s5127p5NM456L3"    ,512,6 ,0.75  ,cms.vint32(1,1,1,0,0,0),True,3)


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

