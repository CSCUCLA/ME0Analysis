import FWCore.ParameterSet.Config as cms
import re

from Configuration.StandardSequences.Eras import eras
process = cms.Process('VALIDATION',eras.Phase2C2)

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
def ME0ReReco(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1),  doRUSegmentAlgo = True, minSegmentLayers = 4,timeWindow = 25, trackName = "generalTracks") :
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
    
    segM_name          = name + "SegmentMatching"
    setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag(seg_name),maxDiffX = cms.double(15.0),maxDiffY = cms.double(15.0),tracksTag = cms.InputTag(trackName)  ))
    seq += getattr(process, segM_name)
      
    me0Muon_name          = name + "Me0Muon"
    setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
    seq += getattr(process, me0Muon_name)
    
def ME0ReRecoOnlyME0Muons(process, seq, name, trackName = "generalTracks") :    
    segM_name          = name + "SegmentMatching"
    setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag("me0Segments"),maxDiffX = cms.double(15.0),maxDiffY = cms.double(15.0),tracksTag = cms.InputTag(trackName)  ))
    seq += getattr(process, segM_name)
      
    me0Muon_name          = name + "Me0Muon"
    setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
    seq += getattr(process, me0Muon_name)
def doAnalysisOnlyME0Muons(process, seq, name, trackName = "generalTracks") :
    ME0ReRecoOnlyME0Muons(process,seq,name,trackName)
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0TrackMatchingTreeMakerForPOG",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string("simMuonME0ReDigis"),
        segmentCollection = cms.string("me0Segments"),
        recHitCollection = cms.string("me0RecHits"),
        muonsTag = cms.string(name+"SegmentMatching"),
        runName           = cms.untracked.string(name+"_"),
        trackCollection   = cms.string(trackName)
        )
    )
    seq += getattr(process, anName)
    
def doAnalysis(process, seq, name, nStrips = 768, nPartitions = 8, neutBKGAcc = 2.0, layerRO =cms.vint32(1,1,1,1,1,1), doRUSegmentAlgo = True, minSegmentLayers = 4,timeWindow = 25, trackName = "generalTracks") :
    ME0ReReco(process,seq,name,nStrips,nPartitions,neutBKGAcc,layerRO,doRUSegmentAlgo,minSegmentLayers,timeWindow,trackName)
    anName = name + "Analysis"
    setattr( process, anName, cms.EDAnalyzer("ME0TrackMatchingTreeMakerForPOG",
        outFileName       = cms.untracked.string(re.sub(r'(.*)\.root',r'\1_'+name+'.root',options.outputFile)),   
        newDigiCollection = cms.string(name+"Digis"),
        segmentCollection = cms.string(name+"Segments"),
        recHitCollection = cms.string(name+"RecHits"),
        muonsTag = cms.string(name+"SegmentMatching"),
        runName           = cms.untracked.string(name+"_"),
        trackCollection   = cms.string(trackName)
        )
    )
    seq += getattr(process, anName)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.selectedVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
    filter = cms.bool(False)
)

process.muonTiming = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("abs(time.timeAtIpInOut) < (12.5 + abs(time.timeAtIpInOutErr))"),
    filter = cms.bool(False)
)

process.bestMuonLooseModExt = cms.EDProducer("MuonTrackProducer",
   muonsTag = cms.InputTag("muonTiming"),
   inputDTRecSegment4DCollection = cms.InputTag("dt4DSegments"),
   inputCSCSegmentCollection = cms.InputTag("cscSegments"),
   vtxTag = cms.InputTag("selectedVertices"),
   selectionTags = cms.vstring('All'),
   trackType = cms.string('globalTrackLooseModExt')
)


import PhysicsTools.RecoAlgos.recoTrackSelector_cfi
process.probeTracks = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.probeTracks.quality = cms.vstring('highPurity')
process.probeTracks.tip = cms.double(3.5)
process.probeTracks.lip = cms.double(30.)
process.probeTracks.ptMin = cms.double(2.0)
process.probeTracks.minRapidity = cms.double(-3.0)
process.probeTracks.maxRapidity = cms.double(3.0)

import SimMuon.MCTruth.MuonAssociatorByHits_cfi
process.trackingParticleRecoTrackAsssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
process.trackingParticleRecoTrackAsssociation.tpTag = 'mix:MergedTrackTruth'
process.trackingParticleRecoTrackAsssociation.tracksTag = 'generalTracks'
process.trackingParticleRecoTrackAsssociation.UseTracker = True
process.trackingParticleRecoTrackAsssociation.UseMuon = False
# import SimMuon.MCTruth.MuonAssociatorByHits_cfi
# process.trackingParticleRecoTrackAsssociation = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone()
# process.trackingParticleRecoTrackAsssociation.tpTag = 'mix:MergedTrackTruth'
# process.trackingParticleRecoTrackAsssociation.tracksTag = 'bestMuonLooseModExt'
# process.trackingParticleRecoTrackAsssociation.UseTracker = True
# process.trackingParticleRecoTrackAsssociation.UseMuon = True
# process.trackingParticleRecoTrackAsssociation.acceptOneStubMatchings = False

process.genMuons = cms.EDFilter("PdgIdCandViewSelector",
    src = cms.InputTag("genParticles"),
    pdgId = cms.vint32(13, -13)
)
 
## filter for sim level
process.genMuonsGEM = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genMuons"),
    cut = cms.string("abs(eta)<3.0 & abs(eta)>1.8"),   
    filter = cms.bool(True)
)

process.newdigiseq  = cms.Sequence()
#run with general tracks
# doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),True)
# doAnalysisOnlyME0Muons(process,process.newdigiseq,"p8s384","generalTracks")
# process.p = cms.Path(process.genMuons* process.genMuonsGEM * process.selectedVertices*process.muonTiming * process.bestMuonLooseModExt* process.trackingParticleRecoTrackAsssociation* process.newdigiseq)


#run with HP tracks:
# doAnalysis(process,process.newdigiseq,"p8s384" ,384,8 ,2.0,cms.vint32(1,1,1,1,1,1),True,4,25,"probeTracks")
doAnalysisOnlyME0Muons(process,process.newdigiseq,"p8s384","probeTracks")

process.trackingParticleRecoTrackAsssociation.tracksTag = cms.InputTag('probeTracks')
process.p = cms.Path(process.probeTracks * process.trackingParticleRecoTrackAsssociation* process.newdigiseq)
