import FWCore.ParameterSet.Config as cms

#Name  
def ME0ReReco(process, seq, name, nStrips, nPartitions) :
    
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    
    process.RandomNumberGeneratorService.simMuonME0NewGeoDigis = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    )
    digi_name          = name + "_digis"
    setattr( process, digi_name, simMuonME0NewGeoDigis.clone(numberOfSrips=cms.uint32(nStrips), numberOfPartitions =cms.uint32(nPartitions)))
    getattr(process.RandomNumberGeneratorService, digi_name, cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(7654326)
    ))
    seq += getattr(process, digi_name)
    
    
    rh_name          = name + "_recHits"
    setattr( process, rh_name, process.me0RecHits.clone(me0DigiLabel = cms.InputTag(digi_name)))
    seq += getattr(process, rh_name)
    
    seg_name          = name + "_segments"
    setattr( process, seg_name, process.me0Segments.clone(algo_type = 1,me0RecHitLabel = cms.InputTag(rh_name)))
    seq += getattr(process, seg_name)
    
    segM_name          = name + "_segmentMatching"
    setattr( process, segM_name, process.me0SegmentMatching.clone(me0SegmentTag = cms.InputTag(seg_name)))
    seq += getattr(process, segM_name)
    
    me0Muon_name          = name + "_me0Muon"
    setattr( process, me0Muon_name, process.me0MuonConverting.clone(me0SegmentTag = cms.InputTag(segM_name)))
    seq += getattr(process, me0Muon_name)