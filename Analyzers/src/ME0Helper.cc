#include "../interface/ME0Helper.h"


namespace ME0Helper{
std::vector<const PSimHit*> getMatchedSimHits(const SimTrack * simTrack, const std::vector<PSimHit>&  simHitH){
	std::vector<const PSimHit*> hitList;
    for (const auto& simHit : simHitH) {
        if(simHit.eventId() != simTrack->eventId()) continue;
        if(simHit.trackId() != simTrack->trackId()) continue;
        hitList.push_back(&simHit);
    }
    return hitList;
}

SimMuons fillSimMuons(const std::vector<SimTrack>& simTracks, const std::vector<PSimHit>&  simHits) {
	SimMuons simMuons;
	for (const auto& simTrack : simTracks) {
		if(TMath::Abs(simTrack.type()) != 13) continue;
		simMuons.emplace_back(&simTrack, getMatchedSimHits(&simTrack,simHits));
	}
	return simMuons;
}


void fillDigiInfoMap(const  ME0DigiPreRecoCollection& digis, const ME0DigiPreRecoMap& digiMap,
		const ME0DigiPreRecoCollection& oldDigis, const std::vector<PSimHit>&  simHits, DigiInfoMap& digiInfoMap){
	digiInfoMap.clear();
	std::vector<int> matchedSHs(simHits.size(),0);
	for(const auto& ch_digis: digis ){
		const ME0DetId detId = ch_digis.first;
		auto ch_oldToNewMap = digiMap.get(detId);
		auto ch_oldDigis = oldDigis.get(detId);
		std::vector<DigiInfo>& ch_infoMap = digiInfoMap[detId];
		ch_infoMap.resize(int(ch_digis.second.second - ch_digis.second.first));
		for (ME0DigiPreRecoCollection::const_iterator idigi = ch_digis.second.first;
				idigi != ch_digis.second.second;idigi++) {
			ch_infoMap[int(idigi - ch_digis.second.first)].digi = &(*idigi);
		}
		for (ME0DigiPreRecoCollection::const_iterator idigi = ch_oldDigis.first;
				idigi != ch_oldDigis.second;idigi++) {
			const int oldIDX = int(idigi-ch_oldDigis.first);
			const int newIDX = ch_oldToNewMap.first[oldIDX];
			if(newIDX < 0) continue;
			ch_infoMap[newIDX].origDigis.push_back(&(*idigi));

			if(!idigi->prompt()) continue;
			double minDR2 = -1;
			int mSH = -1;
			for(unsigned int iSH = 0; iSH < simHits.size(); ++iSH){
				const auto& simHit = simHits[iSH];
				if(matchedSHs[iSH]) continue;
				if(simHit.particleType() != idigi->pdgid() ) continue;
				if(simHit.detUnitId() != detId ) continue;
				if(std::fabs(idigi->tof() - simHit.tof()) > 0.001 ) continue;
				if(std::fabs(simHit.entryPoint().x() - idigi->x()) > .001) continue;
				if(std::fabs(simHit.entryPoint().y() - idigi->y()) > .001) continue;
		          double dr2 = (simHit.entryPoint().x() - idigi->x())*(simHit.entryPoint().x() - idigi->x()) + (simHit.entryPoint().y() - idigi->y())*(simHit.entryPoint().y() - idigi->y());
		          if(minDR2 < 0 || dr2 < minDR2){
		            minDR2 = dr2;
		            mSH = iSH;
		          }
			}
			if(mSH >= 0){
				ch_infoMap[newIDX].simHits.push_back(&simHits[mSH]);
				matchedSHs[mSH] = 1;
			}

		}
	}
}

void associateSimMuons(SimMuons& simMuons, DigiInfoMap& digiInfoMap ) {
	for(unsigned int iM = 0; iM < simMuons.size(); ++iM){
		for(unsigned int iH = 0; iH < simMuons[iM].simHits.size(); ++iH){
			const auto* sh = simMuons[iM].simHits[iH];
			auto digiIt = digiInfoMap.find(sh->detUnitId());
			if(digiIt == digiInfoMap.end()) continue;
			bool matched = false;
			for(unsigned int iD = 0; iD < digiIt->second.size(); ++iD){
				auto& digi = digiIt->second[iD];
				for(unsigned int iH2 = 0; iH2 < digi.simHits.size(); ++iH2){
					if(digi.simHits[iH2] != sh) continue;
					digi.simTracks.push_back( simMuons[iM].track);
					digi.muonIDXs.push_back( iM);
					matched = true;
					break;
				}
				if(matched){
					simMuons[iM].digiInfos.emplace_back(sh->detUnitId(), iD);
					break;
				}
			}
		}
	}
}

int getRecHitIndex(const ME0RecHitCollection& recHits, const ME0RecHit& rh){
	int idx = -1;
	auto ch_RH = recHits.get(rh.me0Id());
	for(auto iR = ch_RH.first; iR != ch_RH.second; ++iR){
		if(iR->me0Id() != rh.me0Id()) continue;
		if(iR->tof() != rh.tof()) continue;
		if(!(iR->localPosition() == rh.localPosition())) continue;
		idx = int(iR -ch_RH.first);
		break;
	}
	return idx;
}

void fillSegmentCategories(const ME0SegmentCollection& segments,
		const ME0RecHitCollection& recHits,const DigiInfoMap& digiInfo,
		SimMuons& muons, SegmentCatMap& segCatMap) {

	segCatMap.clear();
	//First start by filling categorization info with general info
	//and association of muons to segments (many to one)
	for(auto iC = segments.id_begin(); iC != segments.id_end(); ++iC){
		auto ch_segs = segments.get(*iC);
		auto& ch_cats = segCatMap[*iC];
		for(auto iS = ch_segs.first; iS != ch_segs.second; ++iS){
			ch_cats.emplace_back(SEGMENT_ERROR,-1); //default is error and no muon
			const auto& rh = iS->specificRecHits();
			const unsigned int nH = rh.size();
			for(unsigned int iR = 0; iR < nH; ++iR){
				int rhIDX = getRecHitIndex(recHits,rh[iR]);
				const DigiInfo& di =  digiInfo.find(rh[iR].me0Id())->second[rhIDX];
				for(auto i : di.muonIDXs) {
					bool filled = false;
					for(unsigned int iMH = 0; iMH < muons[i].segments.size(); ++iMH){
						if(muons[i].segments[iMH].first.first != *iC ) continue;
						if(muons[i].segments[iMH].first.second != int(iS-ch_segs.first) ) continue;
						filled = true;
						muons[i].segments[iMH].second ++;
						break;
					}
					if(!filled){
						muons[i].segments.emplace_back(std::make_pair(*iC, int(iS-ch_segs.first)),1);
					}
				}
			}

		}
	}

	//Decide one main segment to associate each muon to
	for(unsigned int iM = 0; iM < muons.size(); ++iM){
		auto& msegs = muons[iM].segments;
		if(!msegs.size()) continue;
		std::sort(msegs.begin(),msegs.end(),
				[](const std::pair<std::pair<ME0DetId,int>,int >& a,const std::pair<std::pair<ME0DetId,int>,int >& b){
			return a.second > b.second;
		});
		segCatMap[msegs[0].first.first][msegs[0].first.second].second = iM;
	}

	//Now give a category to each segment
	for(auto iC = segments.id_begin(); iC != segments.id_end(); ++iC){
		auto ch_segs = segments.get(*iC);
		auto& ch_cats = segCatMap[*iC];
		for(auto iS = ch_segs.first; iS != ch_segs.second; ++iS){
			const auto& rh = iS->specificRecHits();
			const unsigned int nH = rh.size();
			auto& cat = ch_cats[iS - ch_segs.first];
			int nMM = 0; //From main muon
			int nOM = 0; //From additional muons
			int nT = 0;  //From prompt tracks
			int nN = 0;  //from neutrons
			int mT = 0;  //keep track of main prompt track (0 for none)
			for(unsigned int iR = 0; iR < nH; ++iR){
				int rhIDX = getRecHitIndex(recHits,rh[iR]);
				const DigiInfo& di =  digiInfo.find(rh[iR].me0Id())->second[rhIDX];
				int nPromptCont = 0;
				int promptPDGID = 0;
				for(const auto* od : di.origDigis){
					if(od->prompt()){
						nPromptCont++;
						promptPDGID = od->pdgid();
					}
				}
				if(di.muonIDXs.size() > 1 ) nOM++;
				else if (di.muonIDXs.size() && di.muonIDXs[0] != cat.second) nOM++;
				else if(di.muonIDXs.size()) nMM++;
				else if (nPromptCont){
					if(nPromptCont > 1) mT = 0;
					else if(nT && promptPDGID != mT ) mT = 0;
					else mT = promptPDGID;
					nT++;
				} else{
					nN++;
				}

			}
			if(nOM) cat.first = FAKE_MUON_MIX;
			else if(nMM){
				bool pure= nMM == int(nH);
				bool complete = nMM == int(muons[cat.second].digiInfos.size());
				if(complete && pure) cat.first = MUON_COMP_PURE;
				else if(complete && nT   ) cat.first = MUON_COMP_DIRTY_TRACK;
				else if(complete         ) cat.first = MUON_COMP_DIRTY_NEUT;
				else if(pure             ) cat.first = MUON_MISS_PURE;
				else if(nT               ) cat.first = MUON_MISS_DIRTY_TRACK;
				else                       cat.first = MUON_MISS_DIRTY_NEUT;
			} else if(nT == int(nH) && mT != 0 ) cat.first =  FAKE_TRACK_PURE;
			else if(nN == int(nH)) cat.first =  FAKE_NEUT_PURE;
			else cat.first =  FAKE_OTHER;
		}
	}
}



ME0DigiList getMatchedDigis(const std::vector<const PSimHit*>& simHits,
		const ME0DigiPreRecoCollection& oldDigis,const  ME0DigiPreRecoCollection * newDigis,const ME0DigiPreRecoMap*  digiMap, ME0DigiList * origDigiList, std::vector<int> * origDigiNewDigiIDX
) {
	if(origDigiList) origDigiList->clear();
	if(origDigiNewDigiIDX) origDigiNewDigiIDX->clear();
	ME0DigiList hitList;
	for(const auto* simHit : simHits){
        const ME0DetId detID = simHit->detUnitId();
        const ME0DigiPreRecoCollection::Range& range = oldDigis.get(detID);
        unsigned int iD = 0;
        double minDR2 = -1;
        const ME0DigiPreReco * digi = 0;
        for (ME0DigiPreRecoCollection::const_iterator idigi = range.first;
            idigi != range.second;idigi++) {
        	if(!idigi->prompt()) continue;
        	if(std::fabs(idigi->tof() - simHit->tof()) > 0.01 ) continue;
          double dr2 = (simHit->entryPoint().x() - idigi->x())*(simHit->entryPoint().x() - idigi->x()) + (simHit->entryPoint().y() - idigi->y())*(simHit->entryPoint().y() - idigi->y());
          if(dr2 >= .001) continue;
          if(minDR2 < 0 || dr2 < minDR2){
            digi = &(*idigi);
            minDR2 = dr2;
            iD = int(idigi-range.first);
          }
        }
        if(minDR2 < 0 ) {//std::cout << "Could not find matching digi!!!" <<std::endl;
        continue;}
        if(newDigis && digiMap){
            const ME0DigiPreRecoCollection::Range& newRange = newDigis->get(detID);
            const ME0DigiPreRecoMap::Range& oldToNewMap = digiMap->get(detID);
            int newIDX =  oldToNewMap.first[iD];
            if(origDigiList) origDigiList->emplace_back(detID,digi);
            if(origDigiList && origDigiNewDigiIDX){
            	if(newIDX < 0) std::cout << simHit->entryPoint().x() <<" "<< simHit->entryPoint().y()<< " "<<simHit->tof()<<" "<< ME0DetId(detID) <<" "<< digi->x()<<" "<<  digi->y() <<" "<< digi->tof() <<std::endl;

            	origDigiNewDigiIDX->push_back(newIDX);
            }
        	if(newIDX >= 0) hitList.emplace_back(detID,&newRange.first[newIDX]);
        }
        else  {
        	hitList.emplace_back(detID,digi);
        }
	}
    return hitList;
}


SimHitProperties getSimTrackProperties( const ME0Geometry* mgeom, const std::vector<const PSimHit* >& simHits) {

	SimHitProperties prop;

	std::map<ME0DetId,int> nChHit;
    std::vector<std::vector<const PSimHit*> > laysHit(6);
    for(const auto& h : simHits){
      laysHit[ME0DetId(h->detUnitId()).layer() - 1].push_back(h);
      nChHit[ME0DetId(h->detUnitId()).chamberId()]++;
    }

    std::vector<const PSimHit *> prunedHits;
    for(auto& l : laysHit){
      if(!l.size()) continue;
      prop.nLaysHit++;

      int iMax = -1;
      int max= 0;
      for(unsigned int iS = 0; iS < l.size(); ++iS){
    	  int nH = nChHit[ME0DetId(l[iS]->detUnitId()).chamberId() ];
    	  if(nH > max ){
    		  max = nH;
    		  iMax = iS;
    	  }
      }
      prunedHits.push_back(l[iMax] );
    }

    if(nChHit.size() > 1) prop.oneChamber = false;



//    std::vector<const PSimHit *> prunedHits;
//    for(auto& l : laysHit){
//      if(l.size() >0) {
//    	  prop.nLaysHit++;
//        prunedHits.push_back(l.front());
//      }
//    }
//
//    for(unsigned int iH = 0; iH < prunedHits.size(); ++iH){
//    	if(!iH) continue;
//    	if(ME0DetId(prunedHits[iH]->detUnitId()).chamber() !=  ME0DetId(prunedHits[0]->detUnitId()).chamber())
//    		prop.oneChamber = false;
//    }


    auto getGlobal = [&](const PSimHit * sh) -> GlobalPoint {
    	return mgeom->etaPartition(sh->detUnitId())->toGlobal(sh->entryPoint());
    };

    if(laysHit[0].size()){
    	prop.inFirst = true;
    	prop.lay1Eta = getGlobal(laysHit[0].front()).eta();
    	prop.lay1Phi = getGlobal(laysHit[0].front()).phi();
    }

    if(prop.nLaysHit > 1){
        auto getProj = [](double zValue, const GlobalPoint& up, const GlobalPoint& down) -> GlobalPoint {

          auto getCenter = [](double v1,double v2,double z1,double z2,double zc)->double{
            if(z2 == z1) return -99;
            double m = (v2 - v1)/(z2 - z1);
            double b =  (z2*v1 - z1*v2)/(z2 - z1);
            return m*zc+b;
          };

          const double zc = up.z() > 0 ? zValue : -1*zValue;
          double xc = getCenter(up.x(),down.x(),up.z(),down.z(),zc);
          double yc = getCenter(up.y(),down.y(),up.z(),down.z(),zc);
          return GlobalPoint(xc,yc,zc);
        };

        const int upInd = floor(prunedHits.size()/2);
        GlobalPoint centerUp = getGlobal(prunedHits[upInd]);
        GlobalPoint centerDown = getGlobal(prunedHits[upInd-1]);
        GlobalPoint centerPt = getProj(centerOfDet, centerUp, centerDown);
        GlobalPoint upPt = getProj(endOfDet, getGlobal(prunedHits[prunedHits.size() - 1]), getGlobal(prunedHits[prunedHits.size() - 2]));
        GlobalPoint downPt = getProj(beginOfDet, getGlobal(prunedHits[1]), getGlobal(prunedHits[0]));

        prop.cenEta = centerPt.eta();
        prop.cenPhi = centerPt.phi();
        prop.dPhi = TVector2::Phi_mpi_pi(upPt.phi() - downPt.phi());
        prop.dEta = upPt.eta() - downPt.eta();

    }
    return prop;
}

ME0Segment * buildSegment(const ME0Geometry* mgeom, const std::vector<std::pair<ME0DetId,const ME0DigiPreReco*> >& hits){
	if(hits.size() < 2) return 0;

	MuonSegFit::MuonRecHitContainer muonRecHits;
	std::vector<const ME0RecHit*> muonRH;

	const ME0DetId refID=hits[0].first;
	const ME0EtaPartition * refPart = mgeom->etaPartition(refID);
	for(const auto& hp : hits ){
	    const ME0EtaPartition * thePartition   =   mgeom->etaPartition(hp.first);

	    GlobalPoint gp = thePartition->toGlobal(LocalPoint(hp.second->x(),hp.second->y(),0));
	    const LocalPoint lp = refPart->toLocal(gp);

	    double errorXX = hp.second->ex()*hp.second->ex();
	    double errorYY = hp.second->ex()*hp.second->ex();
	    double errorXY = hp.second->corr()*hp.second->ex()*hp.second->ey();
	    LocalError le(std::max(errorXX,1.6E-5), std::max(errorXY,1.6E-6), std::max(errorYY,1.6E-5));

	    ME0RecHit* recHit = new ME0RecHit(refID,hp.second->tof(),lp,le);
	    muonRH.push_back(recHit);
	    MuonSegFit::MuonRecHitPtr trkRecHit(recHit);
	    muonRecHits.push_back(trkRecHit);
	}
	  MuonSegFit  * sfit_ = new MuonSegFit(muonRecHits);
	  sfit_->fit();

	  // obtain all information necessary to make the segment:
	  LocalPoint protoIntercept      = sfit_->intercept();
	  LocalVector protoDirection     = sfit_->localdir();
	  AlgebraicSymMatrix protoErrors = sfit_->covarianceMatrix();
	  double protoChi2               = sfit_->chi2();

	  ME0Segment * segment =new ME0Segment(muonRH, protoIntercept, protoDirection, protoErrors, protoChi2, 0.0, 5.0);

	  for (auto rh:muonRecHits) rh.reset();
	  delete sfit_;

	  return segment;
}

SegmentProperties getSegmentProperties(const GeomDet* loc, const ME0Segment * segment){

//	std::cout << segment->me0DetId() <<" "<< ME0DetId(segment->recHits()[0]->rawId()) <<" "<< refPart <<std::endl;
auto ld = segment->localDirection();
	auto lp = segment->localPosition();
	AlgebraicSymMatrix covMatrix = segment->parametersError();

	auto extrap = [&] (double extZ) -> SegmentExtrapPoint {
	    //covmatrix wors like 0: dx/dz dy/dz, x, y
	    double extX = lp.x()+extZ*ld.x()/ld.z();
	    double extY = lp.y()+extZ*ld.y()/ld.z();

	    double exx   =  extZ*extZ *covMatrix [0][0] + 2*extZ*covMatrix[0][2] + covMatrix[2][2];
	    double eyy   =  extZ*extZ *covMatrix [1][1] + 2*extZ*covMatrix[1][3] + covMatrix[3][3];
	    double exy   =  extZ*extZ *covMatrix [0][1] +   extZ*(covMatrix[1][2]+covMatrix[0][3]) + covMatrix[2][3];
	    double edxdx = covMatrix [0][0];
	    double edydy = covMatrix [1][1];
	    SegmentExtrapPoint extP;
	    extP.point =LocalPoint(extX,extY,extZ);
	    extP.pointError =LocalError(exx,exy,eyy);
	    extP.dxdz = ld.x()/ld.z();
	    extP.dydz = ld.y()/ld.z();
	    extP.cov_dxdz_dxdz = edxdx;
	    extP.cov_dydz_dydz = edydy;
	    return extP;
	  };

	  auto layCen =  loc->toGlobal(LocalPoint(0,0,0));
//	  auto locLow  = extrap( beginOfDet  - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());
//	  auto lochigh = extrap( endOfDet    - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());
//	  auto loccen  = extrap( centerOfDet - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z());

	  auto locLow  = extrap( beginOfDet -  (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z()    );
	  auto lochigh = extrap( endOfDet-   (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z()     );
	  auto loccen  = extrap(  centerOfDet - (layCen.z() < 0 ? -1.0 : 1.0) * layCen.z()   );

	  auto globlow  = loc->toGlobal(locLow .point);
	  auto globhigh = loc->toGlobal(lochigh.point);
	  auto globcen  = loc->toGlobal(loccen .point);

	  SegmentProperties properties;

	  properties.segAtCenter = loccen;
	  properties.initialPoint = lp;
	  properties.cenEta   =  globcen.eta();
	  properties.cenPhi   =  globcen.phi();
	  properties.beginEta =  globlow.eta();
	  properties.beginPhi =  globlow.phi();
	  properties.endEta   =  globhigh.eta();
	  properties.endPhi   =  globhigh.phi();
	  properties.dPhi     =  TVector2::Phi_mpi_pi(globhigh.phi() - globlow.phi());
	  properties.dEta     =  globhigh.eta() - globlow.eta();
	  return properties;
}


}


