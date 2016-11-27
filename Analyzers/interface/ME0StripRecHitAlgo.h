#ifndef ME0StripRecHitAlgo_h
#define ME0StripRecHitAlgo_h


#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHit.h"
#include "DataFormats/Common/interface/OwnVector.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "FWCore/Framework/interface/ESHandle.h"

class TrapezoidalStripTopology;

class ME0StripRecHitAlgo  {
 public:
  /// Constructor
  ME0StripRecHitAlgo(const edm::ParameterSet& config);

  /// Destructor
  virtual ~ME0StripRecHitAlgo();

  // Operations

  /// Pass the Event Setup to the algo at each event
  virtual void setES(const edm::EventSetup& setup);


  virtual bool compute(const ME0DigiPreReco& digi,
                       LocalPoint& point,
                       LocalError& error) const;

  /// Build all hits in the range associated to the me0Id, at the 1st step.
  virtual edm::OwnVector<ME0RecHit> reconstruct(const ME0DetId& me0Id,
                                                   const ME0DigiPreRecoCollection::Range& digiRange);

  TrapezoidalStripTopology * buildTopo(std::vector<float>& _p) const;

  edm::ESHandle<ME0Geometry> me0Geom;
  const unsigned int nP;
  const unsigned int nStrips;

};

#endif
