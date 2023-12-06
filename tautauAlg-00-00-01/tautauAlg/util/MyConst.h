#ifndef MY_CONST_H
#define MY_CONST_H

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/IJobOptionsSvc.h"
//#include "MeasuredEcmsSvc/IMeasuredEcmsSvc.h"

#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EventModel/Event.h"


#include "EvTimeEvent/RecEsTime.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "DstEvent/TofHitStatus.h"

#include "MdcRecEvent/RecMdcKalTrack.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "MucRecEvent/RecMucTrack.h"

#include "McTruth/McParticle.h"


#include "TMath.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;

#include "CLHEP/Geometry/Point3D.h"

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom:Point3D<double> HepPoint3D;
#endif

//Vertex
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/BField.h"
//PID
#include "ParticleID/ParticleID.h"


#include "HadronInfo/TrackInfo.h"
#include "HadronInfo/ShowerInfo.h"


#include "McDecayModeSvc/McDecayModeSvc.h"
#include "DTagAlg/LocalPhotonSelector.h"

#include<vector>

#include "BesDChain/CDDecayList.h"
#include "BesDChain/CDPhotonList.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iostream>
#include <fstream>

//mass of e+, mu+, pi+, K+, proton
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double mpi0pdg    = 0.1349766;
const double mkspdg     = 0.497611;
const double metapdg    = 0.547862;
const double momgpdg    = 0.78265;
const double md0pdg     = 1.86483;
const double mdppdg     = 1.86958;
const double mdspdg     = 1.96827;
const double mdsstpdg   = 2.1121;
const double mjpsipdg   = 3.096900;


#endif



