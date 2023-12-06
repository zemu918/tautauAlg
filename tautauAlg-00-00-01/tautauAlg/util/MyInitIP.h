#ifndef MyInitIP_H
#define MyInitIP_H

#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "CLHEP/Vector/ThreeVector.h"
using CLHEP::Hep3Vector;
#include "CLHEP/Geometry/Point3D.h"
#ifdef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom:Point3D<double> HepPoint3D;
#endif
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/KalmanKinematicFit.h"

class MyInitIP
{
  public:
    MyInitIP();
    ~MyInitIP();
    void InitIP(Hep3Vector& xorigin, VertexParameter&  m_privxpar);
  protected:
    IVertexDbSvc*  vtxsvc;

    //IDataProviderSvc* eventSvc;

  private:
};
#endif


