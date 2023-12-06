#ifndef MY_GAMCONV_H
#define MY_GAMCONV_H
#include "tautauAlg/util/MyConst.h"

bool checkGammConversion(HepPoint3D IP_Gam, RecMdcKalTrack* elp1Trk, RecMdcKalTrack* elm1Trk, NTuple::Array<double>& gamcon );
double     Get_Dee(HepVector helixe,HepVector helixp,double r);
HepPoint3D intersections(HepVector helix,double r);
HepPoint3D Get_HelixCenter(HepVector helix);
double     Get_HelixRadius(HepVector helix);

#endif

