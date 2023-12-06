#ifndef MY_ISGOODTRACK_H
#define MY_ISGOODTRACK_H
#include "EvtRecEvent/EvtRecTrack.h"
#include "VertexFit/Helix.h"
using CLHEP::Hep3Vector;
bool isGoodTrack(EvtRecTrack* trk, Hep3Vector xorigin);
bool isSameId(EvtRecTrack* trk, SmartRefVector<EvtRecTrack> showers);
#endif
