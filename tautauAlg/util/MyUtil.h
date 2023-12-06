#ifndef MY_UTIL_H
#define MY_UTIL_H
#include <vector>
#include "GaudiKernel/NTuple.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "EmcRecEventModel/RecEmcShower.h"

void setV4(NTuple::Array<double>& nt_A4, const HepLorentzVector& V4);
void setV4(NTuple::Matrix<double>& nt_M4, const HepLorentzVector& V4, int i);
HepLorentzVector getV4(RecEmcShower* shower);

#endif

