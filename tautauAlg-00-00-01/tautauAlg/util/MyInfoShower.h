#ifndef MY_INFOSHOWER_H
#define MY_INFOSHOWER_H
#include "tautauAlg/util/MyConst.h"
using CLHEP::Hep3Vector;

bool MyInfoShower(RecEmcShower* aShower, NTuple::Matrix<double>& m_shower_info, int sz);
#endif

