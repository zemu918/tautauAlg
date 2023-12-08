#include "tautauAlg/util/MyUtil.h"

void setV4(NTuple::Array<double>& nt_A4, const HepLorentzVector& V4)
{
  nt_A4[0] = V4.x();
  nt_A4[1] = V4.y();
  nt_A4[2] = V4.z();
  nt_A4[3] = V4.e();
  nt_A4[4] = V4.perp();
  nt_A4[5] = V4.m();
  nt_A4[6] = V4.rho();
  nt_A4[7] = V4.cosTheta();
}

void setV4(NTuple::Matrix<double>& nt_M4, const HepLorentzVector& V4, int i)
{
  nt_M4[i][0] = V4.x();
  nt_M4[i][1] = V4.y();
  nt_M4[i][2] = V4.z();
  nt_M4[i][3] = V4.e();
  nt_M4[i][4] = V4.perp();
  nt_M4[i][5] = V4.m();
  nt_M4[i][6] = V4.rho();
  nt_M4[i][7] = V4.cosTheta();
}

HepLorentzVector getV4(RecEmcShower* shower)
{
  double eraw = shower->energy();
  double phi  = shower->phi();
  double the  = shower->theta();
  double vx   = eraw * sin(the) * cos(phi);
  double vy   = eraw * sin(the) * sin(phi);
  double vz   = eraw * cos(phi);

  HepLorentzVector V4(vx, vy, vz, eraw);

  return V4;
}


