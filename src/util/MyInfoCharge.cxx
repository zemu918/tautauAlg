#include "tautauAlg/util/MyInfoCharge.h"
#include "HadronInfo/TrackInfo.h"
void MyInfoCharge(TrackInfo track1Info, NTuple::Array<double>& charge_info)
{
	HepLorentzVector pmu41 (0., 0., 0., 0.);
	pmu41               = track1Info.p4();
	charge_info[0]      = track1Info.dedxchie();
	charge_info[1]      = track1Info.dedxchimu();
	charge_info[2]      = track1Info.dedxchipi();
	charge_info[3]      = track1Info.dedxchik();
	charge_info[4]      = track1Info.dedxchip();
	charge_info[5]      = pmu41.rho();
	charge_info[6]      = pmu41.cosTheta();
	charge_info[7]      = track1Info.emctrke();
	charge_info[8]      = track1Info.emctrkz();
	charge_info[9]      = track1Info.emctrky();
	charge_info[10]     = track1Info.emctrkx();
	charge_info[11]     = track1Info.emctrkphi();
	charge_info[12]     = track1Info.emctrktime();
	charge_info[13]     = track1Info.emctrktheta();
	charge_info[14]     = track1Info.eoptrk();
	charge_info[15]     = track1Info.muLay();
	charge_info[16]     = track1Info.muDep();
	charge_info[17]     = track1Info.muNhits();
	charge_info[18]     = track1Info.muLastLay();
	charge_info[19]     = track1Info.muMaxHitsInLay();
	charge_info[20]     = pmu41.px();
	charge_info[21]     = pmu41.py();
	charge_info[22]     = pmu41.pz();
	charge_info[23]     = pmu41.e();
	charge_info[24]     = pmu41.perp();
	charge_info[25]     = track1Info.charge();
}

