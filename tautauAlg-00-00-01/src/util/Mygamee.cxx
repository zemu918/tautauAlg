#include "tautauAlg/util/Mygamee.h"
#include "tautauAlg/util/MyIsGoodshower.h"
bool checkPi0Dalitz(CDPhotonList& photonList, HepLorentzVector p4lepp, HepLorentzVector p4lepm, double  mpi0masspdg, NTuple::Array<double>& m_info_gamee  ) {
	RecEmcShower* dalitzPhoton_pi0 = NULL;
	double chisq_deltampi0         = 1000.;
	int tmp_pi0_gam_id             = -100;
	HepLorentzVector p4Pi0Dalitz(0., 0., 0., 0.);
	HepLorentzVector p4PhotonDalitz_pi0(0., 0., 0., 0.);
	for ( CDPhotonList::iterator it = photonList.particle_begin(); it != photonList.particle_end(); it++ ) {

		EvtRecTrack* aRecTrack = const_cast<EvtRecTrack*>((*it).particle().photon() );
		if(!isGoodShower(aRecTrack)) continue;
		int pi0_gam_id = -100;
		pi0_gam_id = aRecTrack->trackId();
		RecEmcShower* aPhotonShower = aRecTrack->emcShower();
		double eraw =  aPhotonShower->energy();
		double phi  =  aPhotonShower->phi();
		double the  =  aPhotonShower->theta();
		HepLorentzVector p4photon( eraw * sin(the) * cos(phi),
				eraw * sin(the) * sin(phi),
				eraw * cos(the),
				eraw );

		HepLorentzVector tmp_p4Pi0(0., 0., 0., 0.);
    tmp_p4Pi0 = p4lepp + p4lepm + p4photon;
		double tmp_deltampi0  =  -100.;
    tmp_deltampi0         = fabs(tmp_p4Pi0.m() - mpi0masspdg );
		if (tmp_deltampi0 < chisq_deltampi0 ) {
			dalitzPhoton_pi0    = aPhotonShower;
			chisq_deltampi0     = tmp_deltampi0;
			p4Pi0Dalitz         = tmp_p4Pi0;
			p4PhotonDalitz_pi0  = p4photon;
			tmp_pi0_gam_id      = pi0_gam_id;
		}
	}

	if ( dalitzPhoton_pi0 != NULL && chisq_deltampi0 < 10. ) {
		double gammaE1    = dalitzPhoton_pi0->energy();
		double gammaE2    = (p4lepp + p4lepm).e();
		m_info_gamee[0]   = (gammaE1 - gammaE2)/p4Pi0Dalitz.rho();
		m_info_gamee[1]   = p4Pi0Dalitz.m();
		m_info_gamee[2]   = p4Pi0Dalitz.rho();
		m_info_gamee[3]   = gammaE1;
		m_info_gamee[4]   = cos(dalitzPhoton_pi0->theta() );
		m_info_gamee[5]   = p4PhotonDalitz_pi0.px();
		m_info_gamee[6]   = p4PhotonDalitz_pi0.py();
		m_info_gamee[7]   = p4PhotonDalitz_pi0.pz();
		m_info_gamee[8]   = p4PhotonDalitz_pi0.e();
		m_info_gamee[9]   = p4Pi0Dalitz.px();
		m_info_gamee[10]  = p4Pi0Dalitz.py();
		m_info_gamee[11]  = p4Pi0Dalitz.pz();
		m_info_gamee[12]  = p4Pi0Dalitz.e();
		m_info_gamee[13]  = tmp_pi0_gam_id;
		m_info_gamee[14]  = p4Pi0Dalitz.cosTheta();
		m_info_gamee[15]  = p4Pi0Dalitz.perp();
	}
}
