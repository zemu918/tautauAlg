#include "tautauAlg/util/MyInfoShower.h"

bool MyInfoShower(RecEmcShower* aShower, NTuple::Matrix<double>& m_shower_info, int sz)
{
	double eraw = 0.0; eraw=  aShower->energy();
	double phi  = 0.0; phi =  aShower->phi();
	double the  = 0.0; the =  aShower->theta();
	HepLorentzVector p4_sig_shower( eraw * sin(the) * cos(phi),
			eraw * sin(the) * sin(phi),
			eraw * cos(the),
			eraw );

	m_shower_info[sz][0]=  p4_sig_shower.px();
	m_shower_info[sz][1]=  p4_sig_shower.py();
	m_shower_info[sz][2]=  p4_sig_shower.pz();
	m_shower_info[sz][3]=  p4_sig_shower.e();
	m_shower_info[sz][4]=  aShower->eSeed(); // the energy in the center of crystal
	m_shower_info[sz][5]=  aShower->e3x3();  // energy of 3x3 crystals, totally 9 ones
	m_shower_info[sz][6]=  aShower->e5x5();  // energy of 5x5 crystals, totally 25 ones
	m_shower_info[sz][7]=  aShower->theta();
	m_shower_info[sz][8]=  aShower->phi();

	m_shower_info[sz][9]  = aShower->secondMoment();
	m_shower_info[sz][10] = aShower->latMoment();
	m_shower_info[sz][11] = aShower->a20Moment();
	m_shower_info[sz][12] = aShower->a42Moment();
	m_shower_info[sz][13] = aShower->getEAll();
	m_shower_info[sz][14] = aShower->getELepton();
	m_shower_info[sz][15] = aShower->dE();
	m_shower_info[sz][16] = aShower->eSeed()/aShower->e3x3();
	m_shower_info[sz][17] = aShower->e3x3()/aShower->e5x5();
	m_shower_info[sz][18] = aShower->eSeed()/aShower->e5x5();

  return true;

}

