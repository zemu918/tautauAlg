#include "tautauAlg/util/MyGamConv.h"
#include "/besfs5/groups/tauqcd/lihaolin09/7.0.3/workarea-7.0.3/Reconstruction/GammaConv/GammaConv-00-00-01/GammaConv/GammaConv.h"


bool checkGammConversion(HepPoint3D IP_Gam, RecMdcKalTrack* elp1Trk, RecMdcKalTrack* elm1Trk, NTuple::Array<double>& gamcon ){

	WTrackParameter wvelp1Trk, wvelm1Trk;
	RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
	wvelp1Trk = WTrackParameter(xmass[0], elp1Trk->getZHelixE(), elp1Trk->getZErrorE());
	RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
	wvelm1Trk = WTrackParameter(xmass[0], elm1Trk->getZHelixE(), elm1Trk->getZErrorE());

	HepPoint3D vx(0., 0., 0.);
	HepSymMatrix Evx(3, 0);
	double bx = 1E+6;
	double by = 1E+6;
	double bz = 1E+6;
	Evx[0][0] = bx*bx;
	Evx[1][1] = by*by;
	Evx[2][2] = bz*bz;
	VertexParameter vxparraw;
	vxparraw.setVx(vx);
	vxparraw.setEvx(Evx);

	//2020.1.12, only for the e+ e-
	VertexFit* vtxfit = VertexFit::instance();
	vtxfit->init();
	vtxfit->AddTrack(0,  wvelp1Trk);
	vtxfit->AddTrack(1,  wvelm1Trk);
	vtxfit->AddVertex(0, vxparraw,0, 1);

	if(vtxfit->Fit(0)) {

		HepPoint3D   vx_infit  = vtxfit->vx(0);
		HepSymMatrix Evx_infit = vtxfit->Evx(0);

		GammaConv gconv = GammaConv(elp1Trk->helix(), elm1Trk->helix(),IP_Gam);


		double  Rx_Gam   = -100.;
		double  Ry_Gam   = -100.;
		double  Rvxy_Gam = -100.;

		double  xconv    = -100. ;
		double  yconv    = -100. ;
		double  zconv    = -100. ;
		double  rconv    = -100. ;
		double  Dee_Mdc  = -100. ;
		double  Dee_Bp   = -100. ;

		double  xiep     = -100. ;
		double  deltaxy  = -100.;

		double  mconv    = -100.;
		double  econv    = -100.;
		double  thetaconv= -100.;
		double  phiconv  = -100.;
		double  angEE    = -100.;
		double  costhetaconv= -100.;

		Rx_Gam          = vx_infit.x();
		Ry_Gam          = vx_infit.y();
		Rvxy_Gam        = sqrt(Rx_Gam*Rx_Gam + Ry_Gam*Ry_Gam);


		xconv           = gconv.getRX();
		yconv           = gconv.getRY();
		zconv           = gconv.getRZ();
		rconv           = gconv.getRXY();
		xiep            = gconv.getXiep();
		deltaxy         = gconv.getDeltaXY();

		Dee_Mdc         = Get_Dee(elp1Trk->helix(),elm1Trk->helix(),5.9);
		Dee_Bp          = Get_Dee(elp1Trk->helix(),elm1Trk->helix(),3.15);
		mconv           = ( (elp1Trk->p4(xmass[0])) + (elm1Trk->p4(xmass[0]))).m();
		econv           = ( (elp1Trk->p4(xmass[0])) + (elm1Trk->p4(xmass[0]))).e();
		thetaconv       = ( (elp1Trk->p4(xmass[0])) + (elm1Trk->p4(xmass[0]))).theta();
		costhetaconv    = ( (elp1Trk->p4(xmass[0])) + (elm1Trk->p4(xmass[0]))).cosTheta();
		phiconv         = ( (elp1Trk->p4(xmass[0])) + (elm1Trk->p4(xmass[0]))).phi();

		gamcon[0]       = Rx_Gam;
		gamcon[1]       = Ry_Gam;
		gamcon[2]       = Rvxy_Gam;
		gamcon[3]       = xconv;
		gamcon[4]       = yconv;
		gamcon[5]       = rconv;

		double xconv_gamcon   = -100.0;
		double yconv_gamcon   = -100.0;
		double rconv_gamcon   = -100.0;
		double xconv1_gamcon  = -100.0;
		double yconv1_gamcon  = -100.0;
		double rconv1_gamcon  = -100.0;
		double xconv2_gamcon  = -100.0;
		double yconv2_gamcon  = -100.0;
		double rconv2_gamcon  = -100.0;

		if(xiep<10.) {
			xconv_gamcon  = xconv;
			yconv_gamcon  = yconv;
			rconv_gamcon  = rconv;
			xconv1_gamcon = xconv;
			yconv1_gamcon = yconv;
			rconv1_gamcon = rconv;
		} else {
			xconv_gamcon  = Rx_Gam;
			yconv_gamcon  = Ry_Gam;
			rconv_gamcon  = Rvxy_Gam;
			xconv2_gamcon = Rx_Gam;
			yconv2_gamcon = Ry_Gam;
			rconv2_gamcon = Rvxy_Gam;
		}
		gamcon[6]        = xconv_gamcon;
		gamcon[7]        = yconv_gamcon;
		gamcon[8]        = rconv_gamcon;
		gamcon[9]        = xconv1_gamcon;
		gamcon[10]       = yconv1_gamcon;
		gamcon[11]       = rconv1_gamcon;
		gamcon[12]       = xconv2_gamcon;
		gamcon[13]       = yconv2_gamcon;
		gamcon[14]       = rconv2_gamcon;
		gamcon[15]       = zconv;
		gamcon[16]       = Dee_Mdc;
		gamcon[17]       = Dee_Bp;
		gamcon[18]       = xiep;
		gamcon[19]       = deltaxy;

		gamcon[20]       = mconv;
		gamcon[21]       = econv;
		gamcon[22]       = thetaconv;
		gamcon[23]       = phiconv;
		gamcon[24]       = costhetaconv;
		gamcon[25]       = vtxfit->chisq(0);
	}
}
double Get_Dee(HepVector helixe,HepVector helixp,double r){
	HepPoint3D Inter_e = intersections(helixe,r);
	HepPoint3D Inter_p = intersections(helixp,r);
	double Dee = sqrt(pow((Inter_e.x()-Inter_p.x()),2)+pow((Inter_e.y()-Inter_p.y()),2));
	return Dee;
}

HepPoint3D intersections(HepVector helix,double r){
	double xc, yc,rc;
	xc = Get_HelixCenter(helix).x();
	yc = Get_HelixCenter(helix).y();
	rc = Get_HelixRadius(helix);
	double a,b,c;
	double d0 =helix[0];
	a = r; // radius of Mdc or Beam Pipe 
	b = d0 + rc;
	c = rc;
	double dphi = acos((a*a-b*b-c*c)/(-2.*b*c));
	double fltlen = dphi * rc;
	double phi  = acos((c*c-a*a-b*b)/(-2.*a*b));
	double x    = r*cos(phi+helix[1]);
	double y    = r*sin(phi+helix[1]);
	double z    = helix[3]+fltlen*helix[4];
	return HepPoint3D(x, y, z);
}
HepPoint3D Get_HelixCenter(HepVector helix){

	double bField = VertexFitBField::instance()->getBFieldZ(HepPoint3D(helix[0]*cos(helix[1]), helix[0]*sin(helix[1]), helix[3]));
	int charge = 1;
	const double alpha = -0.00299792458;
	double a = alpha * bField * charge;
	double pxy = charge/helix[2];
	double rad = pxy/a;
	double x = (helix[0] + rad) * cos(helix[1]);
	double y = (helix[0] + rad) * sin(helix[1]);
	double z = 0.0;
	return HepPoint3D(x, y, z);
}
double Get_HelixRadius(HepVector helix){
	double bField = VertexFitBField::instance()->getBFieldZ(HepPoint3D(helix[0]*cos(helix[1]), helix[0]*sin(helix[1]), helix[3]));
	int charge = 1;
	const double alpha = -0.00299792458;
	double a = alpha * bField * charge;
	double pxy = charge/helix[2];
	double radius = pxy/a;
	return radius;
}


