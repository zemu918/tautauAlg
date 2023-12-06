#include "tautauAlg/util/MyInitIP.h"

using namespace std;

MyInitIP::MyInitIP()
{}
MyInitIP::~MyInitIP()
{}

void MyInitIP::InitIP(Hep3Vector& xorigin, VertexParameter&  m_privxpar)
{
	xorigin = Hep3Vector(0,0,0);
	HepPoint3D pvx(0., 0., 0.);
	HepSymMatrix pEvx(3, 0);

	Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
	if (vtxsvc->isVertexValid()) {
		double* vertex = vtxsvc->PrimaryVertex(); //vertex[0]=vx;vertex[1]=vy;vertex[2]=vz;
		double* vertexsigma = vtxsvc->SigmaPrimaryVertex();
		xorigin.setX(vertex[0]);
		xorigin.setY(vertex[1]);
		xorigin.setZ(vertex[2]);
		pvx[0]     = vertex[0];
		pvx[1]     = vertex[1];
		pvx[2]     = vertex[2];
		pEvx[0][0] = vertexsigma[0]*vertexsigma[0];
		pEvx[1][1] = vertexsigma[1]*vertexsigma[1];
		pEvx[2][2] = vertexsigma[2]*vertexsigma[2];
		m_privxpar.setVx(pvx);
		m_privxpar.setEvx(pEvx);
	}
	return;
}

