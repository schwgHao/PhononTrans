#include <iostream>
#include "ioFC.h"
#include <algorithm>
#include <complex>
#include "SurfacePhononGF.h"
#include <fstream>

using std::complex;
using std::fstream;

void CalcPhLead(string slabel, const vector<double>& omg, double delta){
//	if(Isbulk) ssinfo += ".BulkInfo";
//	else	   ssinfo += ".MolInfo";
	
	
	iofc ifc(slabel, true);

	int natoms = ifc.Nat();
/* Read force constant matrix */	
	vector<vector<vector<complex<double> > > > mfc(1,
				 vector<vector<complex<double> > >(3*natoms, 
						  vector<complex<double> >(3*natoms, 0.0)));

//	iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);

	ifc.ReadFC(mfc);
	
	
	vector<vector<double> >kp2d = ifc.kp2d();
	
	vector<vector<int> > atominLayer = ifc.MatCutbyLayer();
	int nk = kp2d.size();
	int nlayers = atominLayer.size();
	double cellLz = ifc.outermost();
	int naAtOnelayer = natoms/nlayers;
	vector<vector<vector<vector<vector<complex<double> > > > > >mfcbylz(nk,
		   vector<vector<vector<vector<complex<double> > > > > (nlayers,
		          vector<vector<vector<complex<double> > > > (nlayers,
				         vector<vector<complex<double> > > (3*naAtOnelayer, 
								vector<complex<double> >(3*naAtOnelayer, 0.0)))));

	ifc.mfcSlice(mfcbylz, mfc, atominLayer);

	surfphGF sphGF(mfcbylz, omg, delta);

	vector<vector<vector<vector<complex<double> > > > > d00r(nk,
		   vector<vector<vector<complex<double> > > > (omg.size(),
		          vector<vector<complex<double> > > (3*naAtOnelayer,
				         vector<complex<double> > (3*naAtOnelayer, 0.0))));

	sphGF.DR00(d00r);

	/* write surface phonon green function to  */
	string ssgflft = "bulksgf.lft";	
	
	sphGF.writeDR00(ssgflft, cellLz, nlayers, kp2d, d00r);
	
	string ssgfrt = "bulksgf.rt";
	/*if the two leads are same, just copy left lead to right*/
	if(true){
		fstream ifs(ssgflft, std::ios::in|std::ios::binary);
		fstream ofs(ssgfrt, std::ios::out|std::ios::binary);
		ofs << ifs.rdbuf();
		ifs.close();
		ofs.close();
	}
}
