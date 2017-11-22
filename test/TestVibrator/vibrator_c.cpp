#include <iostream>
#include <fstream>
#include "ioFC.h"
#include "FuncUtils.h"
#include <algorithm>
#include <complex>
#include "SurfacePhononGF.h"

using std::complex;
using std::cout;
using std::endl;
using std::fstream;

extern "C" void vibrator_c(char* slabel, bool Isbulk);

void vibrator_c(char* slabel_, bool Isbulk){
	vector<string> s1 = strsplit((string)slabel_);
	string slabel = s1[0]; // the difference of string expression between fortran and C. 
//	string ssinfo = slabel; 
//	if(Isbulk) ssinfo += ".BulkInfo";
//	else	   ssinfo += ".MolInfo";
	
	
	iofc ifc(slabel);

	int natoms = ifc.Nat();
/* Read force constant matrix */	
	vector<vector<vector<complex<double> > > > mfc(1,
				 vector<vector<complex<double> > >(3*natoms, 
						  vector<complex<double> >(3*natoms, 0.0)));

//	iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);

	if(Isbulk){
		ifc.ReadFC(mfc);
	}
	
	
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

	vector<double> omg(20);
	for(int i = 0; i < omg.size(); i++)
		omg[i] = 1.0e-4 + 1.0e-1*i/omg.size();
	double delta = 1.0e-6;
	surfphGF sphGF(mfcbylz, omg, delta);

	vector<vector<vector<vector<complex<double> > > > > d00r(nk,
		   vector<vector<vector<complex<double> > > > (omg.size(),
		          vector<vector<complex<double> > > (3*naAtOnelayer,
				         vector<complex<double> > (3*naAtOnelayer, 0.0))));

	sphGF.DR00(d00r);

	/* write surface phonon green function to  */
	string ssgf = slabel + ".SGF";	
	
	sphGF.writeDR00(ssgf, cellLz, nlayers, kp2d, d00r);

	fstream ifs(ssgf, std::ios::in|std::ios::binary);
	double cellLzr;
	ifs.read((char*)&cellLzr, sizeof(double));
	cout << cellLzr << endl;
	int nlayersr;
	ifs.read((char*)&nlayersr, sizeof(int));
	cout << nlayersr << endl;
	auto kp2dr = kp2d;
	for(int i = 0; i < kp2dr.size(); i++)
		ifs.read((char*)&kp2dr[i], kp2dr.size()*sizeof(double));

	auto d00rr = d00r;
	for(int i = 0; i < d00r.size(); i++)
		for(int j = 0; j < d00r[0].size(); j++)
			for(int k = 0; k < d00r[0][0].size(); k++){
		//		for(int l = 0; l < d00r[0][0][0].size(); l++)
		//			d00rr[i][j][k][l] = complex<double>(0.0, 0.0);
				ifs.read(reinterpret_cast<char*>(&d00rr[i][j][k][0]), d00rr[0][0][0].size()*sizeof(complex<double>));
			}
	ifs.close();
	
	
	for(int k = 0; k < nk; k++){
		cout << "k = " << k << endl;
		for(int o = 0; o < omg.size(); o++){
			cout << "omg = " << omg[o] * 519.6 << " cm-1" << endl;
			for(int i = 0; i < 3*naAtOnelayer; i++){
				for(int j = 0; j < 3*naAtOnelayer; j++)
					cout << d00r[k][o][i][j] << " ";
				cout << endl;
				for(int j = 0; j < 3*naAtOnelayer; j++)
					cout << d00rr[k][o][i][j] << " ";
				cout << endl;
				cout << endl;
			}
		}
	}
//	for(int k = 0; k < nk; k++){
//		cout << "kp: " << k << endl;
//		for(int l1= 0; l1< nlayers; l1++){
//			for(int l2= 0; l2< nlayers; l2++){
//				cout << "layer: " << l1<<", "<< l2 << endl;
//				for(int i = 0; i < 3*naAtOnelayer; i++){
//					for(int j = 0; j < 3*naAtOnelayer; j++)
//						cout << mfcbylz[k][l1][l2][i][j] << " ";
//					cout << endl;
//				}
//			}
//		}
//	}
}
