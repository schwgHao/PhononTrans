#include <iostream>
#include <fstream>
#include "ioFC.h"
#include "FuncUtils.h"
#include <algorithm>
#include <complex>

using std::complex;
using std::cout;
using std::endl;
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
	int naAtOnelayer = natoms/nlayers;
	vector<vector<vector<vector<vector<complex<double> > > > > >mfcbylz(nk,
		   vector<vector<vector<vector<complex<double> > > > > (nlayers,
		          vector<vector<vector<complex<double> > > > (nlayers,
				         vector<vector<complex<double> > > (3*naAtOnelayer, 
								vector<complex<double> >(3*naAtOnelayer, 0.0)))));

	ifc.mfcSlice(mfcbylz, mfc, atominLayer);

	for(int k = 0; k < nk; k++){
		cout << "kp: " << k << endl;
		for(int l1= 0; l1< nlayers; l1++){
			for(int l2= 0; l2< nlayers; l2++){
				cout << "layer: " << l1<<", "<< l2 << endl;
				for(int i = 0; i < 3*naAtOnelayer; i++){
					for(int j = 0; j < 3*naAtOnelayer; j++)
						cout << mfcbylz[k][l1][l2][i][j] << " ";
					cout << endl;
				}
			}
		}
	}
}
