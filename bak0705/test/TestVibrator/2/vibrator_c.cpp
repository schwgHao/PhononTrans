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
	string ssinfo = slabel; 
	if(Isbulk) ssinfo += ".BulkInfo";
	else	   ssinfo += ".MolInfo";
	
	std::fstream fs(ssinfo, std::ios::in|std::ios::binary);
	if(!fs){
		std::cerr << "open "<< ssinfo <<" error!" << std::endl;
		abort();
	}
	
	int natoms;
	fs.read((char*)&natoms, sizeof(int));
	
	double xaa[natoms][3];
//	fs.read((char*)&xaa[0][0], sizeof(double)*(natoms*3));
	fs.read(reinterpret_cast<char*>(xaa), sizeof(xaa));
	vector<vector<double> >xa;
	mcopy<double>((double*)xaa, 6, 3, xa);
	
	int xmassa[natoms];
	fs.read((char*)xmassa, sizeof(int)*natoms);
	vector<int> xmass(natoms, 1);
	copy(xmassa, xmassa + natoms, xmass.begin());	
	
	int nsc[3];
	fs.read((char*)nsc, sizeof(int)*3);
	double cella[3][3];
//	fs.read((char*)&cell[0][0], sizeof(double)*9);
	fs.read(reinterpret_cast<char*>(cella), sizeof(cella));
	vector<vector<double> > cell;
	mcopy<double>((double*)cella, 3, 3, cell);
	double scella[3][3];
//	fs.read((char*)&scell[0][0], sizeof(double)*9);
	fs.read(reinterpret_cast<char*>(scella), sizeof(scella));
	vector<vector<double> > scell;
	mcopy<double>((double*)scella, 3, 3, scell);
	
	int nkp;
	fs.read((char*)&nkp, sizeof(int));
	double kpa[nkp][4];
	fs.read(reinterpret_cast<char*>(kpa), sizeof(kpa));
	vector<vector<double> >kp;
	mcopy<double>((double*)kpa, nkp, 4, kp);
	fs.close();

/* Read force constant matrix */	
	vector<vector<vector<complex<double> > > > mfc(nkp,
				 vector<vector<complex<double> > >(3*natoms, 
						  vector<complex<double> >(3*natoms, 0.0)));

	
	iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);

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
