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
//	string ss;
//	getline(fs, ss);
//	std::stringstream sts(ss);
	int natoms;
	fs.read((char*)&natoms, sizeof(int));
	
//	double xaa[natoms][3];
	double xaa[6][3];
//	fs.read((char*)&xaa[0][0], sizeof(double)*(natoms*3));
	fs.read(reinterpret_cast<char*>(xaa), sizeof(xaa));
//	vector<vector<double> >xa(natoms, vector<double>(3, 0.0));
//	copy(&xaa[0][0], &xaa[0][0] + natoms*3, xa.begin());
	vector<vector<double> >xa;
//	mcopy<double>(xaa, natoms, 3, xa);
	mcopy<double>((double*)xaa, 6, 3, xa);
	
	int xmassa[natoms];
	fs.read((char*)xmassa, sizeof(int)*natoms);
	vector<int> xmass(natoms, 1);
	copy(xmassa, xmassa + natoms, xmass.begin());	
	
	int nsc[3];
	fs.read((char*)nsc, sizeof(int)*3);
	double cella[3][3];
	fs.read(reinterpret_cast<char*>(cella), sizeof(cella));
//	fs.read((char*)&cell[0][0], sizeof(double)*9);
	double scella[3][3];
	fs.read(reinterpret_cast<char*>(scella), sizeof(scella));
//	fs.read((char*)&scell[0][0], sizeof(double)*9);
	vector<vector<double> > cell;
	vector<vector<double> > scell;
	mcopy<double>((double*)cella, 3, 3, cell);
	mcopy<double>((double*)scella, 3, 3, scell);
	
	int nkp;
	fs.read((char*)&nkp, sizeof(int));
	double kpa[nkp][4];
//	fs.read((char*)&kpa[0][0], sizeof(double)*(nkp*4));
	fs.read(reinterpret_cast<char*>(kpa), sizeof(kpa));
//	vector<vector<double> >kp(nkp, vector<double>(4, 0.0));
//	copy(&kpa[0][0], &kpa[0][0] + nkp*4, kp.begin());	
	vector<vector<double> >kp;
	mcopy<double>((double*)kpa, nkp, 4, kp);
//	sts >> nlayers >> natoms;
	fs.close();

	for(auto i : kp){
		for(auto j : i) cout << j << " ";
		cout << endl;
	}
	
	vector<vector<vector<complex<double> > > > mfc(nkp,
				 vector<vector<complex<double> > >(3*natoms, 
						  vector<complex<double> >(3*natoms, 0.0)));

	
	iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);
	
	if(Isbulk)
		ifc.ReadFC(mfc);
	for(auto k : mfc)
		for(auto i : k)
			for(auto j: i) std::cout << j << " ";

}
