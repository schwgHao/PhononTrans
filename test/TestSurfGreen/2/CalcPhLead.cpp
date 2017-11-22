#include <iostream>
#include <algorithm>
#include <complex>
#include "SurfacePhononGF.h"
#include <fstream>

using std::complex;
using std::fstream;
using std::vector;

int main(){
	
	std::cout << "Slice Supercell" << std::endl;
	const size_t nk = 1;
	const size_t nlayers = 2;
	const size_t naAtOnelayer = 1;
//	vector<vector<vector<vector<vector<complex<double> > > > > >mfcbylz(nk,
//		   vector<vector<vector<vector<complex<double> > > > > (nlayers,
//		          vector<vector<vector<complex<double> > > > (nlayers,
//				         vector<vector<complex<double> > > (3*naAtOnelayer, 
//								vector<complex<double> >(3*naAtOnelayer)))))
//	vector<vector<vector<vector<vector<complex<double> > > > > > mfcbylz
//	{{{{{0, 0, 0},{0, 0, 0},{0, 0, 1.2}},{{0, 0, 0},{0 , 0, 0},{0, 0, -0.6}}},
//	  {{{0, 0, 0},{0, 0, 0},{0, 0, -0.6}},{{0, 0, 0},{0, 0, 0},{0, 0, 1.2}}}}};
	
	vector<vector<vector<vector<vector<complex<double> > > > > > mfcbylz
	{{{{{1.2, 0, 0},{0, 1.2, 0},{0, 0, 1.2}},{{-0.6, 0, 0},{0 , -0.6, 0},{0, 0, -0.6}}},
	  {{{-0.6, 0, 0},{0, -0.6, 0},{0, 0, -0.6}},{{1.2, 0, 0},{0, 1.2, 0},{0, 0, 1.2}}}}};

	for(size_t i = 0; i < nlayers; ++i){
		for(size_t j = 0; j < nlayers; ++j){
			std::cout << "i, j " << i << " " << j << std::endl;  
			for(auto n1 : mfcbylz[0][i][j]){
				for(auto n2 : n1){
					std::cout << n2 << " ";
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}

	vector<double> omg = {1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8};

	const double delta = 1.0e-8;

	surfphGF sphGF(mfcbylz, omg, delta);
//	std::cout << "CalcPhLead: omg" << std::endl;
//	
//	for(auto i : omg) std::cout << i << std::endl;  
//	
	vector<vector<vector<vector<complex<double> > > > > d00r(nk,
		   vector<vector<vector<complex<double> > > > (omg.size(),
		          vector<vector<complex<double> > > (3*naAtOnelayer,
				         vector<complex<double> > (3*naAtOnelayer, 0.0))));

	sphGF.DR00(d00r);
	std::cout << "Generate surface GF" << std::endl;
	for(size_t j = 0; j < omg.size(); ++j){
		for(auto i : d00r[0][j]){
			for(auto k : i)std::cout << k << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}
