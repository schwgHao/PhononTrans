/*
 *\obtain the surface phonon green function of the electrode iterately.
 * 
 * using Eigen software package to solve the inverse and multiplication of matrice
 * 
 * Ref. Phys. Rev. B 78, 045434 (2008)
 *
 * 2017/06/10 Hao, Wang
 * */
#include <iostream>
#include <vector>
#include <complex>
#include <Dense>
#include <string>

using std::vector;
using std::complex;
using namespace Eigen;
using std::string;

#ifndef SURF_P_GF
#define SURF_P_GF

typedef vector<vector<vector<vector<vector<complex<double> > > > > > v5cd;
typedef vector<vector<vector<vector<complex<double> > > > > v4cd;

class surfphGF{
public:
	surfphGF(v5cd mfcbylz_, vector<double> omgv_, double delta_) 
		: mfcbylz(mfcbylz_), omgv(omgv_), delta(delta_){}
	void DR00(v4cd& d00r);
	void DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd& Tph);
	void writeDR00(string ssgf, double cellLz, int nlayers,
				   const vector<vector<double> >& kp2d, 
				   const v4cd& d00r);
private:
	v5cd mfcbylz;
	vector<double> omgv;
	double delta;
};

#endif
