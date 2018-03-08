#ifndef SURF_P_GF
#define SURF_P_GF
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
#include <queue>

using namespace Eigen;


typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > > v5cd;
typedef std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > v4cd;

class surfphGF{
public:
	surfphGF(v5cd mfcbylz_, std::vector<double> omgv_, double delta_) 
		: mfcbylz(mfcbylz_), omgv(omgv_), delta(delta_){}
	void DR00(v4cd& d00r);
	void DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd Tphnext, MatrixXcd& Tph);
	void DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd t1pre, MatrixXcd t2pre, MatrixXcd fpre, MatrixXcd ftpre, MatrixXcd Tphnext, MatrixXcd& Tph, const double alpha, const size_t exFreq, const size_t nStep, size_t step, std::queue<MatrixXcd>& qt1, std::queue<MatrixXcd>& qt2, std::queue<MatrixXcd>& qf1, std::queue<MatrixXcd>& qf2);
	void writeDR00(std::string ssgf, double cellLz, int nlayers,
				   const std::vector<std::vector<double> >& kp2d, 
				   const v4cd& d00r);
	void DeciTech(MatrixXcd& K0, MatrixXcd& K0s, 
		          MatrixXcd& K1, MatrixXcd& K2, MatrixXcd Iden);
private:
	v5cd mfcbylz;
	std::vector<double> omgv;
	double delta;
};
#endif
