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
#include <map>
#include <utility>
#include <queue>




typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > > v5cd;
typedef std::vector<std::vector<std::vector<std::vector<std::complex<double> > > > > v4cd;

class surfphGF{
public:
	surfphGF(v5cd mfcbylz_, std::vector<double> omgv_, double delta_) 
		: mfcbylz(mfcbylz_), omgv(omgv_), delta(delta_){}
	void DR00(v4cd& d00r);
	void DR00(Eigen::MatrixXcd t0, Eigen::MatrixXcd t0t, Eigen::MatrixXcd Tphnext, Eigen::MatrixXcd& Tph);
	void writeDR00(std::string ssgf, double cellLz, int nlayers,
				   const std::vector<std::vector<double> >& kp2d, 
				   const v4cd& d00r);
	void DeciTech(Eigen::MatrixXcd& K0, Eigen::MatrixXcd& K0s, 
		          Eigen::MatrixXcd& K1, Eigen::MatrixXcd& K2, Eigen::MatrixXcd Iden);
	bool DeciTech(std::map<std::string, Eigen::MatrixXcd>& mM, std::map<std::string, Eigen::MatrixXcd> mMpre, std::map<std::string, Eigen::MatrixXcd> mFpre, std::map<std::string, std::pair<std::queue<Eigen::MatrixXcd>, std::queue<Eigen::MatrixXcd> > >& mq, const Eigen::MatrixXcd& vFrq, const double alpha, const size_t exFreq, const size_t nStep, size_t step);
private:
	v5cd mfcbylz;
	std::vector<double> omgv;
	double delta;
};

#endif
