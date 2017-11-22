#ifndef EMSOLVER_H_
#define EMSOLVER_H_

#include <iostream>
#include <vector>
#include <string>
#include <complex>

using std::vector;
using std::string;
using std::complex;

//typedef vector<vector<double> > vector<vector<double> >;
//typedef vector<vector<complex<double> > > vector<vector<complex<double> > >;
//typedef vector<vector<vector<vector<complex<double> > > > > vector<vector<vector<vector<complex<double> > > > >;

class emsolver{
public:
	emsolver(const vector<vector<vector<complex<double> > > >& mfc_,
			vector<int> klmi_, vector<int> klpi_, vector<int> kci_, 
			vector<vector<double> > kp_,
			vector<double> omg_, double delta_):
			KLmInd(klmi_), KLpInd(klpi_), KCInd(kci_),
			kp(kp_), omg(omg_), delta(delta_){};
	void DRCC(vector<vector<vector<vector<complex<double> > > > >& dccr, 
			  vector<vector<complex<double> > >& VibTRCOverk,
			  vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadL, 
			  vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadR,
			  const vector<vector<vector<vector<complex<double> > > > >& d00rl, 
			  const vector<vector<vector<vector<complex<double> > > > >& d00rr);
	
	void writeTRC(const string& strs, vector<vector<complex<double> > >& vibTRCoverk);
private:
	vector<vector<vector<complex<double> > > > mfc;
	vector<int> KLmInd;
	vector<int> KLpInd;
	vector<int> KCInd;
	vector<vector<double> > kp;
	vector<double> omg;
	double delta;
};


#endif
