#ifndef IOFC_H_
#define IOFC_H_
#include <iostream>
#include <string>
#include <vector>
#include <complex>

using std::string;
using std::vector;
using std::complex;

class iofc{
public:
	iofc(string slabel_, int* nsc_, 
			vector<vector<double> >kp_,
			vector<vector<double> >xa_,
			vector<int> xmass_,
			vector<vector<double> >cell_,
			vector<vector<double> >scell_) 
		: slabel(slabel_), nsc(nsc_), 
		  kp(kp_), xa(xa_), xmass(xmass_),
			cell(cell_), scell(scell_){}
	void ReadFC(vector<vector<vector<complex<double> > > >& mfc);
	void ReadFCX(vector<vector<vector<complex<double> > > >& mfc);
	vector<vector<double> > kp2d(){return kp2D;}
	vector<vector<int> > MatCutbyLayer();
	void mfcSlice(vector<vector<vector<vector<vector<complex<double> > > > > >& mfcbylz, 
				  vector<vector<vector<complex<double> > > > mfc,
				  vector<vector<int> > atominLayer);
private:
	string slabel;
	int* nsc;
	vector<vector<double> > kp;
	vector<vector<double> > xa;
	vector<int> xmass;
	vector<vector<double> > cell;
	vector<vector<double> > scell;
	vector<vector<double> > kp2D;
};


#endif
