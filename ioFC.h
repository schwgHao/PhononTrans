#ifndef IOFC_H_
#define IOFC_H_
#include <iostream>
#include <string>
#include <vector>
#include <complex>

using std::string;
using std::vector;
using std::complex;

struct iofc{
public:
//	iofc(string slabel_, int* nsc_, 
//			vector<vector<double> >kp_,
//			vector<vector<double> >xa_,
//			vector<int> xmass_,
//			vector<vector<double> >cell_,
//			vector<vector<double> >scell_) 
//		: slabel(slabel_), nsc(nsc_), 
//		  kp(kp_), xa(xa_), xmass(xmass_),
//			cell(cell_), scell(scell_){}
	iofc(string slabel_, bool isBulk_);
	void ReadFC(vector<vector<vector<complex<double> > > >& mfc);
	void ReadFCX(vector<vector<vector<complex<double> > > >& mfc);
	vector<vector<double> > kp2d(){return kp2D;}
	void kp2dinit();
	int Nat(){return na;}
	vector<vector<int> > MatCutbyLayer();
	void mfcSlice(vector<vector<vector<vector<vector<complex<double> > > > > >& mfcbylz, 
				  vector<vector<vector<complex<double> > > > mfc,
				  vector<vector<int> > atominLayer);
	double outermost(){return cell[2][2];}
	void setBoundary(double cllft, double clrt, int nlrt, vector<int>& KLm, vector<int>&KLp, vector<int>& KC);
    void setBoundary(int naAtOnelayer_l, int naAtOnelayer_r, vector<int>& KLm, vector<int>&KLp, vector<int>& KC);
public:
	string slabel;
	bool isBulk;
	int na;
	int ia1;
	int ia2;
	vector<int> nsc;
	vector<vector<double> > kp;
	vector<vector<double> > xa;
	vector<int> xmass;
	vector<vector<double> > cell;
	vector<vector<double> > scell;
	vector<vector<double> > kp2D;
};

typedef vector<vector<vector<vector<complex<double> > > > > v4cd;

struct setleads{
	double cellLz;
	int naAtOnelayer;
	vector<vector<double> > kp;
	v4cd d00r;
	setleads(string ssl);
};

#endif
