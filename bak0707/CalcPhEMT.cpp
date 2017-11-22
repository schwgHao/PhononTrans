#include "CalcPhEMT.h"
#include "ioFC.h"
#include "EMsolver.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

using std::cerr;
using std::cout;
using std::endl;
bool doubpredicate(double i, double j){
	return (fabs(i - j) < 1.0e-6);	
}

void CalcPhEMT(string slabel, const vector<double>& omg, double delta, const vector<double>& bathT){
	iofc ifc(slabel, false);
	int natoms = ifc.Nat();
/* Read force constant matrix */	
	vector<vector<vector<complex<double> > > > mfc(1,
				 vector<vector<complex<double> > >(3*natoms, 
						  vector<complex<double> >(3*natoms, 0.0)));

//	iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);
	ifc.ReadFC(mfc);
	
	setleads Llft("bulksgf.lft");
	setleads Lrt("bulksgf.rt");
	
	double cllft = Llft.cellLz;
	double clrt = Lrt.cellLz;
	int nlrt = Lrt.nlayers;
	vector<vector<double> > kp2d = ifc.kp2d();
	vector<vector<double> > kplft = Llft.kp;
	vector<vector<double> > kprt = Lrt.kp;
	
	if(kp2d.size() != kplft.size() || kp2d.size() != kprt.size())
	{
		cout << "Error: different No. of k point in leads and EM.\n";
		abort();
	}
//	for(int i = 0; i < kp2d.size(); i++){
//		if(!equal(kp2d[i].begin(), kp2d[i].end(), kplft[i].begin(), doubpredicate) || 
//		   !equal(kp2d[i].begin(), kp2d[i].end(), kplft[i].begin(), doubpredicate)){
//			cout << "Error: different k point in leads and EM.\n";
//			abort();
//		}
//	}
	
	v4cd d00rl = Llft.d00r;
	v4cd d00rr = Llft.d00r;
	vector<int> KLmInd, KLpInd, KCInd;
	ifc.setBoundary(cllft, clrt, nlrt, KLmInd, KLpInd, KCInd);

	emsolver ems(mfc, KLmInd, KLpInd, KCInd, kp2d, omg, delta, bathT);
	vector<vector<vector<vector<complex<double> > > > > dccr(kp2d.size(),
		   vector<vector<vector<complex<double> > > > (omg.size(),
		          vector<vector<complex<double> > > (3*KCInd.size(),
				         vector<complex<double> > (3*KCInd.size(), 0.0))));

	auto SelfEngLeadL = dccr;
	auto SelfEngLeadR = dccr;
	
	vector<vector<complex<double> > > vibTRC(kp2d.size(),
		   vector<complex<double> > (omg.size(), 0.0));
		
	ems.DRCC(dccr, vibTRC, SelfEngLeadL, SelfEngLeadR, d00rl, d00rr);

	std::fstream ifs("VibTRC.dat", std::ios::out);
	
	for(auto i : vibTRC){
		for(auto j : i){
			ifs << std::setw(14) << std::fixed << std::scientific  << j << " ";
		}
		ifs << std::endl;
	}
	
	ifs.close();

	string strs = slabel + ".PHTRC";
//	cout << string(45,'-') << endl;
//	cout << "		Write phonon TRC." << endl;	
//	cout << string(45,'-') << endl;
	ems.writePhTRC(strs, vibTRC);
}
