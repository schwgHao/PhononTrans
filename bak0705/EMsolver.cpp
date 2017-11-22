#include "EMsolver.h"
#include <fstream>
#include <iomanip>
#include <Dense>

using namespace Eigen;


void emsolver::DRCC(vector<vector<vector<vector<complex<double> > > > >& dccr, 
		            vector<vector<complex<double> > >& VibTRCOverk,
					vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadL, 
					vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadR,
		            const vector<vector<vector<vector<complex<double> > > > >& d00rl, 
					const vector<vector<vector<vector<complex<double> > > > >& d00rr){
	std::cout << string(45,'*')<<std::endl;
	std::cout << "Solve retarded phonon green's function." << std::endl;
	std::cout << string(45,'*')<<std::endl;
	int kcis = KCInd.size();
	int klmis = KLmInd.size();
	int klpis = KLpInd.size();
	for(int k = 0; k < kp.size(); k++){
		for(int o = 0; o < omg.size(); o++){
			Eigen::MatrixXcd Iden = Eigen::MatrixXcd::Identity(3*kcis, 3*kcis);
			complex<double> dtli = complex<double>(o, delta); 
			Iden *= dtli*dtli;
			Eigen::MatrixXcd D00l(3*klmis, 3*klmis);
			Eigen::MatrixXcd D00r(3*klpis, 3*klpis);
			for(int i = 0; i < 3*klmis; i++)
				for(int j = 0; j < 3*klmis; j++){
					D00l(i, j) = d00rl[k][o][i][j];
					D00r(i, j) = d00rr[k][o][i][j];
				}
			Eigen::MatrixXcd KC(3*kcis, 3*kcis);
			Eigen::MatrixXcd KLm(3*klmis, 3*kcis);
			Eigen::MatrixXcd KLp(3*kcis, 3*klpis);
			for(int i = 0; i < kcis; i++){
				for(int j = 0; j < kcis; j++)
				for(int ii = 0; ii < 3; ii++)
				for(int jj = 0; jj < 3; jj++){
					KC(3*i+ii, 3*j+jj) = mfc[k][3*KCInd[i]+ii][3*KCInd[j]+jj];
				}
				for(int j = 0; j < klmis; j++)
				for(int ii = 0; ii < 3; ii++)
				for(int jj = 0; jj < 3; jj++){
					KLm(3*j+jj, 3*i+ii) = mfc[k][3*KLmInd[j]+jj][3*KCInd[i]+ii];
					KLp(3*i+ii, 3*j+jj) = mfc[k][3*KLpInd[j]+jj][3*KCInd[i]+ii];
				}
			}
			Eigen::MatrixXcd SelfEngL = KLm.transpose()*D00l*KLm;
			Eigen::MatrixXcd SelfEngR = KLp*D00r*KLp.transpose();
			Eigen::MatrixXcd TC = (Iden - KC - SelfEngL - SelfEngR).fullPivLu().inverse();
			for(int i = 0; i < 3*kcis; i++)
				for(int j = 0; j < 3*kcis; j++){
					dccr[k][o][i][j] = TC(i, j);
					SelfEngLeadL[k][o][i][j] = SelfEngL(i, j);
					SelfEngLeadR[k][o][i][j] = SelfEngR(i, j);
				}

			Eigen::MatrixXcd vibT =(-1.0)*(SelfEngL-SelfEngL.adjoint())*TC*
							       (SelfEngR-SelfEngR.adjoint())*TC.adjoint();
			VibTRCOverk[k][o] = vibT.trace();
			if(omg.size()%(o+1) == 0) 
				std::cout << std::fixed<<std::setprecision(0)<<o*1.0/(omg.size()-1)*100 << " % completed." << std::endl;
		}
	}
}

void emsolver::writeTRC(const string& strs, vector<vector<complex<double> > >& vibTRCoverk){
//	double unitconv = hbar*sqrt(e/Ang2m/Ang2m/amu2kg); // eVs*sqrt(eV/(Ang^2*amu)) -> eV
	double unitconv = 0.064654147420;
	std::fstream ifs(strs, std::ios::out);
	ifs << "	Omega(eV)	" << "	" << "Vib. Trans. Coeff." << std::endl;  
	ifs << setiosflags(std::ios_base::fixed);
	for(int o = 0; o < omg.size(); o++){
		ifs << std::setprecision(10) << std::setw(15) << omg[o]*unitconv << "		";
		double tk = 0.0;
		for(int k = 0; k < kp.size(); k++)
			tk += vibTRCoverk[k][o].real();
		ifs << std::setw(15) << tk << std::endl;
	}
	ifs.close();
}
