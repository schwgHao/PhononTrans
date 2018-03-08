#include "CalcPhEMT.h"
#include "FuncUtils.h"
#include "constants.h"
#include <cmath>
//#include <iomanip>

extern "C" void vibrator_c(char* slabel, bool Isbulk, 
						   double Enli, double Enlf, int Nenl, double delta);

void vibrator_c(char* slabel_, bool Isbulk,
		        double Enli, double Enlf, int Nenl, double delta){
	vector<string> s1 = strsplit((string)slabel_);
	string slabel = s1[0]; // the difference of string expression between fortran and C. 
//	std::cout << std::fixed << std::setprecision(12) << unitconv << std::endl;
//	double unitconv = hbar*sqrt(e/Ang2m/Ang2m/amu2kg); // eVs*sqrt(eV/(Ang^2*amu)) -> s-1
	double unitconv = 0.064654147420;
	vector<double> omg(Nenl);
	for(int i = 0; i < Nenl; i++){
		if(Nenl == 1) {omg[0] = Enli/unitconv; break;}
		omg[i] = (Enli + 1.0*i/(Nenl-1)*(Enlf - Enli))/unitconv;
	}
//	std::cout <<"Energy range from " << Isbulk <<" ";
//	std::cout << Enli << " to " << Enlf << " totally "<<Nenl << "	" << delta << std::endl;	
	if(Isbulk)
		CalcPhLead(slabel, omg, delta);
	else
		CalcPhEMT(slabel, omg, delta);
	
}
