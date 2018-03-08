#include "CalcPhEMT.h"
#include "FuncUtils.h"
#include "constants.h"
#include <cmath>

extern "C" void vibrator_c(char* slabel, bool Isbulk, 
						   double Enli, double Enlf, int Nenl, double delta);

void vibrator_c(char* slabel_, bool Isbulk,
		        double Enli, double Enlf, int Nenl, double delta){
	vector<string> s1 = strsplit((string)slabel_);
	string slabel = s1[0]; // the difference of string expression between fortran and C. 
	double unitconv = hbar*sqrt(e/Ang2m/Ang2m/amu2kg); // eVs*sqrt(eV/(Ang^2*amu)) -> s-1
	vector<double> omg(Nenl);
	for(int i = 0; i < Nenl; i++){
		if(Nenl == 1) {omg[0] = Enli/unitconv; break;}
		omg[i] = (Enli + 1.0*i/(Nenl-1)*(Enlf - Enli))/unitconv;
	}
	std::cout <<"Energy ranges from "<< Enli << " to ";
	std::cout << Enlf << " including "<<Nenl << " points."<< std::endl;	
	if(Isbulk)
		CalcPhLead(slabel, omg, delta);
	else
		CalcPhEMT(slabel, omg, delta);
	
}
