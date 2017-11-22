#include "CalcPhEMT.h"
#include "FuncUtils.h"
#include "constants.h"
#include <cmath>

extern "C" void vibrator_c(char* slabel, bool Isbulk, 
						   double Enli, double Enlf, int Nenl, double delta,
						   double bathTemp1, double bathTemp2, int Ntemp);

void vibrator_c(char* slabel_, bool Isbulk,
		        double Enli, double Enlf, int Nenl, double delta, 
		        double bathTemp1, double bathTemp2, int Ntemp){
	vector<string> s1 = strsplit((string)slabel_);
	string slabel = s1[0]; // the difference of string expression between fortran and C. 
	double unitconv = hbar*sqrt(e/Ang2m/Ang2m/amu2kg); // eVs*sqrt(eV/(Ang^2*amu)) -> eVs*s-1
	vector<double> omg(Nenl);
	std::cout << "omga: " << std::endl;
	for(int i = 0; i < Nenl; i++){
		if(Nenl == 1) {omg[0] = Enli/unitconv; break;}
		omg[i] = (Enli + 1.0*i/(Nenl-1)*(Enlf - Enli))/unitconv;
		std::cout << omg[i] << std::endl;
	}

	vector<double> bathT(Ntemp);
	for(int i = 0; i < Ntemp; i++){
		if(Ntemp == 1){bathT[0] = bathTemp1; break;}
		bathT[i] = bathTemp1 + 1.0*i/(Ntemp-1)*(bathTemp2 - bathTemp1);
	}
	std::cout <<"Temperature ranges from "<< bathTemp1 << " to ";
	std::cout << bathTemp2 << " including "<< Ntemp << " points."<< std::endl;

	if(Isbulk)
		CalcPhLead(slabel, omg, delta);
	else
		CalcPhEMT(slabel, omg, delta, bathT);

	std::cout << "return from vibrator_c" << std::endl;
	
}
