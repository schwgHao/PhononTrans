#include "CalcPhEMT.h"
#include "FuncUtils.h"

extern "C" void vibrator_c(char* slabel, bool Isbulk);

void vibrator_c(char* slabel_, bool Isbulk){
	vector<string> s1 = strsplit((string)slabel_);
	string slabel = s1[0]; // the difference of string expression between fortran and C. 
	
	vector<double> omg(20);
	for(int i = 0; i < omg.size(); i++)
		omg[i] = 1.0e-4 + 1.0e-1*i/omg.size();
	double delta = 1.0e-6;
	
	if(Isbulk)
		CalcPhLead(slabel, omg, delta);
	else
		CalcPhEMT(slabel, omg, delta);
	
}
