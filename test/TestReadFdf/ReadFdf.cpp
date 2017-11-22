#include "ReadFdf.h"
#include "FuncUtils.h"
#include <vector>
#include <cctype>
#include <algorithm>
#include <fstream>

using namespace std;

string fdf::fdf_string(string str){
	fstream fs(sname);
	string res="";
	string s1;
	while(getline(fs, s1)){
		vector<string> ss = strsplit(s1);
		string strs = strToupper(ss[0]);
		string strt = strToupper(str);
		if(strs == strt){res = ss[1]; break;}
	}
	fs.close();
	return res;
}
