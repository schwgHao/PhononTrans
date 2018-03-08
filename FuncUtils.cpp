#include "FuncUtils.h"


vector<string> strsplit(const string& str, string pat){
	vector<string> res;
	int pos1 = 0, pos2 = str.find(pat, pos1);
	if(pos2 == string::npos){
		res.push_back(str);
		return res;
	}
	while(pos2 != string::npos){
		string s1 = str.substr(pos1, pos2-pos1);
		if(s1 != pat && !s1.empty()) res.push_back(s1);
		pos1 = pos2 + pat.size();
		pos2 = str.find(pat, pos1);
	}
	if(pos1 != str.size()) res.push_back(str.substr(pos1));
	return res;
}

string strToupper(const string& str){
	string s1 = str;
	for(auto i = s1.begin(); i != s1.end(); i++)
		if(*i >= 'a' && *i <= 'z') *i = *i-('a' - 'A');
	return s1;
}

string strTolower(const string& str){
	string s1 = str;
	for(auto i = s1.begin(); i != s1.end(); i++)
		if(*i >= 'A' && *i <= 'Z') *i = *i+('a' - 'A');
	return s1;
}
