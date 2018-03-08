#ifndef FUNCUTIL_H_
#define FUNCUTIL_H_

#include <iostream>
#include <string>
#include <vector>

using std::string;
using std::vector;

vector<string> strsplit(const string& str, string pat=" ");

string strToupper(const string& str);

string strTolower(const string& str);

template<typename T>
void mcopy(T* a, size_t n1, size_t n2, vector<vector<T> >& vec){	
	vector<T> tp;
	for(int i = 0; i < n1; i++){
		for(int j = 0; j < n2; j++)
			tp.push_back(*(a + j + n2*i));
		vec.push_back(tp);
		tp.clear();
	}
//		copy(&a[i][0], &a[i][0]+n2, back_inserter(vec[i]));
}


#endif
