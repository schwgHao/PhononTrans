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
#endif
