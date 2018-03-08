#ifndef READFDF_H_
#define READFDF_H_

#include <iostream>
#include <string>

using std::string;
class fdf{
public:
	fdf(string sname_): sname(sname_){};
	string fdf_string(string str);
private:
	string sname;
};
#endif
