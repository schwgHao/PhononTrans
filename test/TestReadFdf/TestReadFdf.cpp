#include "ReadFdf.h"

using namespace std;
int main(){
	string sname = "Cchain.fdf";
	fdf iofdf(sname);
	string res = iofdf.fdf_string("SystemName");
	cout << res << endl;
}
