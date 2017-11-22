#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct doublessthan{
	bool operator()(double l, double r){
		return (r-l > 1.0e-5);
	}
};

int main(){
	
	vector<double> pos={1.0, 1.4, 2.0, 3.1, 2.0};
	
	map<double, vector<int>, doublessthan> posbylz;
	
	for(int i = 0; i < pos.size(); i++){
		posbylz[pos[i]].push_back(i);
	}
	auto i = posbylz.find(3.1);
	cout << ((*i).second)[0] << endl;
	i = posbylz.begin();
	i++;
	i++;
	i++;
	cout << (*i).first << endl;

	vector<double> pos1;
	pos1 = pos;
	for(auto i : pos1) cout << i << endl;
}
