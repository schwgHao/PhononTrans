#include <iostream>
#include <map>
#include <vector>
#include <Dense>
#include <fstream>


using namespace std;
using namespace Eigen;

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

    //double kp[1][4] = {0,0,0,1};
    vector<vector<double> > kp = {{0, 0, 0 ,1}};
    for(size_t i = 0; i < 4; ++i)
        cout << kp[0][i] << " ";
    cout << endl;
    MatrixXd A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    cout << A << endl;
    EigenSolver<MatrixXd> eigensolver(A);
    cout << eigensolver.eigenvalues() << endl;
    cout << eigensolver.eigenvectors() << endl;
    MatrixXcd B(3, 3);
    std::fstream ifs("a.txt", std::ios::in);
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            ifs >> B(i, j);
    cout << B << endl;
}
