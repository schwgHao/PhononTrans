#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <complex>
using namespace std;

void func(double* arr, int n1, int n2){
	for(int i = 0; i < n1; i++){
		for(int j = 0; j < n2; j++)
			cout << *(arr + j + i*n2) << " ";
//			cout << arr[i][j] << " ";
		cout << endl;
	}
}

int main(){
	fstream ifs("aa.test", ios::out|ios::binary);
	double a[2][3];
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			a[i][j] = (i+1)*10 + j+1;
	ifs.write(reinterpret_cast<char*>(a), sizeof(a));
	ifs.close();

	fstream ofs("aa.test", ios::in|ios::binary);
	double ac[2][3];
	ofs.read(reinterpret_cast<char*>(ac), sizeof(ac));
	ofs.close();
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 3; j++)
			cout << &(ac[i][j]) << " "<< ac[i][j] << " ";
		cout << endl;
	}
	func((double*)ac, 2, 3);
	/////////////////////////
	//
	vector<vector<complex<double> > > avi(2, vector<complex<double> >(3));
	srand(time(0));
//	for(auto i : avi)
//		fill((i.begin()), (i.end()), 1.0);
//		fill((i.begin()), (i.end()), rand()*1.0/RAND_MAX);
	for(int i = 0; i < avi.size(); i++)
		for(int j = 0; j < avi[0].size(); j++)
			avi[i][j] = complex<double>(rand()*1.0/RAND_MAX, 1.0);
	
	for(auto i: avi){
		for(auto j : i) cout << j << " ";
		cout << endl;
	}
	fstream ofsa("ab.test", ios::out|ios::binary);
	ofsa.write((char*)&avi[0], avi.size()*sizeof(complex<double>));
	ofsa.write((char*)&avi[1], avi.size()*sizeof(complex<double>));
	ofsa.close();

	/*read*/

	vector<vector<complex<double> > > avo(2, vector<complex<double> >(3));
	fstream ifsa("ab.test", ios::in|ios::binary);
	ifsa.read((char*)&avo[0], avo.size()*sizeof(complex<double>));
	ifsa.read((char*)&avo[1], avo.size()*sizeof(complex<double>));
	ifsa.close();
	
	for(auto i: avo){
		for(auto j : i) cout << j << " ";
		cout << endl;
	}
}
