#include <iostream>
#include <fstream>
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
	
}
