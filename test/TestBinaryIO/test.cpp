#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

int main(){
	double arr[3][4] ={{1,2,3,4},{5,6,7,8},{9,10,11,12}};
	fstream ofc("test.txt", std::ios::out|std::ios::binary);
	ofc.write((char*)arr, sizeof(arr));
	ofc.close();
	double xy[3][4];
	fstream ifc("test.txt", std::ios::in|std::ios::binary);
	ifc.read((char*)xy, sizeof(xy));
	ifc.close();
	for(size_t i = 0; i < 3; ++i){
		for(size_t j = 0; j < 4; ++j){
			cout << xy[i][j]<<"	";
		}
		cout << endl;
	}
}
