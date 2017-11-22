#include "SurfacePhononGF.h"
#include "FuncUtils.h"
#include "constants.h"
#include <fstream>

using std::fstream;

void surfphGF::DR00(v4cd& d00r){
	int nk = mfcbylz.size();
	int no = omgv.size();
	int nx = mfcbylz[0][0][0].size();
	for(int k = 0; k < nk; k++){
		for(int o = 0; o < no; o++){
			MatrixXcd T0(nx, nx);
			MatrixXcd T0a(nx, nx);
			MatrixXcd T0tilt(nx, nx);
			MatrixXcd T00(nx, nx);
			MatrixXcd T01(nx, nx);
			MatrixXcd Iden = MatrixXcd::Identity(nx, nx);
			complex<double> dtli = complex<double>(o, delta); 
			Iden *= dtli*dtli;
			for(int ia = 0; ia < nx; ia++)
				for(int ja = 0; ja < nx; ja++){
					T00(ia, ja) = mfcbylz[k][0][0][ia][ja];
					T01(ia, ja) = mfcbylz[k][0][1][ia][ja];
				}
			T0a = Iden - T00;
			T0a = T0a.fullPivLu().inverse();
			T0 = T0a*T01.adjoint();
			T0tilt = T0a*T01;
			MatrixXcd Tph = T0tilt;
			DR00(T0, T0tilt, Tph);
			MatrixXcd d00rM = (Iden - T00 - T01*Tph).fullPivLu().inverse();
			d00rM *= hbar/2.0/pi;
			for(int ia = 0; ia < nx; ia++)
				for(int ja = 0; ja < nx; ja++)
					d00r[k][o][ia][ja] = d00rM(ia, ja);

		}
	}
}

void surfphGF::DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd& Tph){
	MatrixXcd Iden = MatrixXcd::Identity(t0.rows(), t0.cols());
	MatrixXcd t0next = (Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0*t0;
	MatrixXcd t0tnext = (Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0t*t0t;
	MatrixXcd Tphnext = Tph*t0tnext;
	if(Tphnext.norm() < 1.0e-8) return; 
	else{
		Tph += Tphnext;
		DR00(t0next, t0tnext, Tph);
	}
}

void surfphGF::writeDR00(string ssgf, double cellLz, int nlayers, 
						 const vector<vector<double> >& kp2d, 
				         const v4cd& d00r){
	fstream ofs(ssgf, std::ios::out|std::ios::binary);
	ofs.write((char*)&cellLz, sizeof(double));
	ofs.write((char*)&nlayers, sizeof(int));
	for(int i = 0; i < kp2d.size(); i++)
		ofs.write((char*)&kp2d[i], kp2d.size()*sizeof(double));
	for(int i = 0; i < d00r.size(); i++)
		for(int j = 0; j < d00r[0].size(); j++)
			for(int k = 0; k < d00r[0][0].size(); k++)
				ofs.write((char*)(&d00r[i][j][k][0]), d00r[0][0][0].size()*sizeof(complex<double>));
//	double kp2da[kp2d.size()][kp2d[0].size()];
//	mcopy(kp2d, (double*)kp2da);
//	ofs.write(reinterpret_cast<char*>(kp2da), sizeof(kp2da));
	ofs.close();
}
