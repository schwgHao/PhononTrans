#include "SurfacePhononGF.h"
#include "FuncUtils.h"
//#include "constants.h"
#include <fstream>
#include <iomanip>

using std::fstream;
using std::vector;
using std::complex;
using std::string;
using std::queue;

void surfphGF::DR00(v4cd& d00r){
	std::cout << string(45,'*')<<std::endl;
	std::cout << "Solve surface phonon green's function iteratively." << std::endl;
	std::cout << string(45,'*')<<std::endl;
	int nk = mfcbylz.size();
	int no = omgv.size();
//	std::cout << "surfphGF: omgv: " << std::endl;
//	for(auto i : omgv) std::cout << i << std::endl;
//	std::cout << "no " << no << std::endl;
//	std::cout << "nk " << nk << std::endl;
	int nx = mfcbylz[0][0][0].size();
//	std::cout << "nx: " << nx << std::endl;
	for(int k = 0; k < nk; k++){
		MatrixXcd T00(nx, nx);
		MatrixXcd T01(nx, nx);
		for(int ia = 0; ia < nx; ia++){
			for(int ja = 0; ja < nx; ja++){
				T00(ia, ja) = mfcbylz[k][0][0][ia][ja];
				T01(ia, ja) = mfcbylz[k][0][1][ia][ja];
			}
		}
		for(int o = 0; o < no; o++){
			MatrixXcd T0(nx, nx);
			MatrixXcd T0a(nx, nx);
			MatrixXcd T0tilde(nx, nx);
			MatrixXcd Tph(nx, nx);
			std::cout << "o, omg[o] " << o << " " << omgv[o] << std::endl;
			complex<double> dtli = complex<double>(omgv[o], delta);
//			complex<double> dtli = complex<double>(0.5, delta);
//			std::cout << "o, omg[o], delta: " << o << " " << omgv[o] << " " << dtli << std::endl;
			MatrixXcd Iden = Eigen::MatrixXcd::Identity(nx, nx) * dtli * dtli;
			T0a = (Iden - T00).householderQr().solve(Eigen::MatrixXcd::Identity(nx, nx));
			T0 = T0a*T01.adjoint();
	//		std::cout << "T0 " << T0 << std::endl;
			T0tilde = T0a*T01;
			Tph = T0tilde;

			const size_t exFreq = 10000;
			const size_t nStep = 1;
			const size_t alpha = 0.1;
			
			queue<MatrixXcd> qt1, qt2, qf1, qf2; 
			
			MatrixXcd mEpty = MatrixXcd::Zero(nx, nx);
			DR00(T0, T0tilde, T0tilde, Tph);
//			DR00(T0, T0tilde, mEpty, mEpty, mEpty, mEpty, T0tilde, Tph, alpha, exFreq, nStep, 0, qt1, qt2, qf1, qf2);
	//		MatrixXcd d00rM = (Iden - T00 - T01*Tph).fullPivLu().inverse();
			MatrixXcd d00rM = (Iden - T00 - T01*Tph).householderQr().solve(Eigen::MatrixXcd::Identity(nx, nx));
	//		d00rM *= hbar/2.0/pi;
			for(int ia = 0; ia < nx; ia++){
				for(int ja = 0; ja < nx; ja++){
					d00r[k][o][ia][ja] = d00rM(ia, ja);
				}
			}
	//		if(no%(o+1) == 0) 
	//			std::cout << std::fixed<<std::setprecision(0) << o*1.0/(no-1)*100 << " % completed." << std::endl;
		}
	}
}

void surfphGF::DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd t1pre, MatrixXcd t2pre, MatrixXcd fpre, MatrixXcd ftpre, MatrixXcd Tphnext, MatrixXcd& Tph, const double alpha, const size_t exFreq, const size_t nStep, size_t step, queue<MatrixXcd>& qt1, queue<MatrixXcd>& qt2, queue<MatrixXcd>& qf1, queue<MatrixXcd>& qf2)
{
	const size_t nr = t0.rows();
	const size_t nc = t0.cols();
	MatrixXcd Iden = MatrixXcd::Identity(nr, nc);
//	MatrixXcd t0next = 0.01*((Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0*t0-t0);
//	MatrixXcd t0tnext = 0.01*((Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0t*t0t-t0t);
	MatrixXcd t0next = (Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0 * t0;
	MatrixXcd t0tnext = (Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0t * t0t;

	
	if(qf1.size() == nStep) {
		qf1.pop();
		qf2.pop();
		qt1.pop();
		qt2.pop();
	}

	auto fnext = t0next - t0;
	auto ftnext = t0tnext - t0t;

	qf1.push(fnext - fpre);
	qf2.push(ftnext - ftpre);
	qt1.push(t0 - t1pre);
	qt2.push(t0t - t2pre);
	
	const size_t sq = qf1.size();
	const size_t sf = nc * sq;
	MatrixXcd Ri1(nr, sf);
	MatrixXcd Ri2(nr, sf);
	MatrixXcd Fi1(nr, sf);
	MatrixXcd Fi2(nr, sf);
	auto cqt1 = qt1;
	auto cqt2 = qt2;
	auto cqf1 = qf1;
	auto cqf2 = qf2;
	for(size_t i = 0; i < sq; ++i) {
		Ri1.block(0, i*nc, nc, nc) = cqt1.front();
		cqt1.pop();
		Ri2.block(0, i*nc, nc, nc) = cqt2.front();
		cqt2.pop();
		Fi1.block(0, i*nc, nc, nc) = cqf1.front();
		cqf1.pop();
		Fi2.block(0, i*nc, nc, nc) = cqf2.front();
		cqf2.pop();
	}

	MatrixXcd C1 = alpha*Eigen::MatrixXcd::Identity(nr, nc);
	MatrixXcd C2 = alpha*Eigen::MatrixXcd::Identity(nr, nc);
	if ((step + 1) % exFreq == 0) {
			C1 = C1 - (Ri1 + alpha*Fi1)*(Fi1.transpose()*Fi1).householderQr().solve(Eigen::MatrixXcd::Identity(sf, sf))*Fi1.transpose();
			C2 = C2 - (Ri2 + alpha*Fi2)*(Fi2.transpose()*Fi2).householderQr().solve(Eigen::MatrixXcd::Identity(sf, sf))*Fi2.transpose();
	}

	t0next = t0 + C1 * fnext; 
	t0tnext = t0t + C2 * ftnext; 
	
	MatrixXcd Tphtnext = Tphnext*t0tnext;
	
	
	double meps= Tphtnext.norm();
	
	if(meps > 1.0e10){
		std::cerr << "Surface Green function can not converge to a finite vaule." << std::endl;
		abort();
	}

	if(meps < 1.0e-8) return; 
	
	if(step % 7 == 0){
		std::cout << "step: " << step << std::endl;
		std::cout << "Tphtnext: "<< std::endl;
		std::cout << Tphtnext << std::endl;
	}

	Tph += Tphtnext;
	
	DR00(t0next, t0tnext, t0, t0t, fnext, ftnext, Tphtnext, Tph, alpha, exFreq, nStep, ++step, qt1, qt2, qf1, qf2);
}

void surfphGF::DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd Tphnext, MatrixXcd& Tph){
	MatrixXcd Iden = MatrixXcd::Identity(t0.rows(), t0.cols());
//	MatrixXcd t0next = 0.01*((Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0*t0-t0);
//	MatrixXcd t0tnext = 0.01*((Iden - t0*t0t - t0t*t0).fullPivLu().inverse()*t0t*t0t-t0t);
	MatrixXcd t0next = t0 + 1 * ((Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0 * t0 - t0);
	MatrixXcd t0tnext = t0t + 1 * ((Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0t * t0t - t0t);
	MatrixXcd Tphtnext = Tphnext*t0tnext;
	std::cout << "Tphtnext: "<< std::endl;
	std::cout << Tphtnext << std::endl;
	if(Tphnext.norm() < 1.0e-8) return; 
	else{
		Tph += Tphtnext;
		DR00(t0next, t0tnext, Tphtnext, Tph);
	}
}


void surfphGF::writeDR00(string ssgf, double cellLz, int nlayers, 
						 const vector<vector<double> >& kp2d, 
				         const v4cd& d00r){
	fstream ofs(ssgf, std::ios::out|std::ios::binary);
	ofs.write((char*)&cellLz, sizeof(double));
	ofs.write((char*)&nlayers, sizeof(int));
	int kpsz1 = kp2d.size();
	int kpsz2 = kp2d[0].size();
	ofs.write((char*)&kpsz1, sizeof(int));
	ofs.write((char*)&kpsz2, sizeof(int));
	for(int i = 0; i < kpsz1; i++)
		ofs.write((char*)&kp2d[i][0], kpsz2*sizeof(double));
	
	int dsz1 = d00r.size();
	int dsz2 = d00r[0].size();
	int dsz3 = d00r[0][0].size();
	int dsz4 = d00r[0][0][0].size();
//	std::cout << "d00r size: "<<dsz1<<"*"<<dsz2<<"*"<<dsz3<<"*"<<dsz4<<"*"<<std::endl;
	ofs.write((char*)&dsz1, sizeof(int));
	ofs.write((char*)&dsz2, sizeof(int));
	ofs.write((char*)&dsz3, sizeof(int));
	ofs.write((char*)&dsz4, sizeof(int));
	for(int i = 0; i < dsz1; i++)
		for(int j = 0; j < dsz2; j++)
			for(int k = 0; k < dsz3; k++)
				ofs.write((char*)(&d00r[i][j][k][0]), dsz4*sizeof(complex<double>));
//	double kp2da[kp2d.size()][kp2d[0].size()];
//	mcopy(kp2d, (double*)kp2da);
//	ofs.write(reinterpret_cast<char*>(kp2da), sizeof(kp2da));
	ofs.close();
}
