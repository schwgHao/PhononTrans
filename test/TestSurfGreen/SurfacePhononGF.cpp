#include "SurfacePhononGF.h"
#include "FuncUtils.h"
//#include "constants.h"
#include <fstream>
#include <iomanip>

using std::vector;
using std::complex;
using namespace Eigen;
using std::string;
using std::fstream;
using std::map;
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
		MatrixXcd K00(nx, nx);
		MatrixXcd K01(nx, nx);
		MatrixXcd K10(nx, nx);
		for(int ia = 0; ia < nx; ia++){
			for(int ja = 0; ja < nx; ja++){
				K00(ia, ja) = mfcbylz[k][0][0][ia][ja];
				K01(ia, ja) = mfcbylz[k][0][1][ia][ja];
				K10(ia, ja) = mfcbylz[k][1][0][ia][ja];
			}
		}
		for(int o = 0; o < no; o++){
			MatrixXcd K0 = K00;
			MatrixXcd K0s = K00;
			MatrixXcd K1 = K01;
			MatrixXcd K2 = K10;
			complex<double> dtli = complex<double>(omgv[o], delta);
			std::cout << "omega " << o << " : " << omgv[o] << std::endl;
			MatrixXcd vFrq = Eigen::MatrixXcd::Identity(nx, nx) * dtli * dtli;
			
			map<string, std::pair<queue<MatrixXcd>, queue<MatrixXcd> > > mq;
			queue<MatrixXcd> qr0s, qr0, qr1, qr2;
			queue<MatrixXcd> qf0s, qf0, qf1, qf2;
			
	//		mq.emplace("K0s", std::make_pair(qr0s, qf0s));
	//		mq.emplace("K0", std::make_pair(qr0, qf0));
			mq.emplace("K1", std::make_pair(qr1, qf1));
			mq.emplace("K2", std::make_pair(qr2, qf2));
			MatrixXcd mEpty = MatrixXcd::Zero(nx, nx);
			
			map<string, MatrixXcd> mMpre;
			mMpre.emplace("K0s", mEpty);
			mMpre.emplace("K0", mEpty);
			mMpre.emplace("K1", mEpty);
			mMpre.emplace("K2", mEpty);
			
			auto mFpre = mMpre;
			map<string, MatrixXcd> mM;
			mM.emplace("K0s", K0s);
			mM.emplace("K0", K0);
			mM.emplace("K1", K1);
			mM.emplace("K2", K2);
			const double alpha = 0.1;
			const size_t exFreq = 10000;
			const size_t nStep = 6;

			DeciTech(mM, mMpre, mFpre, mq, vFrq, alpha, exFreq, nStep, 0);
			K0s = mM["K0s"];
//			DeciTech(K0, K0s, K1, K2, vFrq);
			MatrixXcd d00rM = (vFrq - K0s).householderQr().solve(Eigen::MatrixXcd::Identity(nx, nx));
		//	d00rM = K01 * d00rM * K10;
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

void surfphGF::DeciTech(map<string, MatrixXcd>& mM, map<string, MatrixXcd> mMpre, map<string, MatrixXcd> mFpre, map<string, std::pair<queue<MatrixXcd>, queue<MatrixXcd> > >& mq, const MatrixXcd& vFrq, const double alpha, const size_t exFreq, const size_t nStep, size_t step)
{
	MatrixXcd K0s = mM["K0s"];
	MatrixXcd K0 = mM["K0"];
	MatrixXcd K1 = mM["K1"];
	MatrixXcd K2 = mM["K2"];
	const size_t nx = K0.rows();

	MatrixXcd d00 = (vFrq - K0).householderQr().solve(Eigen::MatrixXcd::Identity(nx, nx));
	
	MatrixXcd K0sNext = K1 * d00 * K2;

	double meps = K0sNext.norm();
	
	if(meps > 1.0e10) {
		std::cerr << "can not converge to finite values!" << std::endl;
		abort();
	}
	
	if(meps < 1.0e-8) return; 
	
	if(step % 100 == 0){
		std::cout << "step, d K0s: " << step + 1 << " : "<< meps << std::endl;
	}
	
//	if(step > 40){	
//		std::cout << "K0sNext " << K0sNext << std::endl;
//	}
	
	if(mq["K1"].first.size() == nStep)	{
		for(auto& i : mq) {
			i.second.first.pop();
			i.second.second.pop();
		}
	}
	
	auto K0Next = K2 * d00 * K1 + K1 * d00 * K2;

	map<string, MatrixXcd> mFi;
//	mFi.emplace("K0s", K0sNext);
//	mFi.emplace("K0", K2 * d00 * K1 + K1 * d00 * K2);
	mFi.emplace("K1", K1 * d00 * K1 - K1);
	mFi.emplace("K2", K2 * d00 * K2 - K2);
//	MatrixXcd fk0s = K0sNext;
	
	
	for(auto& i : mq){
		string ss = i.first;
		auto& dk = i.second.first;
		dk.push(mM[ss] - mMpre[ss]);
		auto& df = i.second.second;
		df.push(mFi[ss] - mFpre[ss]);
	}

	const size_t sq = mq["K1"].first.size();
	const size_t sf = nx * sq;
	
	map<string, MatrixXcd> mR;
	map<string, MatrixXcd> mF;
	for(auto s : mq){
		MatrixXcd Ri(nx, sf);
		MatrixXcd Fi(nx, sf);
		for(size_t i = 0; i < sq; ++i) {
			Ri.block(0, i*nx, nx, nx) = s.second.first.front();
			s.second.first.pop();
			Fi.block(0, i*nx, nx, nx) = s.second.second.front();
			s.second.second.pop();
		}
		mR.emplace(s.first, Ri);	
		mF.emplace(s.first, Fi);	
	}

	map<string, MatrixXcd> mC;
	for(auto s : mq){
		mC.emplace(s.first, alpha * Eigen::MatrixXcd::Identity(nx, nx));
	}

	mMpre = mM;
	for(auto i : mq) {
		string s = i.first;
		if(++step % exFreq == 0) {
			auto& fi = mF[s];
			mC[s] = mC[s] - (mR[s] + alpha * fi)*(fi.transpose() * fi).householderQr().solve(Eigen::MatrixXcd::Identity(sf, sf)) * fi.transpose();
			std::cout << "fiTfi " << fi.transpose() * fi << std::endl;
			std::cout << s << " Ci " << mC[s] << std::endl;

		}
		mM[s] = mM[s] + mC[s] * mFi[s];
	}
	
	mM["K0s"] = mM["K0s"] + K0sNext; 
	mM["K0"] = mM["K0"] + K0Next; 

//	DeciTech(mM, mMpre, mFi, K2, Iden);
	DeciTech(mM, mMpre, mFi, mq, vFrq, alpha, exFreq, nStep, step);
}

void surfphGF::DeciTech(MatrixXcd& K0, MatrixXcd& K0s, 
		                MatrixXcd& K1, MatrixXcd& K2, MatrixXcd Iden)
{	
	const size_t nx = Iden.cols();

	MatrixXcd d00 = (Iden - K0).householderQr().solve(Eigen::MatrixXcd::Identity(nx, nx));
	
	MatrixXcd K0sNext = K1 * d00 * K2;

	double meps = K0sNext.norm();
	if(meps > 1.0e10) {
		std::cerr << "can not converge to finite values!" << std::endl;
		abort();
	}
	if(K0sNext.norm() < 1.0e-8) return; 
	
	K0s = K0s + 0.1*K0sNext;
	K0 = K0 + 0.1 *(K2 * d00 * K1 + K1 * d00 * K2);
	K1 = K1 + 0.1*(K1 * d00 * K1 - K1);
	K2 = K2 + 0.1*(K2 * d00 * K2 - K2);
	DeciTech(K0, K0s, K1, K2, Iden);
}

void surfphGF::DR00(MatrixXcd t0, MatrixXcd t0t, MatrixXcd Tphnext, MatrixXcd& Tph){
	MatrixXcd Iden = MatrixXcd::Identity(t0.rows(), t0.cols());
	MatrixXcd t0next = (Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0 * t0;
	MatrixXcd t0tnext = (Iden - t0*t0t - t0t*t0).householderQr().solve(Iden) * t0t * t0t;
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
