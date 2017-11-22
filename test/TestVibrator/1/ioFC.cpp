#include "ioFC.h"
#include <fstream>
#include <sstream> 
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>

using std::fstream;
using std::stringstream;
using std::set;
using std::pair;
using std::cout;
using std::endl;

/* 
 * Read force from label.FC to construct force constant matrix.
 *
 * mfc[nk][3*na][3*na] has 3 dim.
 * nk is the number of 2d kpoints in transversal direction.
 * 
 * The private variable kp2D is initilized.
 * 2017/6/8, Hao Wang
 * */
void iofc::ReadFC(vector<vector<vector<complex<double> > > >& mfc){
	int na = xa.size();
	string sn = slabel + ".FC";
	vector<vector<vector<vector<vector<double> > > > > pp(nsc[0], 
	 	   vector<vector<vector<vector<double> > > >(nsc[1], 
				  vector<vector<vector<double> > >(nsc[2], 
	                     vector<vector<double> >(na,
							    vector<double>(3, 0.0)))));
	auto pn = pp;
	
	vector<vector<vector<vector<vector<vector<vector<double> > > > > > > phi(na, 
		   vector<vector<vector<vector<vector<vector<double> > > > > > (3, 
				  vector<vector<vector<vector<vector<double> > > > >(nsc[0], 
						 vector<vector<vector<vector<double> > > >(nsc[1], 
								vector<vector<vector<double> > >(nsc[2], 
									   vector<vector<double> >(na,
											  vector<double>(3, 0.0)))))));

	auto phi0 = phi;
	
	fstream ifs(sn);
	string s1;
	getline(ifs, s1);
	stringstream ist;
	for(int i = 0; i < na; i++){
		for(int xy1 = 0; xy1 < 3; xy1++){
			//read positive displacement
			for(int lx = 0 ; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++)
			for(int lz = 0; lz < nsc[2]; lz++)
				for(int j = 0; j < na; j++){
					getline(ifs, s1);
					ist.clear();
					ist.str(s1);
					for(int xy2 = 0; xy2 < 3; xy2++){
						ist >> pp[lx][ly][lz][j][xy2];
					}
				}
			}// end read positive displ.
			// read negative displ.
			for(int lx = 0; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++)
			for(int lz = 0; lz < nsc[2]; lz++)
				for(int j = 0; j < na; j++){
					getline(ifs, s1);
					ist.clear();
					ist.str(s1);
					for(int xy2 = 0; xy2 < 3; xy2++){
						ist >> pn[lx][ly][lz][j][xy2];
					}
				}
			}// end read negative displ.
			// average displ.
			for(int lx = 0; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++)
			for(int lz = 0; lz < nsc[2]; lz++)
				for(int j = 0; j < na; j++)
					for(int xy2 = 0; xy2 < 3; xy2++){
						phi0[i][xy1][lx][ly][lz][j][xy2] = 
							0.5*(pp[lx][ly][lz][j][xy2] +
							     pn[lx][ly][lz][j][xy2]);//average displ.
					}
			}// end average displ.
		}
	}
	ifs.close();
	//make force constant matrix hermitian.
	for(int i = 0; i < na; i++){
		for(int ii = 0; ii < 3; ii++)
			for(int j = 0; j < na; j++)
				for(int jj = 0; jj < 3; jj++)
					for(int lx = 0 ; lx < nsc[0]; lx++)
					for(int ly = 0; ly < nsc[1]; ly++)
					for(int lz = 0; lz < nsc[2]; lz++){
						phi[i][ii][lx][ly][lz][j][jj] =
							0.5 * ( phi0[i][ii][lx][ly][lz][j][jj] +
							phi0[i][ii][nsc[0]-lx-1][nsc[1]-ly-1][nsc[2]-ly-1][j][jj]);
					}
	}


	/* Now form zero(ix, jx, j)*/
	int ncells = 1;
	for(int i = 0; i < 3; i++)
		ncells *= nsc[i];

	vector<vector<vector<double> > > zero(na, 
		   vector<vector<double> >(3,
				  vector<double>(3, 0.0)));
	for(int j = 0; j < na; j++){
		for(int ii = 0; ii < 3; ii++)
			for(int ij = 0; ij < 3; ij++){
				for(int i = 0; i < na; i++)
					for(int lx = 0; lx < nsc[0]; lx++)
					for(int ly = 0; ly < nsc[1]; ly++)
					for(int lz = 0; lz < nsc[2]; lz++)
						zero[j][ii][ij] += phi[j][ii][lx][ly][lz][i][ij];
				zero[j][ii][ij] /= (na*ncells);
			}
	}
	/*Now form zeroo(ix, jx), the sum over j of zero(j, ix, jx)*/
	vector<vector<double> >zeroo(3,
		   vector<double>(3, 0.0));
	for(int ii = 0; ii < 3; ii++)
		for(int ij = 0; ij < 3; ij++){
			for(int j = 0; j < na; j++)
				zeroo[ii][ij] += zero[j][ii][ij];
			zeroo[ii][ij] /= na;
		}
	/*Now form phibar*/
	for(int i = 0; i < na; i++){
		for(int j = 0; j < na; j++)
			for(int ii = 0; ii < 3; ii++)
				for(int ij = 0; ij < 3; ij++){
					int correct = (zeroo[ii][ij] +  zeroo[ij][ii])*0.5 -
								  (zero[j][ii][ij] + zero[i][ij][ii]);
					for(int lx = 0; lx < nsc[0]; lx++)
					for(int ly = 0; ly < nsc[1]; ly++)
					for(int lz = 0; lz < nsc[2]; lz++)
						phi[i][ij][lx][ly][lz][j][ii] += correct;
				}
	}
	//loop over kpoints in transversal direction
	
	set<vector<double> > setkp;
	for(auto i : kp)
		setkp.insert({i[0], i[1], i[3]}); // kp[3] stores the weight of k.
	
	//copy 2d kpoints without duplication;

	vector<vector<double> > kpt;
	for_each(setkp.begin(), setkp.end(), 
			[&kpt](vector<double> e){
				kpt.push_back(e);});
	kp2D = kpt;


//	vector<vector<double> >dc(3*na,
//			   vector<double>(3*na, 0.0)));
	
	if(mfc.size() != kp2D.size())
		mfc.resize(kp2D.size());
	
	for(auto m1 : mfc)
		for(auto m2 : m1)
			fill(m2.begin(), m2.end(), 0.0);

	for(int k = 0; k < kp2D.size(); k++){
		for(int i = 0; i < na; i++)
		for(int j = 0; j < na; j++)
				for(int lx = 0; lx < nsc[0]; lx++)
				for(int ly = 0; ly < nsc[1]; ly++){
	//			for(int lz = 0; lz < nsc[2]; lz++){
					double r[3];
					for(int ii = 0; ii < 3; ii++)
						r[ii] = xa[i][ii] - xa[j][ii] +
							(lx-(nsc[0]-1)/2)*cell[0][ii] +
							(ly-(nsc[1]-1)/2)*cell[1][ii];// +
	//						 lz-(nsc[2]-1)/2)*cell[2][ii]);
					double dmin = 1.0e10;
					int neq = 1;
					vector<double> qr(1);
					for(int llx = -1; llx <=1; llx++)
					for(int lly = -1; lly <=1; lly++){
	//				for(int llz = -1; llz <=1; llz++){
						double R2 = 0.0;
						double r1[3]; 
						for(int ii = 0; ii < 3; ii++){
							r1[ii] = llx*scell[0][ii] + lly*scell[1][ii];// + 
	//							     llz*scell[2][ii];
							R2 += pow(r1[ii]+r[ii], 2.0);
						}
						// find equivalent cell.
						if(fabs(R2-dmin) > 1.0e-4){
							if(R2 < dmin){
								neq = 1;
								dmin = R2;
								qr[0] = kp2D[k][0] * (r1[0] + r[0]) +
							            kp2D[k][1] * (r1[1] + r[1]);// +
						//		        kp2D[k][2] * (r1[2] + r[2]);
							}
						}
						else if(fabs(R2-dmin) < 1.0e-4){
							neq++;
							double t = kp2D[k][0] * (r1[0] + r[0]) +
						               kp2D[k][1] * (r1[1] + r[1]);// +
						//	           kp2D[k][2] * (r1[2] + r[2]);
							qr.push_back(t);
						}
					}
					for(int q = 0; q < neq; q++){
						complex<double> phase = exp(complex<double>(0.0, 1.0)*qr[q]);
						for(int ii = 0; ii < 3; ii++)
						for(int jj = 0; jj < 3; jj++){
							int ix = i*3+ii; 
							int jx = j*3+jj;
							mfc[k][ix][jx] += phi[i][ii][lx][ly][(nsc[2]-1)/2][j][jj]*
											  phase/complex<double>(neq, 0.0)/
											  sqrt(xmass[i] * xmass[j]);
						}
					}
				}
	}
}
