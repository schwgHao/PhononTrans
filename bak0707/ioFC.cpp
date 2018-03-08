#include <fstream>
#include <sstream> 
#include <set>
#include <utility>
#include <algorithm>
#include <cmath>
#include <map>
#include "ioFC.h"
#include "FuncUtils.h"
#include <complex>

using std::fstream;
using std::stringstream;
using std::set;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;
using std::map;
using std::complex;

/* 
 * Read force from label.FC to construct force constant matrix.
 *
 * mfc[nk][3*na][3*na] has 3 dim.
 * nk is the number of 2d kpoints in transversal direction.
 *
 * The unit of elements is eV/ang^2, the conversion of sqrt(K/M) with eV and ang
 * to cm-1 is 519.6
 *
 * The private variable kp2D is initilized in ReadFC().
 * 2017/6/8, Hao Wang
 * */
iofc::iofc(string slabel_, bool isBulk_):slabel(slabel_), isBulk(isBulk_){
	string ssinfo = slabel;
	if(isBulk)
		ssinfo += ".BulkInfo";
	else 
		ssinfo += ".EMInfo";
	std::fstream fs(ssinfo, std::ios::in|std::ios::binary);
	if(!fs){
		std::cerr << "open "<< ssinfo <<" error!" << std::endl;
		abort();
	}
	
	
	fs.read((char*)&na, sizeof(int));
	
	if(!isBulk){
		fs.read((char*)&ia1, sizeof(int));
		fs.read((char*)&ia2, sizeof(int));
	}

	double xaa[na][3];
//	fs.read((char*)&xaa[0][0], sizeof(double)*(natoms*3));
	fs.read(reinterpret_cast<char*>(xaa), sizeof(xaa));
	mcopy<double>((double*)xaa, na, 3, xa);


	int xmassa[na];
	fs.read((char*)xmassa, sizeof(int)*na);
	xmass = vector<int>(na, 1);
	copy(xmassa, xmassa + na, xmass.begin());	
	
	
	int nsca[3];
	fs.read((char*)nsca, sizeof(int)*3);
	copy(nsca, nsca+3, back_inserter(nsc));
	double cella[3][3];
//	fs.read((char*)&cell[0][0], sizeof(double)*9);
	fs.read(reinterpret_cast<char*>(cella), sizeof(cella));
	mcopy<double>((double*)cella, 3, 3, cell);
	double scella[3][3];
//	fs.read((char*)&scell[0][0], sizeof(double)*9);
	fs.read(reinterpret_cast<char*>(scella), sizeof(scella));
	mcopy<double>((double*)scella, 3, 3, scell);
	

	
	int nkp;
	fs.read((char*)&nkp, sizeof(int));
	double kpa[nkp][4];
	fs.read(reinterpret_cast<char*>(kpa), sizeof(kpa));
	
	mcopy<double>((double*)kpa, nkp, 4, kp);
//	for(int i = 0; i < nkp; i++)
//		for(int j = 0; j < 4; j++)
//			cout << kp[i][j] << endl;
	
	fs.close();
}

void iofc::kp2dinit(){
	
//	set<vector<double> > setkp;
//	for_each(setkp.begin(), setkp.end(), 
//			[&kpt](vector<double> e){
//				kpt.push_back(e);});

	map<pair<double, double>, double> mapkp;
	for(auto i : kp)
		mapkp[make_pair(i[0], i[1])] = i[3];
//		setkp.insert({i[0], i[1], i[3]}); // kp[3] stores the weight of k.
	
	//copy 2d kpoints without duplication;

	vector<vector<double> > kpt;
	for(auto i = mapkp.begin(); i != mapkp.end(); i++)
		kpt.push_back({(i->first).first, (i->first).second, i->second});

	double totalw = 0.0;
	for(auto i : kpt)
		totalw += i[2];
	for(auto i = kpt.begin(); i != kpt.end(); i++)
		(*i)[2] /= totalw;
	
	kp2D = kpt;

}

void iofc::ReadFC(vector<vector<vector<complex<double> > > >& mfc){
//	cout << na << endl;
	string sn = slabel + ".FC";
	vector<vector<vector<vector<vector<double> > > > > pp(3, 
	 	   vector<vector<vector<vector<double> > > >(na, 
				  vector<vector<vector<double> > >(nsc[0], 
	                     vector<vector<double> >(nsc[1],
							    vector<double>(nsc[2], 0.0)))));
	auto pn = pp;
	
	vector<vector<vector<vector<vector<vector<vector<double> > > > > > > phi(3, 
		   vector<vector<vector<vector<vector<vector<double> > > > > > (na, 
				  vector<vector<vector<vector<vector<double> > > > >(3, 
						 vector<vector<vector<vector<double> > > >(na, 
								vector<vector<vector<double> > >(nsc[0], 
									   vector<vector<double> >(nsc[1],
											  vector<double>(nsc[2], 0.0)))))));

	auto phi0 = phi;
	
	fstream ifs(sn);
	string s1;
	getline(ifs, s1); // skip the title
	stringstream ist;
	for(int j = 0; j < na; j++){
		if(!isBulk && (j < ia1 || j > ia2)) continue; 
		for(int ij = 0; ij < 3; ij++){
			//read positive displacement
			for(int lx = 0; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++){
			for(int lz = 0; lz < nsc[2]; lz++){
				for(int i = 0; i < na; i++){
					getline(ifs, s1);
					ist.clear();
					ist.str(s1);
					for(int ii = 0; ii < 3; ii++){
						ist >> pp[ii][i][lx][ly][lz];
					}
				}
			}// lz
			}// ly	
			}// lx end read positive displ.
			// read negative displ.
			for(int lx = 0; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++){
			for(int lz = 0; lz < nsc[2]; lz++){
				for(int i = 0; i < na; i++){
					getline(ifs, s1);
					ist.clear();
					ist.str(s1);
					for(int ii = 0; ii < 3; ii++){
						ist >> pn[ii][i][lx][ly][lz];
					}
				}
			}// lz
			}// ly
			}// lx end read negative displ.
			// average displ.
			for(int lx = 0; lx < nsc[0]; lx++){
			for(int ly = 0; ly < nsc[1]; ly++){
			for(int lz = 0; lz < nsc[2]; lz++){
				for(int i = 0; i < na; i++){
					for(int ii = 0; ii < 3; ii++){
						phi0[ii][i][ij][j][lx][ly][lz] = 
							0.5*(pp[ii][i][lx][ly][lz] +
							     pn[ii][i][lx][ly][lz]);//average displ.
					}
				}
			}// lz
			}// ly
			}//	lx end average displ.
		} // ij
	}// j
	ifs.close();
	//make force constant matrix hermitian.
	for(int i = 0; i < na; i++){
		for(int j = 0; j < na; j++){
			for(int ii = 0; ii < 3; ii++){
				for(int jj = 0; jj < 3; jj++){
					for(int lx = 0; lx < nsc[0]; lx++){
					for(int ly = 0; ly < nsc[1]; ly++){
					for(int lz = 0; lz < nsc[2]; lz++){
						phi[ii][i][jj][j][lx][ly][lz] =
							0.5 * ( phi0[ii][i][jj][j][lx][ly][lz] +
							phi0[jj][j][ii][i][nsc[0]-lx-1][nsc[1]-ly-1][nsc[2]-lz-1]);
					}
					}
					}
				}
			}
		}
	}

	/*sum rule*/
	/* Now form zero(ix, jx, j)*/
//	int ncells = 1;
//	for(int i = 0; i < 3; i++)
//		ncells *= nsc[i];
//
//	vector<vector<vector<double> > > zero(3, 
//		   vector<vector<double> >(3,
//				  vector<double>(na, 0.0)));
//	for(int j = 0; j < na; j++){
//		for(int ii = 0; ii < 3; ii++)
//			for(int ij = 0; ij < 3; ij++){
//				for(int i = 0; i < na; i++){
//					for(int lx = 0; lx < nsc[0]; lx++)
//					for(int ly = 0; ly < nsc[1]; ly++)
//					for(int lz = 0; lz < nsc[2]; lz++)
//						zero[ii][ij][j] += phi[ii][i][ij][j][lx][ly][lz];
//				}
//				zero[ii][ij][j] /= (na*ncells);
//			}
//	}
//	/*Now form zeroo(ix, jx), the sum over j of zero(j, ix, jx)*/
//	vector<vector<double> >zeroo(3,
//		   vector<double>(3, 0.0));
//	for(int ii = 0; ii < 3; ii++)
//		for(int ij = 0; ij < 3; ij++){
//			for(int j = 0; j < na; j++)
//				zeroo[ii][ij] += zero[ii][ij][j];
//			zeroo[ii][ij] /= na;
//		}
//	/*Now form phibar*/
//	for(int i = 0; i < na; i++){
//		for(int j = 0; j < na; j++)
//			for(int ii = 0; ii < 3; ii++)
//				for(int ij = 0; ij < 3; ij++){
//					double correct = (zeroo[ii][ij] +  zeroo[ij][ii])*0.5 -
//								  (zero[ii][ij][j] + zero[ij][ii][i]);
//					for(int lx = 0; lx < nsc[0]; lx++)
//					for(int ly = 0; ly < nsc[1]; ly++)
//					for(int lz = 0; lz < nsc[2]; lz++)
//						phi[ii][i][ij][j][lx][ly][lz] += correct;
//				}
//	}
	/* end sum rule */
	/*solve 2d kpoints*/
	kp2dinit();	
	
//	for(int i = 0; i < kp2D.size(); i++)
//		cout << kp2D[i][2] << endl;

//	vector<vector<double> >dc(3*na,
//			   vector<double>(3*na, 0.0)));
	
	if(mfc.size() != kp2D.size())
		mfc.resize(kp2D.size());
	
	for(auto m1 : mfc)
		for(auto m2 : m1)
			fill(m2.begin(), m2.end(), 0.0);
	//loop over kpoints in transversal direction
	for(int k = 0; k < kp2D.size(); k++){
		for(int i = 0; i < na; i++){
			for(int j = 0; j < na; j++){
				for(int lx = 0; lx < nsc[0]; lx++){
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
					for(int llx = -1; llx <=1; llx++){
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
				//	}//end llz
					}//end lly
					}//end llx
					for(int q = 0; q < neq; q++){
						complex<double> phase = exp(complex<double>(0.0, 1.0)*qr[q]);
						for(int ii = 0; ii < 3; ii++)
						for(int jj = 0; jj < 3; jj++){
							int ix = i*3+ii; 
							int jx = j*3+jj;
							mfc[k][ix][jx] += phi[ii][i][jj][j][lx][ly][(nsc[2]-1)/2]*
											  phase/complex<double>(neq, 0.0);
											  
						}
					}
			//	}// end lz
				}//end ly
				}// end lx
			}// end j
		}// end i

		for(int ix = 0; ix < 3*na; ix++){
			for(int jx = 0; jx < 3*na; jx++){
				int i = ix/3.0;
				int j = jx/3.0;
				mfc[k][ix][jx] /= sqrt(1.0*xmass[i]*xmass[j]);
			}
		}
	}// end k
}

struct floatlessthan{
	bool operator()(const double& l, const double& r){
		return ((r - l) > 1.0e-4);
	}
};

/* Return index of atom by z position
 *
 * elements in atlayers, e.g.
 * 0.00000 : 0 1 2 3
 * 1.00000 : 4 5 6 7
 * ...
*/
vector<vector<int> > iofc::MatCutbyLayer(){
	map<double, vector<int>, floatlessthan> atlayers;
	for(int i = 0; i < na; i++)
		atlayers[xa[i][2]].push_back(i);
	vector<vector<int> > res;
	for(auto i = atlayers.begin(); i!= atlayers.end(); i++)
		res.push_back(i->second);
	return res;
}

/*reconstruct lead FC matrix by layer*/
void iofc::mfcSlice(vector<vector<vector<vector<vector<complex<double> > > > > >& mfcbylz, 
				    vector<vector<vector<complex<double> > > > mfc,
					vector<vector<int> > atominLayer){
//	int na = xa.size(); // number of atoms
	int nk = mfc.size(); // number of 2d kpoints
	int nl = mfcbylz[0].size(); // number of layers
//	int napl = mfcbylz[0][0][0].size()/3; //number of atoms per layer
	for(int k = 0; k < nk; k++){
		for(int lz1 = 0; lz1 < nl; lz1++){
			auto aind1 = atominLayer[lz1];
			for(int lz2 = 0; lz2 < nl; lz2++){
				auto aind2 = atominLayer[lz2];
				int c1 = -1;
				for(auto i : aind1){
					c1++;
					int c2 = -1;
					for(auto j : aind2){
						c2++;
						for(int ii = 0; ii < 3; ii++){
						for(int jj = 0; jj < 3; jj++)
							mfcbylz[k][lz1][lz2][c1+ii][c2+jj] = mfc[k][3*i+ii][3*j+jj];
						}
					}
				}
			}
		}
	}
}

/*Idenfity where's the EM area(KCInd) 
*and the interfaces among EM and 
*left lead(KLmInd) and right lead(KLpInd)
*/
void iofc::setBoundary(double cllft, double clrt, int nlrt,
		vector<int>& KLmInd, vector<int>& KLpInd, vector<int>& KCInd){
	/* ia1 = 1 and ia2 = na set the whole unit cell expanding molecule, 
	 * which must be cut into five areas:
	 * left lead, 
	 * left interface to mol., 
	 * mol., 
	 * right interface to mol., 
	 * right lead.
	 * */
	if((ia1 == 1) && (ia2 == na)){
		double Llftoutermost = cllft;
		double Lrtinnermost = cell[2][2] - clrt/nlrt*(nlrt+1);
		double minxaz = xa[0][2];
		for(int i = 1; i < xa.size(); i++)
			if(minxaz > xa[i][2]) minxaz = xa[i][2] ; 
		for(int i = 0; i < xa.size(); i++){
			double xar = xa[i][2] - minxaz;
			if((xar - Llftoutermost < -0.01) || (xar - Lrtinnermost) > 0.01) continue;
			else if(fabs(xar - Llftoutermost) < 1.0e-4)
				KLmInd.push_back(i);
			else if(fabs(xar - Lrtinnermost) < 1.0e-4)
				KLpInd.push_back(i);
			else
				KCInd.push_back(i);
		}
	}
	/*ia1 and ia2 indicate where the interfaces are*/
	else{
		for(int i = ia1; i <= ia2; i++){
			if(fabs(xa[i][2] - xa[ia1][2]) < 1.0e-4)
				KLmInd.push_back(i);
			else if(fabs(xa[i][2] - xa[ia2][2]) < 1.0e-4)
				KLpInd.push_back(i);
			else
				KCInd.push_back(i);
		}
	}
}

setleads::setleads(string ssl){
	fstream ifc(ssl, std::ios::in|std::ios::binary);
	double clz; // length of the unit cell along z;
	ifc.read((char*)&clz, sizeof(double));
	cellLz = clz;
	int nly; // number of layers
	ifc.read((char*)&nly, sizeof(int));
	nlayers = nly;
	
	int kpz1, kpz2;// size of k-points array n*2
	ifc.read((char*)&kpz1, sizeof(int));
	ifc.read((char*)&kpz2, sizeof(int));
	vector<vector<double> > kp2d(kpz1, vector<double>(kpz2));

	for(int k = 0; k < kpz1; k++)
		ifc.read((char*)&kp2d[k][0], kpz2*sizeof(double));

	kp = kp2d;
	
	int dz1, dz2, dz3, dz4;// size of D00R nk*nomga*(3*natomperlayer)*(3*natomperlayer)
	ifc.read((char*)&dz1, sizeof(int));
	ifc.read((char*)&dz2, sizeof(int));
	ifc.read((char*)&dz3, sizeof(int));
	ifc.read((char*)&dz4, sizeof(int));
	
	vector<vector<vector<vector<complex<double> > > > > d00r_(dz1,
		   vector<vector<vector<complex<double> > > > (dz2,
		          vector<vector<complex<double> > > (dz3,
				         vector<complex<double> > (dz4, 0.0))));
	
	for(int i = 0; i < dz1; i++)
		for(int j = 0; j < dz2; j++)
			for(int k = 0; k < dz3; k++)
				ifc.read((char*)(&d00r_[i][j][k][0]), dz4*sizeof(complex<double>));

	d00r = d00r_;
}
