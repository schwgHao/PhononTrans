#include <iostream>
#include "ioFC.h"
#include <algorithm>
#include <complex>
#include "SurfacePhononGF.h"
#include <fstream>

using std::complex;
using std::fstream;

void CalcPhLead(string slabel, const vector<double>& omg, double delta){
//    if(Isbulk) ssinfo += ".BulkInfo";
//    else       ssinfo += ".MolInfo";
    iofc ifc(slabel, true);

    const size_t na = ifc.na;
    const size_t nscz = ifc.nsc.operator[](2);
    
    if(nscz > 3) {
        std::cerr << "cellz = " << nscz <<" , too small size along z."<< std::endl;
        abort();
    }
/* Read force constant matrix */

    vector<vector<vector<complex<double> > > > mfc(1,
                 vector<vector<complex<double> > >(3*na,
                          vector<complex<double> >(3*nscz*na, 0.0)));

//    iofc ifc(slabel, nsc, kp, xa, xmass, cell, scell);

    ifc.ReadFC(mfc);
    std::cout << "Read FC" << std::endl;
    for(size_t i = 0; i < mfc.size(); ++i){
        for(size_t lz = 0; lz < nscz; ++lz){
            for(size_t j1 = 0; j1 < na; ++j1){
                std::cout << "j1, j2 " << j1 << " ";
                for(size_t j2 = 0; j2 < na; ++j2){
                    std::cout << j2 << std::endl;
                    for(size_t k1 = 0; k1 < 3; ++k1){
                        for(size_t k2 = 0; k2 < 3; ++k2){
                            std::cout << mfc[i][3*j1+k1][3*(lz*na + j2) + k2] << "    ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
    }

    vector<vector<double> >kp2d = ifc.kp2d();

    vector<vector<int> > atominLayer = ifc.MatCutbyLayer();
//    std::cout << "Slice Supercell" << std::endl;    
    int nk = kp2d.size();
    const size_t naAtOnelayer = atominLayer[0].size();
//    int nlayers = atominLayer.size();
//    int naAtOnelayer = na/nlayers;
    vector<vector<vector<vector<vector<complex<double> > > > > >mfcbylz(nk,
           vector<vector<vector<vector<complex<double> > > > > (1,
                  vector<vector<vector<complex<double> > > > (2,
                         vector<vector<complex<double> > > (3*naAtOnelayer, 
                                vector<complex<double> >(3*naAtOnelayer, 0.0)))));

    ifc.mfcSlice(mfcbylz, mfc, atominLayer);

    std::cout << "mfcbylz 0, 0: " << std::endl;
    for(size_t i = 0; i < 3*naAtOnelayer; i++){
        for(size_t j = 0; j < 3; ++j){
            std::cout << mfcbylz[0][0][0][i][j] << "    ";
        }
        std::cout << std::endl;
    }

//    std::cout << "Slice Supercell" << std::endl;

    surfphGF sphGF(mfcbylz, omg, delta);
//    std::cout << "CalcPhLead: omg" << std::endl;
//    for(auto i : omg) std::cout << i << std::endl;  
    vector<vector<vector<vector<complex<double> > > > > d00r(nk,
           vector<vector<vector<complex<double> > > > (omg.size(),
                  vector<vector<complex<double> > > (3*naAtOnelayer,
                         vector<complex<double> > (3*naAtOnelayer, 0.0))));

    sphGF.DR00(d00r);
    std::cout << "Generate surface GF" << std::endl;
    for(size_t j = 0; j < omg.size(); ++j){
        for(auto i : d00r[0][j]){
            for(auto k : i)std::cout << k << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    /* write surface phonon green function to  */
    string ssgflft = "bulksgf.lft";    

    for(auto i : kp2d){
        for(auto j : i)
            std::cout << j << " ";
        std::cout << std::endl;
    }
    double cellLz = ifc.cell[2][2];
    sphGF.writeDR00(ssgflft, cellLz, naAtOnelayer, kp2d, d00r);
    std::cout << "Write DR00" << std::endl;
    string ssgfrt = "bulksgf.rt";
    /*if the two leads are same, just copy left lead to right*/
    if(true){
        fstream ifs(ssgflft, std::ios::in|std::ios::binary);
        fstream ofs(ssgfrt, std::ios::out|std::ios::binary);
        ofs << ifs.rdbuf();
        ifs.close();
        ofs.close();
    }
}
