#include "EMsolver.h"
#include <fstream>
#include <iomanip>
#include <Dense>
#include <algorithm>
#include "nr3.h"
#include "interp_1d.h"
#include "qgaus.h"
#include "romberg.h"
#include <iterator>
using namespace Eigen;

void emsolver::vibProj(const vector<vector<vector<vector<complex<double> > > > >& d00rl, const vector<vector<vector<vector<complex<double> > > > >& d00rr)
{
    const double unitconv = 98226947884640.62;
    const double hbar = 6.582119144686551e-16; // in units of eV*s
    
    int kcis = KCInd.size();
    int klmis = KLmInd.size();
    int klpis = KLpInd.size();
    string frvib = slabel + ".EigenModeRate";
    std::fstream orfs(frvib, std::ios::out);
    string feMode = slabel + ".EigenMode";
    std::fstream ivfs(feMode, std::ios::in);
    size_t nkp, nac;
    ivfs >> nkp >> nac;
    orfs << nkp << "    "<<nac << std::endl;
    for(size_t k = 0; k < nkp; k++){
        for(int i = 0; i < 3; ++i) orfs << kp[k][i] << "    ";
        orfs << std::endl;
        Eigen::MatrixXcd KC(3*kcis, 3*kcis);
        Eigen::MatrixXcd KLm(3*klmis, 3*kcis);
        Eigen::MatrixXcd KLp(3*kcis, 3*klpis);
        for(int i = 0; i < kcis; i++){
            for(int j = 0; j < kcis; j++)
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++){
                KC(3*i+ii, 3*j+jj) = mfc[k][3*KCInd[i]+ii][3*KCInd[j]+jj];
            }
            for(int j = 0; j < klmis; j++){
              //  std::cout << "i, j " << i << " " << j << std::endl;
                for(int ii = 0; ii < 3; ii++){
                    for(int jj = 0; jj < 3; jj++){
                        KLm(3*j+jj, 3*i+ii) = mfc[k][3*KLmInd[j]+jj][3*KCInd[i]+ii];
                        KLp(3*i+ii, 3*j+jj) = mfc[k][3*KLpInd[j]+jj][3*KCInd[i]+ii];
                //        std::cout << KLm(3*j+jj, 3*i+ii) << " ";
                    }
                //    std::cout << std::endl;
                }
            //    std::cout << std::endl;
            }
        }
        double tmp;
        for(int i = 0; i < 3; ++i) ivfs >> tmp;
        for(int i = 0; i < nac; ++i) ivfs >> tmp;
        MatrixXcd vmod(1, nac);
        for(int o = 0; o < nac; ++o){
            for(int i = 0; i < nac; ++i){
                ivfs >> vmod(i);
            }
            Eigen::MatrixXcd Iden = Eigen::MatrixXcd::Identity(3*kcis, 3*kcis);
            complex<double> dtli = complex<double>(omg[o], delta); 
            Iden *= dtli*dtli;
            Eigen::MatrixXcd D00l(3*klmis, 3*klmis);
            Eigen::MatrixXcd D00r(3*klpis, 3*klpis);
            for(int i = 0; i < 3*klmis; i++)
                for(int j = 0; j < 3*klmis; j++){
                    D00l(i, j) = d00rl[k][o][i][j];
                    D00r(i, j) = d00rr[k][o][i][j];
                }
            Eigen::MatrixXcd SelfEngL = KLm.transpose()*D00l*KLm;
            Eigen::MatrixXcd SelfEngR = KLp*D00r*KLp.transpose();
            orfs << std::scientific << hbar * unitconv * 1000.0 * omg[o] << "        ";//meV
            orfs << -1.0/omg[o]/omg[o] * vmod.array().real().matrix() * SelfEngL.array().imag().matrix() * vmod.transpose().array().real().matrix() << "        ";
            orfs << -1.0/omg[o]/omg[o] * vmod.array().real().matrix() * SelfEngR.array().imag().matrix() * vmod.transpose().array().real().matrix() << std::endl;
        }
    }
    orfs.close();
    ivfs.close();
}


void emsolver::DRCC(vector<vector<vector<vector<complex<double> > > > >& dccr, 
                    vector<vector<complex<double> > >& VibTRCOverk,
                    vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadL, 
                    vector<vector<vector<vector<complex<double> > > > >& SelfEngLeadR,
                    const vector<vector<vector<vector<complex<double> > > > >& d00rl, 
                    const vector<vector<vector<vector<complex<double> > > > >& d00rr){
    const double unitconv = 98226947884640.62;
    const double hbar = 6.582119144686551e-16; // in units of eV*s
    
    std::cout << string(60,'*')<<std::endl;
    std::cout << string(60,'*')<<std::endl;
    std::cout << "siesta: Solve retarded phonon green's function." << std::endl;
    int kcis = KCInd.size();
    int klmis = KLmInd.size();
    int klpis = KLpInd.size();
//    std::cout << "kcis klmis klpis: " << kcis << " " << klmis <<" " << klpis << std::endl;  
    std::cout << "EMsolver: KCInd" << std::endl;
    for(auto i : KCInd) std::cout << i << " ";
    std::cout << std::endl;
    std::cout << "EMsolver: KLmInd" << std::endl;
    for(auto i : KLmInd) std::cout << i << " ";
    std::cout << std::endl;
    std::cout << "EMsolver: KLpInd" << std::endl;
    for(auto i : KLpInd) std::cout << i << " ";
    std::cout << std::endl;
    
    // write DOS, eigenmode and eigenchannel
    string fbDos = slabel + ".BulkDos";
    std::fstream ofs(fbDos, std::ios::out);
    ofs << "# Energy(meV) " << "     bulk.DOS (1 / meV) " << " EM.DOS (1 / meV)        " <<std::endl;
    
    string feMode = slabel + ".EigenMode";
    std::fstream ovfs(feMode, std::ios::out);
    ovfs << "# kp.size()  EM.size()" <<std::endl;
    ovfs << "# kp" <<std::endl;
    ovfs << "# eigenvalues: 1 * " << 3*kcis << std::endl;
    ovfs << "# eigenmodes: " << 3*kcis << " * " << 3*kcis << std::endl;
    
    string fchl = slabel + ".EigenChannel";
    std::fstream ocfs(fchl, std::ios::out);
    ocfs << "# kp.size()  energy.size()  EM.size()" <<std::endl;
    ocfs << "# kp" <<std::endl;
    ocfs << "# Energy(meV)      eigenchannel1 ... eigenchannel" << 3*kcis << std::endl;
   
    ofs  << kp.size() << "     " << omg.size();
    ovfs << kp.size() << "     " << 3*kcis << std::endl;
    ocfs << kp.size() << "     " << omg.size() << "   "  << 3*kcis << std::endl;
    for(size_t k = 0; k < kp.size(); k++){
        Eigen::MatrixXcd KC(3*kcis, 3*kcis);
        Eigen::MatrixXcd KLm(3*klmis, 3*kcis);
        Eigen::MatrixXcd KLp(3*kcis, 3*klpis);
        for(int i = 0; i < kcis; i++){
            for(int j = 0; j < kcis; j++)
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++){
                KC(3*i+ii, 3*j+jj) = mfc[k][3*KCInd[i]+ii][3*KCInd[j]+jj];
            }
            for(int j = 0; j < klmis; j++){
//                std::cout << "i, j " << i << " " << j << std::endl;
                for(int ii = 0; ii < 3; ii++){
                    for(int jj = 0; jj < 3; jj++){
                        KLm(3*j+jj, 3*i+ii) = mfc[k][3*KLmInd[j]+jj][3*KCInd[i]+ii];
                        KLp(3*i+ii, 3*j+jj) = mfc[k][3*KLpInd[j]+jj][3*KCInd[i]+ii];
//                        std::cout << KLm(3*j+jj, 3*i+ii) << " ";
                    }
//                    std::cout << std::endl;
                }
//                std::cout << std::endl;
            }
        }
        //for(auto i : kp[k]) ovfs << i << "  ";
        //ovfs << std::endl;
        std::copy(kp[k].begin(), kp[k].end(), std::ostream_iterator<double>(ovfs, " "));
        ovfs << std::endl;
        std::copy(kp[k].begin(), kp[k].end(), std::ostream_iterator<double>(ocfs, " "));
        ocfs << std::endl;
        std::copy(kp[k].begin(), kp[k].end(), std::ostream_iterator<double>(ofs, " "));
        ocfs << std::endl;
        SelfAdjointEigenSolver<MatrixXcd> eigensolver(KC);
        ovfs << eigensolver.eigenvalues().transpose().array().abs().sqrt().matrix() << std::endl;
        ovfs << eigensolver.eigenvectors().transpose() << std::endl;
        //        std::cout << "KLm: " << KLm << std::endl;
        //std::cout << "KLp: " << KLp << std::endl;
        //std::cout << "Kc: " << KC << std::endl;
        for(size_t o = 0; o < omg.size(); o++){
            Eigen::MatrixXcd Iden = Eigen::MatrixXcd::Identity(3*kcis, 3*kcis);
            complex<double> dtli = complex<double>(omg[o], delta); 
            Iden *= dtli*dtli;
            Eigen::MatrixXcd D00l(3*klmis, 3*klmis);
            Eigen::MatrixXcd D00r(3*klpis, 3*klpis);
            for(int i = 0; i < 3*klmis; i++)
                for(int j = 0; j < 3*klmis; j++){
                    D00l(i, j) = d00rl[k][o][i][j];
                    D00r(i, j) = d00rr[k][o][i][j];
                }
            ofs << std::scientific << hbar * unitconv * 1000.0 * omg[o] << "        " << -0.636619772367581 / (unitconv * hbar * 1000.0)  *  omg[o] * D00l.trace().imag() << "        "; // -2/pi*omg*Im[Dr(omg^2)]
            Eigen::MatrixXcd SelfEngL = KLm.transpose()*D00l*KLm;
            Eigen::MatrixXcd SelfEngR = KLp*D00r*KLp.transpose();
        //    std::cout << "Iden: " << std::scientific << Iden << std::endl;
        //    std::cout << "KC: " << std::scientific << KC << std::endl;
//            ofs << std::scientific << SelfEngL.real() << std::endl;
//            ofs << std::scientific << SelfEngL.imag() << std::endl;
        //    std::cout << "SL: " << std::scientific << SelfEngL << std::endl;
        //    std::cout << "SR: " << std::scientific << SelfEngR << std::endl;
//            Eigen::MatrixXcd TC = (Iden - KC - SelfEngL - SelfEngR).fullPivLu().inverse();
            Eigen::MatrixXcd TC = (Iden - KC - SelfEngL - SelfEngR).householderQr().solve(Eigen::MatrixXcd::Identity(3 * kcis, 3 * kcis));
            ofs << std::scientific <<  -0.636619772367581 / (unitconv * hbar * 1000.0)  *  omg[o] * TC.trace().imag() << std::endl;
            for(int i = 0; i < 3*kcis; i++)
                for(int j = 0; j < 3*kcis; j++){
                    dccr[k][o][i][j] = TC(i, j);
                    SelfEngLeadL[k][o][i][j] = SelfEngL(i, j);
                    SelfEngLeadR[k][o][i][j] = SelfEngR(i, j);
                }
        //    Eigen::MatrixXcd vibT =(-1.0)*(SelfEngL-SelfEngL.adjoint())*TC.adjoint()*
    //                               (SelfEngR-SelfEngR.adjoint())*TC;
            Eigen::MatrixXcd vibT =(-1.0)*TC*(SelfEngL-SelfEngL.adjoint())*TC.adjoint()*(SelfEngR-SelfEngR.adjoint());
            Eigen::SelfAdjointEigenSolver<MatrixXcd> eigenchl(vibT);
            ocfs << std::scientific <<  unitconv * hbar * 1000.0 * omg[o] << "  "; //meV
            ocfs << eigenchl.eigenvalues().transpose() << std::endl;
            
            VibTRCOverk[k][o] = vibT.trace();
    //        if(omg.size()%(o+1) == 0) 
//                std::cout << std::fixed<<std::setprecision(0)<<o*1.0/(omg.size()-1)*100 << " % completed." << std::endl;
        }
    }
    ofs.close();
    ovfs.close();
    ocfs.close();
}



struct SpectraConductance{
    static constexpr double Kb = 1.380650524e-23; //  J/K
    static constexpr double pi = 3.141592653589793;
    static constexpr double hbar = 6.6260695729e-34/2.0/pi; // in units of J*s
    double T;
    VecDoub omga;
    VecDoub Tomga;
    SpectraConductance(){}
    SpectraConductance(double T_, VecDoub omga_, 
        VecDoub Tomga_): T(T_), omga(omga_), Tomga(Tomga_){}
    double operator()(double omgi){

        Poly_interp pinterp = Poly_interp(omga, Tomga, 4);

        auto dfBEdT = [this](double x){return 1.0/pow((exp(hbar*x/Kb/T)-1.0), 2.0)*exp(hbar*x/Kb/T)*(hbar*x/Kb/T/T);};
//        std::cout << "ff: " << pinterp.interp(omgi) << std::endl;
        return pinterp.interp(omgi)*hbar*omgi*dfBEdT(omgi)/2.0/pi;
    }
};

double funcInt(double omgi, double T){
    const double Kb = 1.380650524e-23; //  J/K
    const double pi = 3.141592653589793;
    const double hbar = 6.6260695729e-34/2.0/pi; // in units of eV*s
    return 1.0/pow(exp(hbar*omgi/Kb/T)-1.0, 2.0)
           *exp(hbar*omgi/Kb/T)*hbar*omgi/Kb/T/T*hbar*omgi/2.0/pi;
    
}

double quadsimp(vector<double> y, double h){
    int n = y.size();
    vector<double> co(n);
    double sum = 0.0;
    
    if(n == 1) return 0.0;
    else if(n == 2){return h*(y[0]+y[1])*0.5;}
    else if(n == 3){return h*(1.0/3.0*(y[0]+y[2])+4.0/3.0*y[1]);}
    else{
        co[0] = 1.0/3.0;
        co[n-1] = 1.0/3.0;
        for(int i = 1; i < (int)(n-2)/2+1; i++){
            co[2*i-1] = 4.0/3.0;
            co[2*i] = 2.0/3.0;
        }
        if(n & 1) co[n-2] = 4.0/3.0;
        for(int i = 0; i < n; i++)
            sum += co[i]* y[i]*h;
    }
    return sum;
}

void emsolver::writePhTRC(const string& slabel, const vector<vector<complex<double> > >& vibTRCoverk){
//    double unitconv = sqrt(e/Ang2m/Ang2m/amu2kg); // sqrt(eV/(Ang^2*amu)) -> s-1
//    math.sqrt(1.602e-19/1.66e-27/1.0e-10/1.0e-10) = 98226947884640.62
    const double unitconv = 98226947884640.62;
    const double hbar = 6.582119144686551e-16; // in units of eV*s
    int No = omg.size();
    VecDoub Omgv(No, 0.0);
    VecDoub Tomg(No, 0.0);    

    string sTrc = slabel + ".PhTrc";
    std::fstream ifsTrc(sTrc, std::ios::out);
    ifsTrc << setiosflags(std::ios_base::fixed);
    ifsTrc << "    Energy(meV)    " << "    " << "Phonon Transmission Coeff." << std::endl;  
    
    for(int o = 0; o < No; o++){
        Omgv[o] = omg[o]*unitconv; // now elements in Omgv in units of s-1
        for(size_t k = 0; k < kp.size(); k++)
            Tomg[o] += vibTRCoverk[k][o].real();
        ifsTrc << Omgv[o] * hbar * 1000.0 << "            " << std::scientific << Tomg[o] << std::endl;  
    }

    ifsTrc.close();

    std::cout << string(45,'*')<<std::endl;
    std::cout << "    Write Thermal conductance." << std::endl;
    std::cout << string(45,'*')<<std::endl;
    
    string strs = slabel + ".THERMPH";
    
    std::fstream ifs(strs, std::ios::out);
    ifs << setiosflags(std::ios_base::fixed);
    ifs << "    Temperature(K)    " << "    " << "Thermal Condunctance(W/K)" << std::endl;  
/*    try{    
    // interpolate and integrate instead
        SpectraConductance AGph; 
        for(int tp = 0; tp < bathT.size(); tp++){
            ifs << std::setprecision(10) << std::setw(15)  << bathT[tp] << "        ";
            AGph = SpectraConductance(bathT[tp], Omgv, Tomg);
            double Gph = qtrap(AGph, Omgv[0], Omgv[No-1], 1.0e-5);
            ifs << std::setw(15)<< std::scientific << Gph << std::endl;
            if(bathT.size()%(tp+1) == 0){
                std::cout << std::fixed<<std::setprecision(0);
                std::cout <<(tp+1)*1.0/(bathT.size())*100;
                std::cout << " % completed." << std::endl;
            }
        }
    }
    catch(char*){ 
*/
    /*integrate thermal conductance directly*/
        double h = Omgv[1] - Omgv[0];
        for(size_t tp = 0; tp < bathT.size(); tp++){
            ifs << std::setprecision(10) << std::setw(15)  << bathT[tp] << "        ";
            vector<double> fInt(No);
            for(int o = 0; o < No; o++){
                fInt[o] = abs(Tomg[o])*funcInt(Omgv[o], bathT[tp]);
            }
            double Gph = quadsimp(fInt, h);
            ifs << std::setw(15)<< std::scientific << Gph << std::endl;
            if(bathT.size()%(tp+1) == 0){
                std::cout << std::fixed<<std::setprecision(0);
                std::cout <<(tp+1)*1.0/(bathT.size())*100;
                std::cout << " % completed." << std::endl;
            }
        }
/* }
*/
    ifs.close();
}
