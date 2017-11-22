#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include "FuncUtils.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::setw;

extern "C" void writensc(char* slabel, int nsc[], bool isbulk) ;
/*
 * Read slabel.XV after coordinates optimating calculations
 * and generate the atom coordinates in supercell
 *
 * int[] nsc : number of supercells in each direction
 *
 * The coor. info. is save in 'FC.fdf', one can use it in 
 * input.fdf by adding %include 'FC.fdf' at the end of file.
 * */
void writensc(char* slabel_, int nsc[], bool isbulk){
	vector<string> s0 = strsplit(string(slabel_));
	string slabel = s0[0];
	string coorSorceName = slabel + ".XV";
//	cout << "slabel.XV" << coorSorceName << endl;
//	cout << "nsc :" << nsc[0] <<" " <<nsc[1] <<" " << nsc[2] << endl;
	std::ifstream ifs(coorSorceName);
	string si;
	std::stringstream st;
	
	/*read cell*/	
	double cell[3][3];
	for(int i = 0; i < 3; i++){
		getline(ifs, si);
		st.clear();
		st.str(si);
		for(int j = 0; j < 3; j++)
			st >> cell[i][j];
	}
	/*read number of atom*/
	int na;
	getline(ifs, si);
	st.clear();
	st.str(si);
	st >> na;
	/*read coordinates*/
	double xyz[na][3];
	int nlabel[na];
	int xmass[na];

	for(int i = 0; i < na; i++){
		getline(ifs, si);
		st.clear();
		st.str(si);
		st >> nlabel[i] >> xmass[i];
		for(int j = 0; j < 3; j++)
			st >> xyz[i][j];
	}
	ifs.close();

/* Define number of unit cells*/
	int ncells = 1;
	for(int i = 0; i < 3; i++)
		ncells *= nsc[i];

/* Build lattice vector of the supercell*/
	double scell[3][3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			scell[i][j] = nsc[i]*cell[i][j];

/* Build atomic coordinates in the supercell*/
	int nsa = na*ncells;
	double xa[nsa][3];
	vector<vector<int> > cellIndex(nsa, vector<int>(3)); 
	// knows exactly what cells one atom belongs to.
	int nslabel[nsa];
	int Snumber[nsa];
	int iatom = -1;
	//loop over unit cells
	int lx = (nsc[0]-1)/2;
	int ly = (nsc[1]-1)/2;
	int lz = (nsc[2]-1)/2;
	for(int n1 = -lx; n1 <= lx; n1++)
		for(int n2 = -ly; n2 <= ly; n2++)
			for(int n3 = -lz; n3 <= lz; n3++){
				double r[3];
				for(int i = 0; i < 3; i++)
					r[i] = n1*cell[0][i] + n2*cell[1][i] + n3*cell[2][i];
				//loop over atoms in unit cell;
				for(int i = 0; i < na; i++){
					iatom++;
					nslabel[iatom] = nlabel[i];
					Snumber[iatom] = iatom+1;
					for(int j = 0; j < 3; j++){
						xa[iatom][j] = xyz[i][j] + r[j];
						cellIndex[iatom]={n1, n2, n3};
					}
				}
			}
	if(nsa != ++iatom) cout << "Error: computing atom numbers" << endl;

/* Determine the indices of the atoms in the centra cell(n1 = n2 = n3 =0)*/
	int i1 = na * (4*lx*ly*lz +
				   2*(lx*ly + lx*lz + ly*lz) +
				   lx + ly + lz) + 1;
	int i2 = i1 + na -1 ;

/* Read kpoints */
	string skp = slabel + ".KP";
	std::ifstream ifk(skp);
	getline(ifk, si);
	st.clear();
	st.str(si);
	int nkp;
	st >> nkp;
	double kp[nkp][4];
	for(int i = 0; i < nkp; i++){
		getline(ifk, si);
		st.clear();
		st.str(si);
		int indk;
		st >> indk;
		for(int j = 0; j < 4; j++)
			st >> kp[i][j];
	}
	ifk.close();
/*----------------------WRITE FDF FORMAT FILE-----------------------------*/
	// CellIndex info. will be used again in layered lead calc.
	string ssinfo = slabel;
	if(isbulk)
		ssinfo += ".BulkInfo";
	else
		ssinfo += ".EMInfo";
	std::fstream ofc(ssinfo, std::ios::out|std::ios::binary);
	ofc.write((char*)&na, sizeof(int));
//	ofc.write((char*)&xyz[0][0], sizeof(double)*(na*3));
	ofc.write(reinterpret_cast<char*>(xyz), sizeof(xyz));
	ofc.write((char*)xmass, sizeof(int)*na);
	ofc.write((char*)nsc, sizeof(int)*3);
//	ofc.write((char*)&cell[0][0], sizeof(double)*9);
//	ofc.write((char*)&scell[0][0], sizeof(double)*9);
	ofc.write(reinterpret_cast<char*>(cell), sizeof(cell));
	ofc.write(reinterpret_cast<char*>(scell), sizeof(scell));
	ofc.write((char*)&nkp, sizeof(int));
//	ofc.write((char*)&kp[0][0], sizeof(double)*(nkp*4));
	ofc.write(reinterpret_cast<char*>(kp), sizeof(kp));
	ofc.close();
	
	/*test binary in and out*/
//	std::fstream ifc(ssinfo, std::ios::in|std::ios::binary);
//	int nar;
//	ifc.read((char*)&nar, sizeof(int));
//	double xyzra[nar][3];
//	ifc.read(reinterpret_cast<char*>(xyzra), sizeof(xyzra));
//	
//	vector<vector<double> >xyzr;
//	mcopy<double>((double*)xyzra, nar, 3, xyzr);
//	
//	int xmassr[nar];
//	ifc.read((char*)xmassr, sizeof(int)*nar);
//	int nscr[3];
//	ifc.read((char*)nscr, sizeof(int)*3);
//	double cellr[3][3];
//	double scellr[3][3];
//	ifc.read(reinterpret_cast<char*>(cellr), sizeof(cellr));
//	ifc.read(reinterpret_cast<char*>(scellr), sizeof(scellr));
//	int nkpr;
//	ifc.read((char*)&nkpr, sizeof(int));
//	
//	double kpr[nkpr][4];
//	ifc.read(reinterpret_cast<char*>(kpr), sizeof(kpr));
//	
//	ifc.close();
////	ofc.write((char*)&kp[0][0], sizeof(double)*(nkp*4));
//	ofc.write(reinterpret_cast<char*>(kp), sizeof(kp));
//	cout << "Test binary write and read:\n";
//	cout << nar << endl;
//	for(int i = 0; i < nar; i++)
//		for(int j = 0; j < 3; j++)
//			cout << xyzr[i][j] << endl;
//	for(int i = 0; i < nar; i++)
//		cout << xmassr[i] << " ";
//	cout << endl;
//	for(int j = 0; j < 3; j++)
//		cout << nscr[j] << " ";
//	cout << endl;
//	for(int i = 0; i < 3; i++)
//		for(int j = 0; j < 3; j++){
//			cout << cellr[i][j] << " ";
//			cout << scellr[i][j] << " "; 
//		}
//	cout << endl;
//	cout << nkpr << endl;
//	for(int i = 0; i < nkpr; i++){
//		for(int j = 0; j < 4; j++)
//			cout << kpr[i][j] << " ";
//		cout << endl;
//	}

	std::ofstream ofs("FC.fdf");
	ofs  << "NumberOfAtoms		" << nsa << endl;
	
	ofs << "LatticeConstant";
	ofs << setiosflags(std::ios_base::fixed);
	ofs << std::setprecision(10) << setw(15) << 1.0 << "	Bohr" << endl;

	ofs << "%block LatticeVectors" << endl;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++)
			ofs << setw(15) << scell[i][j];
		ofs << endl;
	}
	ofs << "%endblock LatticeVectors" << endl;

	ofs << "AtomicCoordinatesFormat		NotScaledCartesianBohr" << endl;

	ofs << "%block AtomicCoordinatesAndAtomicSpecies" << endl;
	for(int i = 0; i < nsa; i++){
		for(int j = 0; j < 3; j++)
			ofs << setw(15) << xa[i][j];
		ofs << setw(6) << nslabel[i] << setw(6) << Snumber[i] << endl;
	}
	ofs << "%endblock AtomicCoordinatesAndAtomicSpecies" << endl;
	
	ofs << "%block kgrid_Monkhorst_Pack" << endl;
	ofs << "		1		0		0"<< endl;
	ofs << "		0		1		0"<< endl;
	ofs << "		0		0		1"<< endl;
	ofs << "%endblock kgrid_Monkhorst_Pack" << endl;

/* Write MD optins: Force Constant Calculations*/	
	ofs << "CalcPhononTRC		 T" << endl;
	if(isbulk)
		ofs << "IsBulkFC 		     T" << endl;
	else
		ofs << "IsBulkFC 		     F" << endl;

	ofs << "MD.TypeOfRun		FC" << endl;
	ofs << "MD.FCfirst			" << i1 << endl;
	ofs << "MD.FClast			" << i2 << endl;
	ofs << "MD.FCdispl			" << 0.04 << "	Bohr" << endl;
	ofs.close();

	cout << string(50,'*') << endl;
	cout << "		SuperCell info. has been written in FC.fdf." << endl;
	cout << string(50,'*') << endl;
}
