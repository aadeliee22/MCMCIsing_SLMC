#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>
using namespace std;

random_device rd;
trng::yarn2 gen;
trng::uniform01_dist<> dis;

void initialize(vector<double>& v, int size) //initial -random- state
{
	for (int i = 0; i < size * size; i++) {
		v[i] = dis(gen) < 0.5 ? 1 : -1;
	}
}
void neighbor(vector < vector <double> >& na, int size)
{
	int sizes = size * size;
	for (int i = 0; i < size * size; i++) {
		int na_2 = (i - 1 + size) % size + (i / size) * size;
		int na_3 = (i + 1) % size + (i / size) * size;
		na[i][0] = (i + size * (size - 1)) % sizes;
		na[i][1] = (i + size) % sizes;
		na[i][2] = na_2;
		na[i][3] = na_3;
		na[i][4] = (na_2 + size * (size - 1)) % sizes;
		na[i][5] = (na_3 + size * (size - 1)) % sizes;
		na[i][6] = (na_2 + size) % sizes;
		na[i][7] = (na_3 + size) % sizes;
		na[i][8] = (i + size * (size - 2)) % sizes;
		na[i][9] = (i + 2 * size) % sizes;
		na[i][10] = (i - 2 + size) % size + (i / size) * size;
		na[i][11] = (i + 2) % size + (i / size) * size;
	}
}
double Magnet(vector<double>& v, int size)
{
	double m = 0;
	for (int i = 0; i < size * size; i++){
		m = m + v.at(i);
	}
	m = m / (size*size);
	return abs(m);
}
double nnnEne(vector<double> v, int size, vector < vector <double> > na, int ith)
{
	double nnn = 0;
	for (int i = 0; i < size * size; i++){
		nnn = nnn - v[i] * (v[na[i][4*ith-3]] + v[na[i][4*ith-1]]);
	}
	return nnn;
}
double Energy(vector<double> v, int size, vector < vector <double> > na, vector<double> J, int nth)
{
	double e = J[0];
	for (int i = 1; i < nth + 1; i++){
		e = e + J[i] * nnnEne(v, size, na, i);
	}
	return e;
}
double delU(vector<double> v, int size, int i, vector < vector <double> > na, int ith)
{
	double E = 2 * v[i] * (v[na[i][4*ith-4]] + v[na[i][4*ith-3]] + v[na[i][4*ith-2]] + v[na[i][4*ith-1]]);
	return E;
}
double exp_delU(double E, double* expE)
{
	double result1;
	if (E == 8) result1 = *(expE);
	else if (E == 4) result1 = *(expE + 1);
	else if (E == -4) result1 = *(expE + 2);
	else if (E == -8) result1 = *(expE + 3);
	else result1 = 1;
	return result1;
}
void MC_1step(vector<double>& v, int size, double* expE, vector < vector <double> > na, vector<double> J)
{
	double Ediff, Ediff1, Ediff2, Ediff3;
	gen.seed(rd);
	for (int i = 0; i < size; i=i+2) {
		for (int j = 0; j < size; j=j+2) {
			Ediff1 = delU(v, size, size * i + j, na, 1);
			Ediff2 = delU(v, size, size * i + j, na, 2);
			Ediff3 = delU(v, size, size * i + j, na, 3);
			Ediff = Ediff1 * J[1] + Ediff2 * J[2] + Ediff3 * J[3];
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= pow(exp_delU(Ediff1, expE), J[1])
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) 
				{v[size * i + j] = -v[size * i + j];}
			}
		}
	}
	for (int i = 0; i < size; i=i+2) {
		for (int j = 1; j < size; j=j+2) {
			Ediff1 = delU(v, size, size * i + j, na, 1);
			Ediff2 = delU(v, size, size * i + j, na, 2);
			Ediff3 = delU(v, size, size * i + j, na, 3);
			Ediff = Ediff1 * J[1] + Ediff2 * J[2] + Ediff3 * J[3];
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= pow(exp_delU(Ediff1, expE), J[1])
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) 
				{v[size * i + j] = -v[size * i + j];}
			}
		}
	}
	for (int i = 1; i < size; i=i+2) {
		for (int j = 0; j < size; j=j+2) {
			Ediff1 = delU(v, size, size * i + j, na, 1);
			Ediff2 = delU(v, size, size * i + j, na, 2);
			Ediff3 = delU(v, size, size * i + j, na, 3);
			Ediff = Ediff1 * J[1] + Ediff2 * J[2] + Ediff3 * J[3];
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= pow(exp_delU(Ediff1, expE), J[1]) 
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) 
				{v[size * i + j] = -v[size * i + j];}
			}
		}
	}
	for (int i = 1; i < size; i=i+2) {
		for (int j = 1; j < size; j=j+2) {
			Ediff1 = delU(v, size, size * i + j, na, 1);
			Ediff2 = delU(v, size, size * i + j, na, 2);
			Ediff3 = delU(v, size, size * i + j, na, 3);
			Ediff = Ediff1 * J[1] + Ediff2 * J[2] + Ediff3 * J[3];
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= pow(exp_delU(Ediff1, expE), J[1])
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3]))
				{v[size * i + j] = -v[size * i + j];}
			}
		}
	}
}
void met_cycle(int size, double T, vector < vector <double> > na, vector<double> J, double& ene, int nth, int step2)
{
	int step1 = 2500;
	int trash_step =  size*size/50;
	if (T > 2.4 && T < 2.7) trash_step = trash_step * 2;

	double expE[4] = { exp(-8 / T), exp(-4 / T), exp(4 / T), exp(8 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na, J); } 
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			MC_1step(array, size, &expE[0], na, J);
		}
		//magnet.at(k) = Magnet(array, size);
		ene = ene + Energy(array, size, na, J, nth);
	}
	ene = ene / step2;
}
int main()
{
	random_device rd; gen.seed(rd);
	int size, nth, step2;
	vector<double> J(4, 0);
	double K = 0.2; double Tc = 2.493; double temp;
	double ene = 0;

	//filein.txt format: nth \n temperature \n E0 \n J1 \n J2 \n J3
	ifstream Filein; Filein.open("filein_met.txt"); 
	Filein >> size;
	Filein >> nth; 
	Filein >> step2; 
	for (int i = 0; i < nth + 1; i++){ Filein >> J[i]; }

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	
	neighbor(near, size);

	clock_t start = clock();	

	ofstream Fileout;
	Fileout.open("plot_met_a.txt");
	cout << "File open: " << endl;
	Fileout << "J: " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << " " << endl;
	Fileout << "sizes nth step2 temp ene " << endl;
	for (int k = 300; k < 500; k++) {
		temp = 0.005 * k;
		for (int h = 0; h < 25; h++) {
			met_cycle(size, temp, near, J, ene, nth, step2);
			Fileout << size << " " << nth << " " << step2 << " " << temp << " " << ene << " " << endl;
		}
	}
	Fileout.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
