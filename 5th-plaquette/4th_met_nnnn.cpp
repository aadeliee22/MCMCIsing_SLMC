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
	for (vector<int>::size_type i=0; i<v.size(); i++){
		m = m + v.at(i);
	}
	m = m / (v.size());
	return abs(m);
}
double nnnEne(vector<double>& v, int size, vector < vector <double> >& na, int ith)
{
	double nnn = 0;
	for (int i = 0; i < size * size; i++){
		nnn = nnn - v[i] * (v[na[i][4*ith-3]] + v[na[i][4*ith-1]]);
	}
	return nnn;
}
double Energy(vector<double>& v, int size, vector < vector <double> >& na, vector<double>& J, int nth)
{
	double e = J[0];
	for (int i = 1; i < nth + 1; i++){
		e = e + J[i] * nnnEne(v, size, na, i);
	}
	return e;
}
double delU(vector<double>& v, int size, int i, vector < vector <double> >& na, int ith)
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
void MC_1step(vector<double>& v, int size, double* expE, vector < vector <double> >& na, vector<double>& J)
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
				if (dis(gen) <= exp_delU(Ediff1, expE) 
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) v[size * i + j] = -v[size * i + j];
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
				if (dis(gen) <= exp_delU(Ediff1, expE) 
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) v[size * i + j] = -v[size * i + j];
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
				if (dis(gen) <= exp_delU(Ediff1, expE) 
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) v[size * i + j] = -v[size * i + j];
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
				if (dis(gen) <= exp_delU(Ediff1, expE) 
				* pow(exp_delU(Ediff2, expE), J[2]) 
				* pow(exp_delU(Ediff3, expE), J[3])) v[size * i + j] = -v[size * i + j];
			}
		}
	}
}
void MC_1cycle(int size, double T, double& mag, double& mag_sus, double& mag2, double& mag4, double& ene, double& sp_heat, vector < vector <double> >& na, vector<double>& J, int nth)
{
	int step1 = 2500, step2 = 10000;
	int trash_step =  size+5;
	if (T > 2.3 && T < 2.7) trash_step = trash_step * 2;

	double expE[4] = { exp(-8 / T), exp(-4 / T), exp(4 / T), exp(8 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na, J); }

	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			MC_1step(array, size, &expE[0], na, J);
		}
		magnet.at(k) = Magnet(array, size);
		energy.at(k) = Energy(array, size, na, J, nth);
	}

	double Mag = 0, Mag2 = 0, Mag4 = 0, Ene = 0, Ene2 = 0;
	double a, b;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		a = magnet.at(i);
		b = energy.at(i);
		Mag = Mag + a;
		Mag2 = Mag2 + a*a;
		Mag4 = Mag4 + a*a*a*a;
		Ene = Ene + b;
		Ene2 = Ene2 + b*b;
	}

	mag = Mag / step2;
	ene = Ene / step2;
	mag_sus = size*size * (Mag2 / step2 - (Mag/step2)*(Mag/step2)) / T;
	mag2 = Mag2 / step2;
	mag4 = Mag4 / step2;
	sp_heat = (Ene2 / step2 - (Ene/step2)*(Ene/step2)) / (size*size*T*T);
}

int main()
{
	random_device rd; gen.seed(rd);
	int size = 10; int nth = 3;
	vector<double> J(nth + 1, 0);
	J[0] = 0; J[1] = 1;

	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0, Ene = 0, sp_heat = 0;
	double temp = 0;

	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(4*nth, 0));
	neighbor(near, size);

	ofstream File;
	File.open("met_3.txt");
	cout << "File met3 open: " << endl;
	File << "J2 J3 temperature m m2 m4 ms e sp_h " << endl;
	for (int m = 0; m < 5; m++){
		J[2] = 0.1 - 0.05 * m;
		for (int n = 0; n < 7; n++){
			J[3] = 0.075 - 0.025 * n;
			for (int k = 200; k < 700; k++) {
				temp = 0.005 * k;
				for (int h = 0; h < 2; h++) {
					MC_1cycle(size, temp, Mag, mag_sus, Mag2, Mag4, Ene, sp_heat, near, J, nth);
					File << J[2] << " " << J[3] << " " << temp << " " << Mag << " " << Mag2 << " " << Mag4 
					<< " " << mag_sus << " " << Ene << " " << sp_heat << " " << endl;
				}
			}
			cout << J[3] << " end * " << endl;
		}
		cout << J[2] << " end " << endl;
	}
	File.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
