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
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if ((i + j) % 2 == 0) v[size * i + j] = 1;
			else v[size * i + j] = -1;
		}
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
double originalEnergy(vector<double>& v, int size, vector < vector <double> >& na, double K)
{
	double e = 0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]] + K * v[na[i][1]] * v[na[i][3]] * v[na[i][7]]);
	}
	return e;
}
double nnnEne(vector<double>& v, int size, vector < vector <double> >& na, int ith)
{
	double nnn = 0;
	for (int i = 0; i < size * size; i++){
		nnn = nnn - v[i] * (v[na[i][4*ith-3]] + v[na[i][4*ith-1]]);
	}
	return nnn;
}
double delU1(vector<double>& v, int size, int i, vector < vector <double> >& na)
{
	double E = 2 * v[i] * (v[na[i][0]] + v[na[i][1]] + v[na[i][2]] + v[na[i][3]]);
	return E;
}
double delU2(vector<double>& v, int size, int i, vector < vector <double> >& na)
{
	double E = 2 * v[i] * ((v[na[i][4]]*v[na[i][2]] + v[na[i][5]]*v[na[i][3]])*v[na[i][0]] 
							+ (v[na[i][6]]*v[na[i][2]] + v[na[i][7]]*v[na[i][3]])*v[na[i][1]] );
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
void MC_1step(vector<double>& v, int size, double* expE, vector < vector <double> >& na, double K)
{
	double Ediff;
	double Ediff1;
	double Ediff2;
	gen.seed(rd);
	for (int i = 0; i < size; i=i+2) {
		for (int j = 0; j < size; j=j+2) {
			Ediff1 = delU1(v, size, size * i + j, na);
			Ediff2 = delU2(v, size, size * i + j, na);
			Ediff = Ediff1 + Ediff2 * K;
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= exp_delU(Ediff1, expE)*pow(exp_delU(Ediff2, expE), K)) v[size * i + j] = -v[size * i + j];
			}
		}
	}
	for (int i = 0; i < size; i=i+2) {
		for (int j = 1; j < size; j=j+2) {
			Ediff1 = delU1(v, size, size * i + j, na);
			Ediff2 = delU2(v, size, size * i + j, na);
			Ediff = Ediff1 + Ediff2 * K;
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= exp_delU(Ediff1, expE)*pow(exp_delU(Ediff2, expE), K)) v[size * i + j] = -v[size * i + j];
			}
		}
	}
	for (int i = 1; i < size; i=i+2) {
		for (int j = 0; j < size; j=j+2) {
			Ediff1 = delU1(v, size, size * i + j, na);
			Ediff2 = delU2(v, size, size * i + j, na);
			Ediff = Ediff1 + Ediff2 * K;
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= exp_delU(Ediff1, expE)*pow(exp_delU(Ediff2, expE), K)) v[size * i + j] = -v[size * i + j];
			}
		}
	}
	for (int i = 1; i < size; i=i+2) {
		for (int j = 1; j < size; j=j+2) {
			Ediff1 = delU1(v, size, size * i + j, na);
			Ediff2 = delU2(v, size, size * i + j, na);
			Ediff = Ediff1 + Ediff2 * K;
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= exp_delU(Ediff1, expE)*pow(exp_delU(Ediff2, expE), K)) v[size * i + j] = -v[size * i + j];
			}
		}
	}

}
void met_cycle(int size, double T, vector < vector <double> >& na, double K, 
vector<double>& magnet)
{
	int step1 = 2500, step2 = 100000;
	int trash_step =  size;
	if (T > 2.4 && T < 2.7) trash_step = trash_step * 2;

	double expE[4] = { exp(-8 / T), exp(-4 / T), exp(4 / T), exp(8 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na, K); }
	for (int k = 0; k < step2; k++) {
		MC_1step(array, size, &expE[0], na, K);
		magnet.at(k) = Magnet(array, size);
	}
}

int main()
{
	random_device rd;
	gen.seed(rd);
	double K = 0.2;
	int size = 10;
	double temp = 2.493;

	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	vector<double> magnet(100000, 0); 
	neighbor(near, size);

	ofstream Fileout;
	Fileout.open("10_p_at.txt");
	cout << "Fileout open: " << size << endl;
	Fileout << "size magnet " << endl;
	met_cycle(size, temp, near, K, magnet);
	for (int i = 0; i < 100000; i++){
		Fileout << size << " " << magnet.at(i) << endl;
	}
	Fileout.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
