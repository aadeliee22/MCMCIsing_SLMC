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
	/*for (int i = 0; i < size * size; i++) {
		v[i] = dis(gen) < 0.5 ? 1 : -1;
	}*/
}
void color(vector<double>& v, int size) //graphing state
{
	for (int i = 0; i < size * size; i++) {
		if (i % size == 0) cout << endl;
		if (v[i] == 1) cout << "* ";
		if (v[i] == -1) cout << ". ";
	}
}
void neighbor(vector < vector <double> >& na, int size)
{
	int sizes = size * size;
	for (int i = 0; i < size * size; i++) {
		na[i][0] = (i + size * (size - 1)) % sizes;
		na[i][1] = (i + size) % sizes;
		na[i][2] = (i - 1 + size) % size + (i / size) * size;
		na[i][3] = (i + 1) % size + (i / size) * size;
	}
}
double Magnet(vector<double>& v, int size)
{
	double m = 0;
	for (vector<int>::size_type i = 0; i < v.size(); i++) {
		m = m + v.at(i);
	}
	m = abs(m) / (v.size()); //absolute value of average spin
	return m;
}
double delU(vector<double>& v, int size, int i, vector < vector <double> >& na)
{
	double E = 2 * v[i] * (v[na[i][0]] + v[na[i][1]] + v[na[i][2]] + v[na[i][3]]);
	return E;
}
double exp_delU(double E, double* expE)
{
	double result;
	if (E == 8) result = *(expE);
	else if (E == 4) result = *(expE + 1);
	else result = 1;
	return result;
}
void MC_1step(vector<double>& v, int size, double* expE, vector < vector <double> >& na)
{
	gen.seed(rd);
	for (int i = 0; i < size; i++) {
		int j = i % 2;
		for (int k = 0; k<int(size / 2); k++) {
			double Ediff = delU(v, size, size * i + j, na);
			if (Ediff <= 0) v[size * i + j] = -v[size * i + j];
			else {
				if (dis(gen) <= exp_delU(Ediff, expE)) v[size * i + j] = -v[size * i + j];
			}
			j = j + 2;
		}
	}
	for (int i = 1; i < size + 1; i++) {
		int j = i % 2;
		for (int k = 0; k<int(size / 2); k++) {
			double Ediff = delU(v, size, size * (i - 1) + j, na);
			if (Ediff <= 0) v[size * (i - 1) + j] = -v[size * (i - 1) + j];
			else {
				if (dis(gen) <= exp_delU(Ediff, expE)) v[size * (i - 1) + j] = -v[size * (i - 1) + j];
			}
			j = j + 2;
		}
	}

}
void MC_1cycle(int size, double T, vector < vector <double> >& na, vector<double>& magnet, vector<double>& magsus)
{
	int step1 = 5000, step2 = 150000;
	// int trash_step = 5 + size;
	// if (T > 2.0 && T < 2.5) trash_step = trash_step * 2;

	double expE[2] = { exp(-8 / T), exp(-4 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na); }

	double M=0, Mag = 0, Mag2 = 0;
	for (int k = 0; k < step2; k++) {
		MC_1step(array, size, &expE[0], na);
		M = Magnet(array, size);
		//Mag = Mag + M; Mag2 = Mag2 + pow(M, 2);
		magnet.at(k) = M;
		//magsus.at(k) = pow(size, 2)*(Mag2 / (k+1) - pow(Mag / (k+1), 2))/T;
	}
}

int main()
{
	int size=128;
	//cout << "Size: "; cin >> size;
	double temperature = 2.2692;

	clock_t start = clock();

	ofstream File;
	File.open("at_met_trng.txt");
	cout << "(trng) File open: " << size << endl;
	File << "trial magnet size: " << size << endl;

	vector < vector <double> > near(size * size, vector<double>(4, 0));
	vector<double> magnet(150000, 0); vector<double> magsus(10000,0);
	neighbor(near, size);
	for (int h = 0; h < 20; h++) {
		gen.seed(rd);
		MC_1cycle(size, temperature, near, magnet, magsus);
		for (int i=0;i<150000;i++){
			File << h << " " << magnet.at(i) << endl;
		}
		cout << h+1 << " trial end" << endl;
	}
	
	File.close();
	cout << "File closed " << endl;

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
