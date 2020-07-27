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
	int step1 = 5000, step2 = 100000;
	// int trash_step = 5 + size;
	// if (T > 2.0 && T < 2.5) trash_step = trash_step * 2;

	double expE[2] = { exp(-8 / T), exp(-4 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na); }

	double M = 0, Mag = 0, Mag2 = 0, magsum=0;
	for (int k = 0; k < step2; k++) {
		MC_1step(array, size, &expE[0], na);
		M = Magnet(array, size);
		magsum = magsum + M;
		magnet.at(k) = M;
	}
	for (int k=0;k<step2;k++){
		magsus.at(k)=pow(size, 2) * pow(magnet.at(k)-magsum/step2, 2)/T;
	}
}
double jack_error(int size, vector<double>& target)
{
	int Nall = target.size();
	int binsize = 5 * pow(size, 0.5);
	int binnum = Nall / binsize;
	double target_avg;
	for (int i = 0; i < Nall; i++){ target_avg = target_avg + target.at(i); }
	target_avg = target_avg / Nall;

	vector<double> Mbin(binnum,0);
	for (int j = 0; j < binnum; j++){
		double temp = 0;
		for (int k = 0; k < binsize; k++){
			temp = temp + target.at(j * binsize + k);
		}
		Mbin.at(j) = temp / binsize;
	}
	double Mbinsum = 0;
	for (int h = 0; h < binnum; h++){
		Mbinsum = Mbinsum + Mbin.at(h);
	}

	double inside = 0;
	for (int k = 0; k < binnum; k++){
		inside = inside + pow( (Mbinsum - Mbin.at(k)) / (binnum - 1) - target_avg, 2);
	}
	inside = sqrt(inside * (binnum - 1) / binnum);
	return inside;
}

int main()
{
	int size;
	//cout << "Size: "; cin >> size;
	//double temperature = 2.2692;

	clock_t start = clock();

	ofstream File;
	File.open("met_error.txt");
	//cout << "(cluster) File open: " << size << endl;
	File << "size temperature err_m err_ms" << endl;

	for (int s = 0; s < 5; s++){
		size = 16 * (1 + s);
		vector < vector <double> > near(size * size, vector<double>(4, 0));
		vector<double> magnet(100000, 0); vector<double> magsus(100000,0);
		neighbor(near, size);
		for (int t = 2200; t < 2500; t++) {
			for (int i=0;i<1;i++){
				double err_m=0, err_ms=0;
				gen.seed(rd);
				MC_1cycle(size, 0.001*t, near, magnet, magsus);
				err_m = jack_error(size, magnet);
				err_ms = jack_error(size, magsus);
				File << size << " " << 0.001*t << " " << err_m << " " << err_ms << endl;
				}
			cout << 0.001*t << " end" << endl;
		}
		cout << "size " << size << " finished" << endl;
	}
	cout << "File closed " << endl;

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
