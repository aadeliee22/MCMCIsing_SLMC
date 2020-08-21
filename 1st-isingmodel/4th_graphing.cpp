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
void color(vector<double> v, int size) //graphing state
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
double Magnet(vector<double> v, int size)
{
	double m = 0;
	for (vector<int>::size_type i = 0; i < v.size(); i++) {
		m = m + v.at(i);
	}
	m = abs(m) / (v.size()); //absolute value of average spin
	return m;
}
double delU(vector<double> v, int size, int i, vector < vector <double> > na)
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
void MC_1step(vector<double>& v, int size, double* expE, vector < vector <double> > na)
{
	gen.seed(rd);
	int point;
	for (int i = 0; i < size; i++) {
		int j = i % 2;
		for (int k = 0; k<int(size / 2); k++) {
			point = size * i + j;
			double Ediff = delU(v, size, point, na);
			if (Ediff <= 0) v[point] = -v[point];
			else {
				if (dis(gen) <= exp_delU(Ediff, expE)) v[point] = -v[point];
			}
			j = j + 2;
		}
	}
	for (int i = 1; i < size + 1; i++) {
		int j = i % 2;
		for (int k = 0; k<int(size / 2); k++) {
			point = size * (i - 1) + j;
			double Ediff = delU(v, size, point, na);
			if (Ediff <= 0) v[point] = -v[point];
			else {
				if (dis(gen) <= exp_delU(Ediff, expE)) v[point] = -v[point];
			}
			j = j + 2;
		}
	}

}
void MC_1cycle_graphing(int size, double T, vector<double>& array, vector < vector <double> > na)
{
	int step1 = 5000, step2 = 10000;
	int trash_step = size*size/50;
	if (T>2.2 && T <2.4) trash_step = trash_step*2;
	double expE[2] = { exp(-8 / T), exp(-4 / T) };

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na); }
	//for (int k = 0; k < step2*trash_step; k++) {MC_1step(array, size, &expE[0], na);}

	cout << "Final state with temperature " << T << " :" << endl;
	color(array, size);
}
int main()
{
	random_device rd;
	gen.seed(rd);
	int size=20;
	double temp;
	cout << "T: "; cin >> temp;
	vector<double> array(size * size, 0);
	initialize(array, size);

	clock_t start = clock();

	ofstream File;
	File.open("graph.txt");
	cout << "File open: " << size << endl;
	File << "size: " << size << endl;
	File << "temp: " << temp << endl;
	File << "s " << endl;

	vector < vector <double> > near(size * size, vector<double>(4, 0));
	neighbor(near, size);

	MC_1cycle_graphing(size, temp, array, near);
	for (int i = 0; i < size*size; i++){
		File << array[i] << endl;
	}
	
	File.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
