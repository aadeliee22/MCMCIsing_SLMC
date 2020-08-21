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
void MC_1cycle(int size, double T, double& mag, double& mag_sus, double& mag2, double& mag4, vector < vector <double> >& na)
{
	int step1 = 2500, step2 = 10000;
	int trash_step =  size+8;
	if (T > 2.1 && T < 2.4) trash_step = trash_step * 2;

	double expE[2] = { exp(-8 / T), exp(-4 / T) };
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { MC_1step(array, size, &expE[0], na); }

	vector<double> magnet(step2, 0);
	//vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			MC_1step(array, size, &expE[0], na);
		}
		magnet.at(k) = Magnet(array, size);
		//energy.at(k) = Energy(array, size, na);
	}

	double Mag = 0, Mag2 = 0, Mag4 = 0;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		Mag = Mag + magnet.at(i);
		Mag2 = Mag2 + pow(magnet.at(i), 2);
		//Mag4 = Mag4 + pow(magnet.at(i), 4);
		/*Ene = Ene + energy.at(i);
		Ene2 = Ene2 + pow(energy.at(i), 2);*/
	}

	mag = Mag / step2;
	//ene = Ene / step2;
	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
	//mag2 = Mag2 / step2;
	//mag4 = Mag4 / step2;
	//sp_heat = (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2);
	//return Mag, Ene, mag_sus, sp_heat;
}

int main()
{
	random_device rd;
	gen.seed(rd);
	int size=80;
	//cout << "Center, What size?: ";	cin >> size;
	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0;

	clock_t start = clock();

	ofstream File;
	File.open("met_center_80.txt");
	cout << "(met_center) File open: " << size << endl;
	File << "size temperature m m^2 m^4 mag_sus" << endl;

	vector < vector <double> > near(size * size, vector<double>(4, 0));
	neighbor(near, size);
	// for (int k = 1; k < 75; k++) { // 0.02~1.48
	// 	for (int h = 0; h < 5; h++) {
	// 		MC_1cycle(size, 0.02 * k, Mag, mag_sus, Mag2, Mag4, near);
	// 		File << size << " " << 0.02 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
	// 	}
	// 	cout << 0.02 * k << "end" << endl;
	// }
	// for (int k = 150; k < 200; k++) { // 1.50~1.99
	// 	for (int h = 0; h < 5; h++) {
	// 		MC_1cycle(size, 0.01 * k, Mag, mag_sus, Mag2, Mag4, near);
	// 		File << size << " " << 0.01 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
	// 	}
	// 	cout << 0.01 * k << "end" << endl;
	// }
	for (int k = 2200; k < 2500; k++) { // 2.000~2.499
		for (int h = 0; h < 5; h++) {
			MC_1cycle(size, 0.001 * k, Mag, mag_sus, Mag2, Mag4, near);
			File << size << " " << 0.001 * k << " " << Mag << " "  << mag_sus << " " << endl;
		}
		cout << 0.001 * k << " end" << endl;
	}
	// for (int k = 250; k < 300; k++) { // 2.50~2.99
	// 	for (int h = 0; h < 5; h++) {
	// 		MC_1cycle(size, 0.01 * k, Mag, mag_sus, Mag2, Mag4, near);
	// 		File << size << " " << 0.01 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
	// 	}
	// 	cout << 0.01 * k << "end" << endl;
	// }
	// for (int k = 150; k < 226; k++) { // 3.00~4.5
	// 	for (int h = 0; h < 5; h++) {
	// 		MC_1cycle(size, 0.02 * k, Mag, mag_sus, Mag2, Mag4, near);
	// 		File << size << " " << 0.02 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
	// 	}
	// 	cout << 0.02 * k << "end" << endl;
	// }
	File.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
