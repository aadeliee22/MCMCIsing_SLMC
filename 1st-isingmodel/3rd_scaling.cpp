#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <time.h>
using namespace std;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dis(0, 1);

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
/*double Energy(vector<double>& v, int size, vector < vector <double> >& na)
{
	double e = 0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]]);
	}
	return e;
}*/
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
	int step1 = 2000, step2 = 10000;
	int trash_step = 5 + size;
	if (T > 2.0 && T < 2.5) trash_step = trash_step * 2;

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
		Mag4 = Mag4 + pow(magnet.at(i), 4);
		/*Ene = Ene + energy.at(i);
		Ene2 = Ene2 + pow(energy.at(i), 2);*/
	}

	mag = Mag / step2;
	//ene = Ene / step2;
	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
	mag2 = Mag2 / step2;
	mag4 = Mag4 / step2;
	//sp_heat = (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2);
	//return Mag, Ene, mag_sus, sp_heat;
}

int main()
{
	clock_t start = clock();

	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0;
	int size = 0;
	ofstream File;
	File.open("ising.txt");
	cout << "File open" << endl;
	File << "size temperature m m^2 m^4 mag_sus" << endl;

	for (int s = 0; s < 4; s++) {
		size = 16 * (1 + s);
		vector < vector <double> > near(size * size, vector<double>(4, 0));
		neighbor(near, size);
		for (int k = 75; k < 151; k++) {
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.02 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.02 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.1 * k << "K end" << endl;
		}
		cout << "size: " << size << " finished :)" << endl;
	}
	
	File.close();
	cout << "File closed " << endl;

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
