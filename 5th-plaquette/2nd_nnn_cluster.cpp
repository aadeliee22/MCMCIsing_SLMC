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
double Energy(vector<double>& v, int size, vector < vector <double> >& na, double J0, double J1, double J2)
{
	double e = J0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (J1 * (v[na[i][1]] + v[na[i][3]]) + J2 * v[na[i][1]] * v[na[i][3]] * v[na[i][7]]);
	}
	return e;
}
void Cluster_1step(vector<double>& v, int size, double padd1, double padd2, vector < vector <double> >& na)
{
	gen.seed(rd);
	int i = size * size * dis(gen);
	vector<int> stack(1, i);
	double oldspin = v[i];
	double newspin = -v[i];
	v[i] = newspin;
	int sp = 0;
	while (1) {
		if (v[na[i][0]] == oldspin && dis(gen) < padd1) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
		if (v[na[i][1]] == oldspin && dis(gen) < padd1) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
		if (v[na[i][2]] == oldspin && dis(gen) < padd1) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
		if (v[na[i][3]] == oldspin && dis(gen) < padd1) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
		if (v[na[i][4]] == oldspin && dis(gen) < padd2) { stack.push_back(na[i][4]); v[na[i][4]] = newspin; }
		if (v[na[i][5]] == oldspin && dis(gen) < padd2) { stack.push_back(na[i][5]); v[na[i][5]] = newspin; }
		if (v[na[i][6]] == oldspin && dis(gen) < padd2) { stack.push_back(na[i][6]); v[na[i][6]] = newspin; }
		if (v[na[i][7]] == oldspin && dis(gen) < padd2) { stack.push_back(na[i][7]); v[na[i][7]] = newspin; }
		sp++;
		if (sp >= stack.size()) break;
		i = stack.at(sp);
	}
}
void MC_1cycle(int size, double T, double& mag, double& mag_sus, double& mag2, double& mag4, vector < vector <double> >& na, double J0, double J1, double J2)
{
	int step1 = 2000, step2 = 10000;
	int scale=1;
	double slope = (double(size)*size/3.0)/2.7;
	if (T>2.3) {
		scale = slope * (T - 2.3);
		if (scale == 0) scale = 1;
	}
	int trash_step = scale*(sqrt(size));

	double padd1 = 1 - exp(-2 * J1 / T);
	double padd2 = 1 - exp(-2 * J2 / T);
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1*scale; k++) { Cluster_1step(array, size, padd1, padd2, na); }

	vector<double> magnet(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd1, padd2, na);
		}
		magnet.at(k) = Magnet(array, size);
	}

	double Mag = 0, Mag2 = 0, Mag4 = 0;
	double a;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		a = magnet.at(i);
		Mag = Mag + a;
		Mag2 = Mag2 + a*a;
		Mag4 = Mag4 + a*a*a*a;
	}

	mag = Mag / step2;
	mag_sus = size*size * (Mag2 / step2 - (Mag/step2)*(Mag/step2)) / T;
	mag2 = Mag2 / step2;
	mag4 = Mag4 / step2;
}
void MC_1cycle_graphing(int size, double T, vector < vector <double> >& na, double J0, double J1, double J2)
{
	int step1 = 1000, step2 = 5000;
	int scale=1;
	double slope = (double(size)*size/3.0)/2.7;
	if (T>2.3) {
		scale = slope * (T - 2.3);
		if (scale == 0) scale = 1;
	}
	int trash_step = scale*(sqrt(size));

	double padd1 = 1 - exp(-2 * J1 / T);
	double padd2 = 1 - exp(-2 * J2 / T);
	vector<double> array(size * size, 0);

	initialize(array, size);
	cout << "Initial state: ";
	color(array, size);
	cout << endl;
	cout << "Magnetization: " << Magnet(array, size) << endl;
	cout << "Energy (H): " << Energy(array, size, na, J0, J1, J2) << endl;

	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd1, padd2, na); }

	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++){
		Cluster_1step(array, size, padd1, padd2, na);
		magnet.at(k) = Magnet(array, size);
		energy.at(k) = Energy(array, size, na, J0, J1, J2);
	}

	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
	double a, b;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		a = magnet.at(i);
		b = energy.at(i);
		Mag = Mag + a;
		Mag2 = Mag2 + a*a;
		Ene = Ene + b;
		Ene2 = Ene2 + b*b;
	}

	cout << endl << "Final state with temperature " << T <<" :";
	color(array, size);
	cout << endl;
	cout << "Magnetization: " << Mag/step2 << endl;
	cout << "Energy (H): " << Ene/step2 << endl;
	cout << "Magnetic susceptibility: " << size*size * (Mag2 / step2 - (Mag/step2)*(Mag/step2)) / T << endl;
	cout << "Specific heat: " << (Ene2 / step2 - (Ene/step2)*(Ene/step2)) / (size*size*T*T) << endl;
}
int main()
{
	random_device rd;
	gen.seed(rd);
	double J0 = 0, J1 = 1, J2 = 0.1;
	int size; double temp;
	cout << "What size?: ";	cin >> size;
	cout << "What T?: ";	cin >> temp;

	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0;

	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(8, 0));
	neighbor(near, size);
	MC_1cycle_graphing(size, temp, near, J0, J1, J2);


	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
