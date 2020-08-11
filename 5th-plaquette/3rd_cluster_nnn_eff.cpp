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
void Cluster_1step(vector<double>& array, int size, vector<double>& padd, 
vector < vector <double> >& na, double T, vector<double>& J, int nth)
{
	gen.seed(rd);
	vector<double> v = array;
	int i = size * size * dis(gen);
	vector<int> stack(size*size, 0);
	stack[0] = i;
	double oldspin = v[i]; double newspin = -v[i]; v[i] = newspin;
	int sp = 0, sh = 0;
	while (1) {
		if (v[na[i][0]] == oldspin && dis(gen) < padd[0]) { sh++; stack.at(sh) = na[i][0]; v[na[i][0]] = newspin; }
		if (v[na[i][1]] == oldspin && dis(gen) < padd[0]) { sh++; stack.at(sh) = na[i][1]; v[na[i][1]] = newspin; }
		if (v[na[i][2]] == oldspin && dis(gen) < padd[0]) { sh++; stack.at(sh) = na[i][2]; v[na[i][2]] = newspin; }
		if (v[na[i][3]] == oldspin && dis(gen) < padd[0]) { sh++; stack.at(sh) = na[i][3]; v[na[i][3]] = newspin; }
		sp++;
		if (sp > sh) break;
		i = stack.at(sp);
	}
	double ediff = (Energy(v, size, na, J, nth) - Energy(v, size, na, J, 1)) - (Energy(array, size, na, J, nth) - Energy(array, size, na, J, 1));
	ediff = ediff / T;
	if (ediff <= 0) { array = v; }
	else {
		if (dis(gen) < exp(-ediff)) { array = v; }
	}
}
void wolff_cycle(int size, double T, double& mag, double& mag_sus, double& mag2, double& mag4,
double& ene, double& sp_heat, vector < vector <double> >& na, vector<double>& J,
vector<double>& padd, double Tstart, double clsizef, int nth)
{
	int step1 = 2500, step2 = 10000;
	int scale=1; 
	double slope = (double(size)*size/clsizef)/(5-Tstart);
	if (T>Tstart) { scale = slope * (T - Tstart); if (scale == 0) scale = 1; }
	int trash_step = scale*(sqrt(size));

	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1*scale; k++) { Cluster_1step(array, size, padd, na, T, J, nth); }
	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd, na, T, J, nth);
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
	int size = 10; int nth = 2;
	vector<double> J(nth + 1, 0);
	J[0] = 0; J[1] = 1;

	double Tstart = 2.3 * J[1], clsizef = 1.86 * J[1] * J[1] + 1;
	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0, Ene = 0, sp_heat = 0;
	double temp = 0;
	
	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(8, 0));
	neighbor(near, size);

	ofstream Fileout; 
	Fileout.open("cluster_eff_2.txt");
	cout << "File eff open: " << J[2] << endl;
	Fileout << "J2 temperature m m2 m4 ms e sp_h " << endl;
	vector<double> padd(1, 0);
	for (int m = 0; m < 9; m++){
		J[2] = 0.2 - 0.05 * m;
		for (int k = 200; k < 700; k++) {
			temp = 0.005 * k;
			padd[0] = 1 - exp(-2 * J[1] / temp);
			for (int h = 0; h < 2; h++) {
				wolff_cycle(size, temp, Mag, mag_sus, Mag2, Mag4, Ene, sp_heat, near, J, padd, Tstart, clsizef, nth);
				Fileout << J[2] << " " << temp << " " << Mag << " " << Mag2 << " " << Mag4 
				<< " " << mag_sus << " " << Ene << " " << sp_heat << " " << endl;
			}
		}
		cout << J[2] << " end" << endl;
	}
	Fileout.close();

	cout << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}