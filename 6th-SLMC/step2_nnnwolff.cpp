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
double originalEnergy(vector<double>& v, int size, vector < vector <double> >& na, double K)
{
	double e = 0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]] + K * v[na[i][1]] * v[na[i][3]] * v[na[i][7]]);
	}
	return e;
}
double effEnergy(vector<double>& v, int size, vector < vector <double> >& na, vector<double>& J, int nth)
{
	double e = J[0];
	for (int i = 1; i < nth + 1; i++){
		e = e + J[i] * nnnEne(v, size, na, i);
	}
	return e;
}
void Cluster_1step(vector<double>& array, int size, vector<double>& padd, 
vector < vector <double> >& na, double T, double K, vector<double>& J, int nth)
{
	gen.seed(rd);
	vector<double> v = array;
	int i = size * size * dis(gen);
	vector<int> stack(1, i);
	double oldspin = v[i]; double newspin = -v[i]; v[i] = newspin;
	int sp = 0;
	while (1) {
		if (padd[0]!=0){
			if (v[na[i][0]] == oldspin && dis(gen) < padd[0]) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
			if (v[na[i][1]] == oldspin && dis(gen) < padd[0]) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
			if (v[na[i][2]] == oldspin && dis(gen) < padd[0]) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
			if (v[na[i][3]] == oldspin && dis(gen) < padd[0]) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
		}
		if (padd[1]!=0){
			if (v[na[i][4]] == oldspin && dis(gen) < padd[1]) { stack.push_back(na[i][4]); v[na[i][4]] = newspin; }
			if (v[na[i][5]] == oldspin && dis(gen) < padd[1]) { stack.push_back(na[i][5]); v[na[i][5]] = newspin; }
			if (v[na[i][6]] == oldspin && dis(gen) < padd[1]) { stack.push_back(na[i][6]); v[na[i][6]] = newspin; }
			if (v[na[i][7]] == oldspin && dis(gen) < padd[1]) { stack.push_back(na[i][7]); v[na[i][7]] = newspin; }
		}
		if (padd[2]!=0){
			if (v[na[i][8]] == oldspin && dis(gen) < padd[2]) { stack.push_back(na[i][8]); v[na[i][8]] = newspin; }
			if (v[na[i][9]] == oldspin && dis(gen) < padd[2]) { stack.push_back(na[i][9]); v[na[i][9]] = newspin; }
			if (v[na[i][10]] == oldspin && dis(gen) < padd[2]) { stack.push_back(na[i][10]); v[na[i][10]] = newspin; }
			if (v[na[i][11]] == oldspin && dis(gen) < padd[2]) { stack.push_back(na[i][11]); v[na[i][11]] = newspin; }
		}
		sp++;
		if (sp >= stack.size()) break;
		i = stack.at(sp);
	}
	double ediff = (originalEnergy(v, size, na, K) - effEnergy(v, size, na, J, nth)) - (originalEnergy(array, size, na, K) - effEnergy(array, size, na, J, nth));
	ediff = ediff / T;
	if (ediff <= 0) { array = v; }
	else {
		if (dis(gen) < exp(-ediff)) { array = v; }
	}
}
void wolff_cycle(int size, double T, vector < vector <double> >& na, double K, vector<double>& J,
vector<double>& energy, vector<double>& nn, vector<double>& nnn, vector<double>& nnnn, int nth)
{
	int step1 = 2500, step2 = 10000;
	int scale=1; double Tstart = 2.3 * J[1], clsizef = 1.86 * J[1] * J[1] + 1;
	if (nth > 1 && J[2]>0) { Tstart = Tstart + 0.2 * (J[2]/0.05); clsizef = clsizef + 1.6*J[2]/0.02; }
	double slope = (double(size)*size/clsizef)/(5-Tstart);
	if (T>Tstart) { scale = slope * (T - Tstart); if (scale == 0) scale = 1; }
	int trash_step = scale*(sqrt(size));

	vector<double> padd(3, 0);
	for (int i = 0; i < nth; i++){ padd.at(i) = 1 - exp(-2 * J[i+1] / T); }
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1*scale; k++) { Cluster_1step(array, size, padd, na, T, K, J, nth); }
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd, na, T, K, J, nth);
		}
		//magnet.at(k) = Magnet(array, size);
		energy.at(k) = originalEnergy(array, size, na, K);
		nn.at(k) = nnnEne(array, size, na, 1);
		nnn.at(k) = nnnEne(array, size, na, 2);
		nnnn.at(k) = nnnEne(array, size, na, 3);
	}
}
int main()
{
	random_device rd; gen.seed(rd);
	double K = 0.2;
	int nth; //How much nearest neighbors are we considering?
	int size = 20; double temp;
	//filein.txt format: temperature \n E0 \n J1 \n J2 \n J3
	ifstream Filein; Filein.open("filein.txt"); 
	Filein >> nth;
	vector<double> J(nth + 1, 0);
	for (int i = 0; i < nth + 1; i++){ Filein >> J[i]; }
	Filein >> temp;
	temp = temp - 0.2;

	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	vector<double> energy(10000, 0); 
	vector<double> nn(10000,0);
	vector<double> nnn(10000,0);
	vector<double> nnnn(10000,0);
	neighbor(near, size);

	ofstream Fileout; 
	Fileout.open("fileout.txt");
	cout << "Fileout open: " << temp << ", " << nth << endl;
	Fileout << "nth temp ene nn nnn nnnn " << endl;
	wolff_cycle(size, temp, near, K, J, energy, nn, nnn, nnnn, nth);
	for (int i = 0; i < 10000; i++){
		Fileout << nth << " " << temp << " " << energy.at(i) << " " << nn.at(i) << " " << nnn.at(i) << " " << nnnn.at(i) << endl;
	}
	Fileout.close();

	cout << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
