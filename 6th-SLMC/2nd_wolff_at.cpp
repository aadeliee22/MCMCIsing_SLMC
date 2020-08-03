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
double effEnergy(vector<double>& v, int size, vector < vector <double> >& na, vector<double>& J)
{
	double e = J[0];
	for (int i = 1; i < 4; i++){
		e = e + J[i] * nnnEne(v, size, na, i);
	}
	return e;
}
void Cluster_1step(vector<double>& array, int size, vector<double>& padd, 
vector < vector <double> >& na, double T, double K, vector<double>& J)
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
	double ediff = (originalEnergy(v, size, na, K) - effEnergy(v, size, na, J)) - (originalEnergy(array, size, na, K) - effEnergy(array, size, na, J));
	ediff = ediff / T;
	if (ediff <= 0) { array = v; }
	else {
		if (dis(gen) < exp(-ediff)) { array = v; }
	}
}
void wolff_cycle(int size, double T, vector < vector <double> >& na, double K, vector<double>& J,
vector<double>& padd, vector<double>& magnet)
{
	int step1 = 2500, step2 = 50000;
	int scale=1;
	int trash_step = scale*(sqrt(size));

	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1*scale; k++) { Cluster_1step(array, size, padd, na, T, K, J); }
	vector<double> magnet(step2, 0);
	for (int k = 0; k < step2; k++) {
		Cluster_1step(array, size, padd, na, T, K, J);
		magnet.at(k) = Magnet(array, size);
	}
}
int main()
{
	random_device rd; gen.seed(rd);
	double K = 0.2;
	vector<double> J(4);
	int size = 10; double temp=2.493;
	ifstream Filein; Filein.open("filein.txt"); 
	for (int i = 0; i < 4; i++){ Filein >> J[i]; }

	for (int k = 300; k < 700; k++) {
		vector<double> padd(3);
		for (int i = 0; i < 3; i++){ padd.at(i) = 1 - exp(-2 * J[i+1] / temp); }
	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	vector<double> magnet(50000, 0); 
	neighbor(near, size);

	ofstream Fileout; 
	Fileout.open("10_c_at.txt");
	cout << "File open: " << size << endl;
	Fileout << "size magnet " << endl;
	wolff_cycle(size, temp, near, K, J, padd, magnet);
	for (int i = 0; i < 50000; i++){
		Fileout << size << " " << magnet.at(i) << endl;
	}
	Fileout.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
