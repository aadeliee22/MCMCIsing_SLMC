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
void Cluster_1step(vector<double>& v, int size, double padd1, double padd2, vector < vector <double> >& na, double& clustersize)
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
	clustersize = stack.size();
}
void MC_1cycle_clustersize(int size, double T, vector < vector <double> >& na, double K, double& cl_size)
{
	int step1 = 5000, step2 = 5000;
	double padd1 = 1 - exp(-2 * K / T);
	double padd2 = 0;
	double clustersize;
	vector<double> array(size * size, 0);
	initialize(array, size);
	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd1, padd2, na, clustersize); }
	for (int k = 0; k < step2; k++) {
		Cluster_1step(array, size, padd1, padd2, na, clustersize);
		cl_size = clustersize + cl_size;
	}
	cl_size = cl_size / step2;
	//cout << "Final state with temperature " << T << " :" << endl;
	//color(array, size);
}
int main()
{
	random_device rd;
	gen.seed(rd);
	double K;
	int size=20;
	//cout << "What size?: ";	cin >> size;

	clock_t start = clock();

	ofstream File;
	File.open("clustersize.txt");
	File << "temp cl_size size: " << size << " " << endl;
	vector < vector <double> > near(size * size, vector<double>(8, 0));
	neighbor(near, size);
	for (int j = 0; j < 9; j++){
		K = 0.8 + 0.05*double(j);
		for (int k = 10; k < 51; k++) {
			double cl_size = 0;
			MC_1cycle_clustersize(size, 0.1*k, near, K, cl_size);
			File << 0.1 * k << " " << cl_size << endl;
		}
		cout << K << " done" << endl;
	}
	File.close();

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
