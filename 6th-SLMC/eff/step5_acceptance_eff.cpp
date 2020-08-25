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
	for (int i = 0; i < size * size; i++) {v[i] = dis(gen) < 0.5 ? 1 : -1;}
}
void neighbor(vector < vector <double> >& na, int size)
{
	int sizes = size * size;
	int na_2, na_3;
	for (int i = 0; i < size * size; i++) {
		na_2 = (i - 1 + size) % size + (i / size) * size;
		na_3 = (i + 1) % size + (i / size) * size;
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
double Magnet(vector<double> v)
{
	double m = 0;
	for (vector<int>::size_type i = 0; i < v.size(); i++) { m = m + v.at(i);}
	m = abs(m) / (v.size()); //absolute value of average spin
	return m;
}
double nnnEne(vector<double> v, int size, vector < vector <double> > na, int ith)
{
	double nnn = 0;
	for (int i=0; i<size*size; i++){
		nnn = nnn - v[i] * (v[na[i][4*ith-3]] + v[na[i][4*ith-1]]);
	}
	return nnn;
}
double originalEnergy(vector<double> v, int size, vector < vector <double> > na, double K)
{
	double e = 0;
	for (int i=0; i<size*size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]] + K * v[na[i][1]] * v[na[i][3]] * v[na[i][7]]);
	}
	return e;
}
double effEnergy(vector<double> v, int size, vector < vector <double> > na, vector<double> J, int nth)
{
	double e = J[0];
	for (int i = 1; i < nth + 1; i++){e = e + J[i] * nnnEne(v, size, na, i);}
	return e;
}
class Cluster
{
	public:
	Cluster(vector < vector <double> > na, vector<double> J, double K, double padd, 
		int size, double T, int nth) 
	{ size_ = size; na_ = na; padd_ = padd; J_ = J; K_ = K; T_ = T; nth_ = nth; }
	void flip(vector<double>& array_, int& ac_1, int& ac_a)
	{
		vector<double> v = array_;
		int i = size_*size_ * dis(gen);
		int sp = 0, sh = 0; 
		double oldspin = v[i], newspin = -v[i]; 
		vector<int> stack(size_*size_, 0);
		stack[0] = i; v[i] = newspin;
		while (1) {
		if (v[na_[i][0]] == oldspin && dis(gen) < padd_) { sh++; stack.at(sh) = na_[i][0]; v[na_[i][0]] = newspin; }
		if (v[na_[i][1]] == oldspin && dis(gen) < padd_) { sh++; stack.at(sh) = na_[i][1]; v[na_[i][1]] = newspin; }
		if (v[na_[i][2]] == oldspin && dis(gen) < padd_) { sh++; stack.at(sh) = na_[i][2]; v[na_[i][2]] = newspin; }
		if (v[na_[i][3]] == oldspin && dis(gen) < padd_) { sh++; stack.at(sh) = na_[i][3]; v[na_[i][3]] = newspin; }
		sp++;
		if (sp > sh) break;
		i = stack.at(sp);
		}
		double ediff1, ediff2;
		ediff1 = (effEnergy(v, size_, na_, J_, nth_) - effEnergy(v, size_, na_, J_, 1)) 
			- (effEnergy(array_, size_, na_, J_, nth_) - effEnergy(array_, size_, na_, J_, 1));
		ediff1 = ediff1 / T_;
		if ((ediff1 <= 0) || (dis(gen) < exp(-ediff1))) { 
			ac_1++;
			ediff2 = (originalEnergy(v, size_, na_, K_) - effEnergy(v, size_, na_, J_, nth_)) 
				- (originalEnergy(array_, size_, na_, K_) - effEnergy(array_, size_, na_, J_, nth_));
			ediff2 = ediff2 / T_;
			if ((ediff2 <= 0) || (dis(gen) < exp(-ediff2))) { ac_a++; array_ = v; }
		}
	}

	private:
	int size_, nth_;
	double K_, T_, padd_; 
	vector<double> J_; 
	vector < vector <double> > na_;
};
void wolff_cycle(int size, double T, int step2, int nth, 
vector < vector <double> > na, double K, vector<double> J, double padd, 
int& ac_1, int& ac_a)
{
	int step1 = 2500;
	int scale=1; double Tstart = 2.3 * J[1], clsizef = 1.86 * J[1] * J[1] + 1;
	double slope = (double(size)*size/clsizef)/(5-Tstart);
	if (T>Tstart) { scale = slope * (T - Tstart); if (scale == 0) scale = 1; }
	int trash_step = scale*(sqrt(size));

	
	vector<double> array(size * size, 0);
	initialize(array, size);
	Cluster c(na, J, K, padd, size, T, nth);

	for (int k = 0; k < step1*scale; k++) { gen.seed(rd); c.flip(array, ac_1, ac_a); }
	ac_1 = 0; ac_a = 0;
	for (int k = 0; k < step2; k++) { gen.seed(rd); c.flip(array, ac_1, ac_a); }
}
int main()
{
	random_device rd; gen.seed(rd);
	double K = 0.2; double temp = 2.493;
	int nth = 1; int size = 128; int step2 = 512;
	ifstream Filein; Filein.open("filein.txt");
	vector<double> J(4, 0);
	for (int i = 0; i < nth + 1; i++){ Filein >> J[i]; }
	double padd = 1 - exp(-2 * J[1] / temp);

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	neighbor(near, size);

	clock_t start = clock();

	ofstream Fileout; 
	Fileout.open("fileout_eff.txt");
	cout << "Fileout open: " << size<< " " <<nth << endl;
	Fileout << "s nth step2 ac1 aca " << endl;
	for (int i = 0; i < 3; i++){
		int ac_1 = 0, ac_a = 0;
		step2 = step2 * 2;
		wolff_cycle(size, temp, step2, nth, near, K, J, padd, ac_1, ac_a);
		Fileout << size << " " << nth << " " << step2 << " " << ac_1 << " " << ac_a << endl;
	}
	Fileout.close();

	cout << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
