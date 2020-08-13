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
double nnnEne(vector<double> v, int size, vector < vector <double> > na, int ith)
{
	double nnn = 0;
	for (int i = 0; i < size * size; i++){
		nnn = nnn - v[i] * (v[na[i][4*ith-3]] + v[na[i][4*ith-1]]);
	}
	return nnn;
}
double originalEnergy(vector<double> v, int size, vector < vector <double> > na, double K)
{
	double e = 0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]] + K * v[na[i][1]] * v[na[i][3]] * v[na[i][7]]);
	}
	return e;
}
double effEnergy(vector<double> v, int size, vector < vector <double> > na, vector<double> J, int nth)
{
	double e = J[0];
	for (int i = 1; i < nth + 1; i++){
		e = e + J[i] * nnnEne(v, size, na, i);
	}
	return e;
}
class Cluster
{
	public:
	Cluster(vector < vector <double> > na, vector<double> padd, int size) 
	{ size_ = size; na_ = na; padd_ = padd; }
	void flip(vector<double>& v_)
	{
		int i = size_*size_ * dis(gen);
		int sp = 0, sh = 0, n2 = 0, n3 = 0; 
		double point; double prob;
		double oldspin = v_[i], newspin = -v_[i]; 
		vector<int> stack(size_*size_, 0); vector<int> search(size_*size_, 0);
		stack[0] = i; search[i] = 1; v_[i] = newspin;
		while (1) {
			for (int k = 0; k < 4; k++){
				point = na_[i][k];
				n2 = search[na_[point][4]] + search[na_[point][5]] + search[na_[point][6]] + search[na_[point][7]];
				n3 = search[na_[point][8]] + search[na_[point][9]] + search[na_[point][10]] + search[na_[point][11]];
				prob = padd_[5*n2+n3];
				if (prob>0){
					if (v_[point] == oldspin && dis(gen) < prob) { 
						sh++; 
						stack.at(sh) = point; search[point] = 1;
						v_[point] = newspin; 
					}
				}
			}
			sp++;
			if (sp > sh) break;
			i = stack.at(sp);
		}
	}

	private:
	int size_;
	vector<double> padd_; vector < vector <double> > na_;
};
void wolff_cycle(int size, double T, 
vector < vector <double> > na, double K, vector<double> J, vector<double> padd,
vector<double>& energy, vector<double>& nn, vector<double>& nnn, vector<double>& nnnn, int nth)
{
	int step1 = 2500, step2 = 10000;
	int scale=1; double Tstart = 2.3 * J[1], clsizef = 1.86 * J[1] * J[1] + 1;
	double slope = (double(size)*size/clsizef)/(5-Tstart);
	if (T>Tstart) { scale = slope * (T - Tstart); if (scale == 0) scale = 1; }
	int trash_step = scale*(sqrt(size));

	vector<double> array(size * size, 0);
	initialize(array, size);
	Cluster c(na, padd, size);

	for (int k = 0; k < step1*scale; k++) { gen.seed(rd); c.flip(array); }
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			gen.seed(rd); c.flip(array);
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
	int nth; int size = 10; double temp;
	//filein.txt format: nth \n temperature \n E0 \n J1 \n J2 \n J3
	ifstream Filein; Filein.open("filein_srch.txt"); 
	Filein >> nth; 
	Filein >> temp;
	vector<double> J(4, 0);
	for (int i = 0; i < nth + 1; i++){ Filein >> J[i]; }
	temp = temp - 0.2;

	clock_t start = clock();

	vector < vector <double> > near(size * size, vector<double>(12, 0));
	neighbor(near, size);
	vector<double> energy(10000, 0); 
	vector<double> nn(10000,0);
	vector<double> nnn(10000,0);
	vector<double> nnnn(10000,0);
	vector<double> padd(25, 0);
	for (int i = 0; i < 5; i++){ 
		for (int j = 0; j < 5; j++){ 
			padd[5*i + j] = 1 - exp(-2 * (J[1] + i * J[2] + j * J[3]) / temp); 
		}
	}

	ofstream Fileout; 
	Fileout.open("fileout_srch.txt");
	cout << "Fileout open: " << temp << ", " << nth << endl;
	Fileout << "nth temp ene nn nnn nnnn " << endl;
	wolff_cycle(size, temp, near, K, J, padd, energy, nn, nnn, nnnn, nth);
	for (int i = 0; i < 10000; i++){
		Fileout << nth << " " << temp << " " << energy.at(i) << " " << nn.at(i) << " " << nnn.at(i) << " " << nnnn.at(i) << endl;
	}
	Fileout.close();

	cout << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
