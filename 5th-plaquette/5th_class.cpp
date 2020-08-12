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

void neighbor(vector < vector <double> >& na, int size)
{
	int sizes = size * size;
	int na_2, na_3;
	for (int i = 0; i < sizes; i++) {
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
	}
}
class Sites
{
	public:
	Sites(vector<double>& v, vector < vector <double> > na, vector<double> J, int size, int nth)
	{ size_ = size; J_ = J; na_ = na; nth_ = nth; initialize(v); }
	double Magnet(vector<double> v_)
	{
		double m = 0;
		for (int i = 0; i < size_*size_; i++) { m = m + v_.at(i); }
		m = abs(m) / (size_*size_);
		return m;
	}
	double nnnEne(vector<double> v_, int ith)
	{
		double nnn = 0;
		for (int i = 0; i < size_* size_; i++){
			nnn = nnn - v_[i] * (v_[na_[i][4*ith-3]] + v_[na_[i][4*ith-1]]);
		}
		return nnn;
	}
	double Energy(vector<double> v_)
	{
		double e = J_[0];
		for (int i = 1; i < nth_ + 1; i++){ e = e + J_[i] * nnnEne(v_, i); }
		return e;
	}

	private:
	int size_, nth_; 
	vector < vector <double> > na_;
	vector<double> J_;
	void initialize(vector<double>& v_)
	{ for (int i = 0; i < size_*size_; i++) { v_[i] = dis(gen) < 0.5 ? 1 : -1; } }
};
class Cluster
{
	public:
	Cluster(vector < vector <double> > na, vector<double> padd, int size) 
	{ size_ = size; na_ = na; padd_ = padd; }
	void flip(vector<double>& v_)
	{
		int i = size_*size_ * dis(gen);
		int sp = 0, sh = 0, nnn = 0; 
		double point;
		double oldspin = v_[i], newspin = -v_[i]; 
		vector<int> stack(size_*size_, 0); vector<int> search(size_*size_, 0);
		stack[0] = i; search[i] = 1; v_[i] = newspin;
		while (1) {
			for (int k = 0; k < 4; k++){
				point = na_[i][k];
				nnn = search[na_[point][4]] + search[na_[point][5]] + search[na_[point][6]] + search[na_[point][7]];
				if (v_[point] == oldspin && dis(gen) < padd_[nnn]) { 
					sh++; 
					stack.at(sh) = point; search[point] = 1;
					v_[point] = newspin; 
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
void wolff_cycle(int size, double T, vector < vector <double> > na, vector<double> J,
vector<double> padd, vector<double>& magnet, vector<double>& energy, 
double Tstart, double clsizef, int nth, int step2)
{
	int step1 = 2500, scale = 1; 
	double slope = (double(size)*size/clsizef)/(5-Tstart);
	if (T>Tstart) { scale = slope * (T - Tstart); if (scale == 0) scale = 1; }
	int trash_step = scale*(sqrt(size));

	vector<double> array(size*size, 0);
	gen.seed(rd); 
	Sites s(array, na, J, size, nth);
	Cluster c(na, padd, size);

	for (int k = 0; k < step1*scale; k++) { gen.seed(rd); c.flip(array); }
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) { gen.seed(rd); c.flip(array); }
		magnet.at(k) = s.Magnet(array);
		energy.at(k) = s.Energy(array);
	}

}
void cal_variable(vector<double> magnet, vector<double> energy, int sizes, double T,
double& mag, double& mag_sus, double& mag2, double& mag4, double& ene, double& sp_heat, int step2)
{
	double Mag = 0, Mag2 = 0, Mag4 = 0, Ene = 0, Ene2 = 0;
	double a, b;
	for (int i = 0; i < step2; i++) {
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
	mag_sus = sizes * (Mag2 / step2 - (Mag/step2)*(Mag/step2)) / T;
	mag2 = Mag2 / step2;
	mag4 = Mag4 / step2;
	sp_heat = (Ene2 / step2 - (Ene/step2)*(Ene/step2)) / (sizes*T*T);
}
int main()
{
	int size = 10, step2 = 10000, nth = 2;
	vector<double> J(nth + 1, 0); J[1] = 1;
	vector<double> padd(5, 0);
	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	double Tstart = 2.3*J[1], clsizef = 1.86*J[1]*J[1]+1;
	double mag = 0, ms = 0, mag2 = 0, mag4 = 0, ene = 0, sp_h = 0, temp;
	vector < vector <double> > near(size * size, vector<double>(4*nth, 0));
	neighbor(near, size);

	clock_t start = clock();

	ofstream Fileout; 
	Fileout.open("class_test.txt");
	cout << "File test open: " << J[2] << endl;
	Fileout << "J2 temperature m m2 m4 ms e sp_h " << endl;
	for (int m = 0; m < 4; m++){
		J[2] = 0.15 - 0.1 * m;
		for (int k = 200; k < 700; k++) {
			temp = 0.005 * k;
			for (int i = 0; i < 5; i++){ padd[i] = 1 - exp(-2 * (J[1] + i * J[2]) / temp); }
			for (int h = 0; h < 2; h++) {
				wolff_cycle(size, temp, near, J, padd, magnet, energy, Tstart, clsizef, nth, step2);
				cal_variable(magnet, energy, size*size, temp, mag, ms, mag2, mag4, ene, sp_h, step2);
				Fileout << J[2] << " " << temp << " " << mag << " " << mag2 << " " << mag4 
				<< " " << ms << " " << ene << " " << sp_h << " " << endl;
			}
		}
	}
	Fileout.close();

	cout << "total time: " << (double(clock())-double(start))/CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}