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
double Energy(vector<double>& v, int size, vector < vector <double> >& na)
{
	double e = 0;
	for (int i = 0; i < size * size; i++) {
		e = e - v[i] * (v[na[i][1]] + v[na[i][3]]);
	}
	return e;
}
void Cluster_1step(vector<double>& v, int size, double padd, vector < vector <double> >& na)
{
	int i = size * size * dis(gen);
	vector<int> stack(1, i);
	double oldspin = v[i];
	double newspin = -v[i];
	v[i] = newspin;
	int sp = 0;
	while (1) {
		if (v[na[i][0]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
		if (v[na[i][1]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
		if (v[na[i][2]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
		if (v[na[i][3]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
		sp++;
		if (sp >= stack.size()) break;
		i = stack.at(sp);
	}
	
}
void MC_1cycle(int size, double T, double& mag, double& mag_sus, double& mag2, double& mag4, vector < vector <double> >& na)
{
	int step1 = 2000, step2 = 10000;
	int scale;
	if (T <= 2) scale = 1;
	if (T > 2 && T <= 2.6) scale = int(size * size / 40);
	if (T > 2.6 && T <= 3.2) scale = int(size * size / 10);
	if (T > 3.2 && T <= 3.8) scale = int(size * size / 5);
	if (T > 3.8 && T <= 4.4) scale = int(size * size / 4);
	if (T > 4.4) scale = int(size * size / 2);
	step1 = step1 * scale; step2 = step2 * scale;

	int trash_step = int(size / 4);
	if (T > 2.0 && T <= 2.5) trash_step = trash_step * 2;

	double padd = 1 - exp(-2 / T);
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }

	vector<double> magnet(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd, na);
		}
		magnet.at(k) = Magnet(array, size);
	}

	double Mag = 0, Mag2 = 0, Mag4 = 0;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		Mag = Mag + magnet.at(i);
		Mag2 = Mag2 + pow(magnet.at(i), 2);
		Mag4 = Mag4 + pow(magnet.at(i), 4);
	}

	mag = Mag / step2;
	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
	mag2 = Mag2 / step2;
	mag4 = Mag4 / step2;
}
int main()
{
	int size;
	double temperature;
	double Mag = 0, mag_sus = 0, Mag2 = 0, Mag4 = 0;

	clock_t start = clock();

	ofstream File;
	File.open("wolff_scale.txt");
	cout << "File open '^'" << endl;
	File << "size temperature m m^2 m^4 mag_sus" << endl;
	for (int s = 0; s < 4; s++) {
		size = 16 * (1 + s);
		vector < vector <double> > near(size * size, vector<double>(4, 0));
		neighbor(near, size);
		for (int k = 1; k < 75; k++) { // 0.02~1.48
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.02 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.02 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.02 * k << "end" << endl;
		}
		for (int k = 150; k < 200; k++) { // 1.50~1.99
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.01 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.01 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.01 * k << "end" << endl;
		}
		for (int k = 2000; k < 2500; k++) { // 2.000~2.499
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.001 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.001 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.001 * k << "end" << endl;
		}
		for (int k = 250; k < 300; k++) { // 2.50~2.99
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.01 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.01 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.01 * k << "end" << endl;
		}
		for (int k = 150; k < 251; k++) { // 3.00~5.00
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.02 * k, Mag, mag_sus, Mag2, Mag4, near);
				File << size << " " << 0.02 * k << " " << Mag << " " << Mag2 << " " << Mag4 << " " << mag_sus << " " << endl;
			}
			cout << 0.02 * k << "end" << endl;
		}
		cout << "size: " << size << " finished :)" << endl;
	}
	File.close();


	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
