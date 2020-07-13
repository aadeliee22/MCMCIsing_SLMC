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
	for (int sp = 0; sp < size * size; sp++) {
		try {
			if (v[na[i][0]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
			if (v[na[i][1]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
			if (v[na[i][2]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
			if (v[na[i][3]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
			i = stack.at(sp+1);
		}
		catch (out_of_range) { break; }
	}
	/*int sp = 1;
	while (sp!=0) {
		int j = stack[--sp];
		if (v[na[j][0]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][0]); v[na[j][0]] = newspin; sp++; } }
		if (v[na[j][1]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][1]); v[na[j][1]] = newspin; sp++; } }
		if (v[na[j][2]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][2]); v[na[j][2]] = newspin; sp++; } }	
		if (v[na[j][3]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][3]); v[na[j][3]] = newspin; sp++; } }
	}*/
	
}
void MC_1cycle(int size, double T, double& mag, double& ene, double& mag_sus, double& sp_heat, vector < vector <double> >& na)
{
	int step1 = 2000, step2 = 10000;
	int trash_step = int(size/5);
	if (T > 2.0 && T < 2.5) trash_step = trash_step*2;

	double padd = 1 - exp(-2 / T);
	vector<double> array(size * size, 0);
	initialize(array, size);

	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }

	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd, na);
		}
		magnet.at(k) = Magnet(array, size);
		energy.at(k) = Energy(array, size, na);
	}

	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		Mag = Mag + magnet.at(i);
		Mag2 = Mag2 + pow(magnet.at(i), 2);
		Ene = Ene + energy.at(i);
		Ene2 = Ene2 + pow(energy.at(i), 2);
	}

	mag = Mag / step2;
	ene = Ene / step2;
	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
	sp_heat = (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2);
	//return Mag, Ene, mag_sus, sp_heat;
}
void MC_1cycle_graphing(int size, double T, vector < vector <double> >& na)
{
	int step1 = 2000, step2 = 2000;
	int trash_step = 2;
	if (T > 2.0 && T < 2.5) trash_step = 10;
	double padd = 1 - exp(-2 / T);
	vector<double> array(size * size, 0);

	initialize(array, size);
	cout << "Initial state: " << endl;
	color(array, size);
	cout << endl;

	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }

	vector<double> magnet(step2, 0);
	vector<double> energy(step2, 0);
	for (int k = 0; k < step2; k++) {
		for (int h = 0; h < trash_step; h++) {
			Cluster_1step(array, size, padd, na);
		}
		magnet.at(k) = Magnet(array, size);
		energy.at(k) = Energy(array, size, na);
	}

	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
		Mag = Mag + magnet.at(i);
		Mag2 = Mag2 + pow(magnet.at(i), 2);
		Ene = Ene + energy.at(i);
		Ene2 = Ene2 + pow(energy.at(i), 2);
	}

	cout << "Final state with temperature " << T << " :" << endl;
	color(array, size);
	cout << endl;
	cout << "Magnetization: " << Mag / step2 << endl;
	cout << "Energy (H): " << Ene / step2 << endl;
	cout << "Magnetic susceptibility: " << pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T << endl;
	cout << "Specific heat: " << (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2) << endl;
}

int main()
{
	int size;
	double temperature;

	//cout << "Size: "; cin >> size;
	//cout << "Temperature: "; cin >> temperature;

	clock_t start = clock();

	double Mag = 0, Ene = 0, mag_sus = 0, sp_heat = 0;
	ofstream File;
	File.open("wolff.txt");
	cout << "File open '^'" << endl;
	File << "size temperature mag energy mag_sus sp_heat" << endl;
	for (int s = 0; s < 3; s++) {
		size = 10 * (1 + s);
		vector < vector <double> > near(size * size, vector<double>(4, 0));
		neighbor(near, size);
		for (int k = 75; k < 151; k++) {
			for (int h = 0; h < 5; h++) {
				MC_1cycle(size, 0.02 * k, Mag, Ene, mag_sus, sp_heat, near);
				File << size << " " << 0.02 * k << " " << Mag << " " << Ene << " " << mag_sus << " " << sp_heat << " " << endl;
			}
			cout << 0.02*k << "end"<<endl;
		}
		cout << "size: " << size << " finished :)" << endl;
	}
	File.close();

	//MC_1cycle_graphing(size, temperature, near);

	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
	return 0;
}
