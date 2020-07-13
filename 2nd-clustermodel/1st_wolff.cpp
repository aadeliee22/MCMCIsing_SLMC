//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <random>
//#include <vector>
//#include <time.h>
//using namespace std;
//
//random_device rd;
//mt19937 gen(rd());
//uniform_real_distribution<double> dis(0, 1);
//
//void initialize(vector<double>& v, int size) //initial -random- state
//{
//	for (int i = 0; i < size; i++) {
//		for (int j = 0; j < size; j++) {
//			if ((i + j) % 2 == 0) v[size * i + j] = 1;
//			else v[size * i + j] = -1;
//		}
//	}
//	/*for (int i = 0; i < size * size; i++) {
//		v[i] = dis(gen) < 0.5 ? 1 : -1;
//	}*/
//}
//void color(vector<double>& v, int size) //graphing state
//{
//	for (int i = 0; i < size * size; i++) {
//		if (i % size == 0) std::cout << endl;
//		if (v[i] == 1) std::cout << "* ";
//		if (v[i] == -1) std::cout << ". ";
//	}
//}
//void neighbor(vector < vector <double> >& na, int size)
//{
//	int sizes = size * size;
//	for (int i = 0; i < size * size; i++) {
//		na[i][0] = (i + size * (size - 1)) % sizes;
//		na[i][1] = (i + size) % sizes;
//		na[i][2] = (i - 1 + size) % size + (i / size) * size;
//		na[i][3] = (i + 1) % size + (i / size) * size;
//	}
//}
//double Magnet(vector<double>& v, int size)
//{
//	double m = 0;
//	for (vector<int>::size_type i = 0; i < v.size(); i++) {
//		m = m + v.at(i);
//	}
//	m = abs(m) / (v.size()); //absolute value of average spin
//	return m;
//}
//double Energy(vector<double>& v, int size, vector < vector <double> >& na)
//{
//	double e = 0;
//	for (int i = 0; i < size * size; i++) {
//		e = e - v[i] * (v[na[i][1]] + v[na[i][3]]);
//	}
//	return e;
//}
//void Cluster_1step(vector<double>& v, int size, double padd, vector < vector <double> >& na)
//{
//	int i = size * size * dis(gen);
//	vector<int> stack(1, i);
//	double oldspin = v[i];
//	double newspin = -v[i];
//	v[i] = newspin;
//	/*for (int sp = 0; sp < size * size; sp++) {
//		try {
//			if (v[na[i][0]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
//			if (v[na[i][1]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
//			if (v[na[i][2]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
//			if (v[na[i][3]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
//			i = stack.at(sp+1);
//		}
//		catch (out_of_range) { break; }
//	}*/ //Because it is so slow
//	/*int sp = 1;
//	while (sp){
//		int j = stack[--sp];
//		if (v[na[j][0]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][0]); v[na[j][0]] = newspin; sp++; } }
//		if (v[na[j][3]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][3]); v[na[j][3]] = newspin; sp++; } }
//		if (v[na[j][1]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][1]); v[na[j][1]] = newspin; sp++; } }
//		if (v[na[j][2]] == oldspin) { if (dis(gen) < padd) { stack.push_back(na[j][2]); v[na[j][2]] = newspin; sp++; } }
//	}*/ //Because it does not work
//	int sp = 0;
//	while (1) {
//		if (v[na[i][0]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][0]); v[na[i][0]] = newspin; }
//		if (v[na[i][1]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][1]); v[na[i][1]] = newspin; }
//		if (v[na[i][2]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][2]); v[na[i][2]] = newspin; }
//		if (v[na[i][3]] == oldspin && dis(gen) < padd) { stack.push_back(na[i][3]); v[na[i][3]] = newspin; }
//		sp++;
//		if (sp >= stack.size()) break;
//		i = stack.at(sp);
//	}
//
//}
//void MC_1cycle(int size, double T, double& mag, double& ene, double& mag_sus, double& sp_heat, vector < vector <double> >& na)
//{
//	int step1 = 2000, step2 = 10000;
//	int scale; 
//	if (T <= 2) scale = 1;
//	if (T > 2 && T <= 2.6) scale = int(size*size / 40);
//	if (T > 2.6 && T <= 3.2) scale = int(size * size / 10);
//	if (T > 3.2 && T <= 3.8) scale = int(size * size / 5);
//	if (T > 3.8 && T <= 4.4) scale = int(size * size / 4);
//	if (T > 4.4) scale = int(size * size / 2);
//	step1 = step1 * scale; step2 = step2 * scale;
//
//	int trash_step = int(size / 4);
//	if (T > 2.0 && T <= 2.5) trash_step = trash_step * 2;
//
//	double padd = 1 - exp(-2 / T);
//	vector<double> array(size * size, 0);
//	initialize(array, size);
//
//	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }
//
//	vector<double> magnet(step2, 0);
//	vector<double> energy(step2, 0);
//	for (int k = 0; k < step2; k++) {
//		for (int h = 0; h < trash_step; h++) {
//			Cluster_1step(array, size, padd, na);
//		}
//		magnet.at(k) = Magnet(array, size);
//		energy.at(k) = Energy(array, size, na);
//	}
//
//	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
//	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
//		Mag = Mag + magnet.at(i);
//		Mag2 = Mag2 + pow(magnet.at(i), 2);
//		Ene = Ene + energy.at(i);
//		Ene2 = Ene2 + pow(energy.at(i), 2);
//	}
//
//	mag = Mag / step2;
//	ene = Ene / step2;
//	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
//	sp_heat = (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2);
//	//return Mag, Ene, mag_sus, sp_heat;
//}
//void MC_1cycle_graphing(int size, double T, vector < vector <double> >& na)
//{
//	int step1 = 2000, step2 = 10000;
//	int scale;
//	double slope = (double(size) * size / 2.8) / 2.7;
//	if (T > 2.3) {
//		scale = slope * (T - 2.3);
//		if (scale == 0) scale = 1;
//	}
//	step1 = step1 * scale; step2 = step2 * scale;
//
//	int trash_step = int(size / 4);
//	if (T > 2.0 && T <= 2.5) trash_step = trash_step * 2;
//
//	double padd = 1 - exp(-2 / T);
//	vector<double> array(size * size, 0);
//
//	initialize(array, size);
//	std::cout << "Initial state: " << endl;
//	color(array, size);
//	std::cout << endl;
//
//	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }
//
//	vector<double> magnet(step2, 0);
//	vector<double> energy(step2, 0);
//	for (int k = 0; k < step2; k++) {
//		for (int h = 0; h < trash_step; h++) {
//			Cluster_1step(array, size, padd, na);
//		}
//		magnet.at(k) = Magnet(array, size);
//		energy.at(k) = Energy(array, size, na);
//	}
//
//	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
//	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
//		Mag = Mag + magnet.at(i);
//		Mag2 = Mag2 + pow(magnet.at(i), 2);
//		Ene = Ene + energy.at(i);
//		Ene2 = Ene2 + pow(energy.at(i), 2);
//	}
//
//	std::cout << "Final state with temperature " << T << " :" << endl;
//	color(array, size);
//	std::cout << endl;
//	std::cout << "Magnetization: " << Mag / step2 << endl;
//	std::cout << "Energy (H): " << Ene / step2 << endl;
//	std::cout << "Magnetic susceptibility: " << pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T << endl;
//	std::cout << "Specific heat: " << (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2) << endl;
//}
//
//int main()
//{
//	int size;
//	double temperature;
//
//	//cout << "Size: "; cin >> size;
//	cout << "Temperature: "; cin >> temperature;
//
//	clock_t start = clock();
//	vector < vector <double> > near(size * size, vector<double>(4, 0));
//	neighbor(near, size);
//
//	double Mag = 0, Ene = 0, mag_sus = 0, sp_heat = 0;
//	/*ofstream File;
//	File.open("ising.txt");
//	for (int k = 1; k < 51; k++) {
//		for (int h = 0; h < 10; h++) {
//			MC_1cycle(size, 0.1 * k, Mag, Ene, mag_sus, sp_heat, near);
//			cout << 0.1 * k << " " << Mag << " " << Ene << " " << mag_sus << " " << sp_heat << " " << endl;
//		}
//		cout << 0.1 * k << " finished" << endl;
//	}
//	File.close();*/
//
//	MC_1cycle_graphing(size, temperature, near);
//
//	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
//	return 0;
//}