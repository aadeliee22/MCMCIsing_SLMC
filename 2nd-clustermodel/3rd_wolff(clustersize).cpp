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
//		if (i % size == 0) cout << endl;
//		if (v[i] == 1) cout << "* ";
//		if (v[i] == -1) cout << ". ";
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
//void Cluster_1step(vector<double>& v, int size, double padd, vector < vector <double> >& na, double& clustersize)
//{
//	int i = size * size * dis(gen);
//	vector<int> stack(1, i);
//	double oldspin = v[i];
//	double newspin = -v[i];
//	v[i] = newspin;
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
//	clustersize = stack.size();
//}
////void MC_1cycle(int size, double T, double& mag, double& ene, double& mag_sus, double& sp_heat, vector < vector <double> >& na)
////{
////	int step1 = 2000, step2 = 10000;
////	int trash_step = 5;
////	if (T > 2.0 && T < 2.5) trash_step = 10;
////
////	double padd = 1 - exp(-2 / T);
////	vector<double> array(size * size, 0);
////	initialize(array, size);
////
////	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }
////
////	vector<double> magnet(step2, 0);
////	vector<double> energy(step2, 0);
////	for (int k = 0; k < step2; k++) {
////		for (int h = 0; h < trash_step; h++) {
////			Cluster_1step(array, size, padd, na);
////		}
////		magnet.at(k) = Magnet(array, size);
////		energy.at(k) = Energy(array, size, na);
////	}
////
////	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
////	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
////		Mag = Mag + magnet.at(i);
////		Mag2 = Mag2 + pow(magnet.at(i), 2);
////		Ene = Ene + energy.at(i);
////		Ene2 = Ene2 + pow(energy.at(i), 2);
////	}
////
////	mag = Mag / step2;
////	ene = Ene / step2;
////	mag_sus = pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T;
////	sp_heat = (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2);
////	//return Mag, Ene, mag_sus, sp_heat;
////}
////void MC_1cycle_graphing(int size, double T, vector < vector <double> >& na)
////{
////	int step1 = 2000, step2 = 10000;
////	int trash_step = 2;
////	if (T > 2.0 && T < 2.5) trash_step = 4;
////	double padd = 1 - exp(-2 / T);
////	vector<double> array(size * size, 0);
////
////	initialize(array, size);
////	cout << "Initial state: " << endl;
////	color(array, size);
////	cout << endl;
////
////	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na); }
////
////	vector<double> magnet(step2, 0);
////	vector<double> energy(step2, 0);
////	for (int k = 0; k < step2; k++) {
////		for (int h = 0; h < trash_step; h++) {
////			Cluster_1step(array, size, padd, na);
////		}
////		magnet.at(k) = Magnet(array, size);
////		energy.at(k) = Energy(array, size, na);
////	}
////
////	double Mag = 0, Mag2 = 0, Ene = 0, Ene2 = 0;
////	for (vector<int>::size_type i = 0; i < magnet.size(); i++) {
////		Mag = Mag + magnet.at(i);
////		Mag2 = Mag2 + pow(magnet.at(i), 2);
////		Ene = Ene + energy.at(i);
////		Ene2 = Ene2 + pow(energy.at(i), 2);
////	}
////
////	cout << "Final state with temperature " << T << " :" << endl;
////	color(array, size);
////	cout << endl;
////	cout << "Magnetization: " << Mag / step2 << endl;
////	cout << "Energy (H): " << Ene / step2 << endl;
////	cout << "Magnetic susceptibility: " << pow(size, 2) * (Mag2 / step2 - pow(Mag / step2, 2)) / T << endl;
////	cout << "Specific heat: " << (Ene2 / step2 - pow(Ene / step2, 2)) / pow(size * T, 2) << endl;
////}
//
//void MC_1cycle_clustersize(int size, double T, vector < vector <double> >& na, double& cl_size)
//{
//	int step1 = 5000, step2=5000;
//	double padd = 1 - exp(-2 / T);
//	double clustersize;
//	vector<double> array(size * size, 0);
//	initialize(array, size);
//	for (int k = 0; k < step1; k++) { Cluster_1step(array, size, padd, na, clustersize); }
//	for (int k = 0; k < step2; k++) {
//		Cluster_1step(array, size, padd, na, clustersize);
//		cl_size = clustersize + cl_size;
//	}
//	cl_size = cl_size / step2;
//	//cout << "Final state with temperature " << T << " :" << endl;
//	//color(array, size);
//}
//
//int main()
//{
//	int size;
//	//double temperature;
//	
//
//	cout << "Size: "; cin >> size;
//	//cout << "Temperature: "; cin >> temperature;
//
//	clock_t start = clock();
//	vector < vector <double> > near(size * size, vector<double>(4, 0));
//	neighbor(near, size);
//
//	ofstream File;
//	File.open("clustersize.txt");
//	for (int k = 1; k < 51; k++) {
//		double cl_size = 0;
//		MC_1cycle_clustersize(size, 0.1*k, near, cl_size);
//		File << 0.1 * k << " " << cl_size << endl;
//	}
//	File.close();
//
//	cout << endl << "total time: " << (double(clock()) - double(start)) / CLOCKS_PER_SEC << " sec" << endl;
//	return 0;
//}