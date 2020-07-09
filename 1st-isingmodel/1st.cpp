//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <random>
//#include <vector>
//#include <time.h>
//using namespace std;
//random_device rd;
//mt19937 gen(rd());
//uniform_real_distribution<double> dis(0, 1);
//
//void initialize(vector < vector <double> >& v, int size) //initial -random- state
//{
//	for (int i = 0; i < size; i++)
//	{
//		for (int j = 0; j < size; j++)
//		{
//			v[i][j] = dis(gen) < 0.5 ? 1 : -1;
//		}
//	}
//}
//void color(vector < vector <double> >& v, int size) //graphing state
//{
//	for (int i = 0; i < size; i++)
//	{
//		for (int j = 0; j < size; j++)
//		{
//			if (v[i][j] == 1) cout << "* ";
//			else cout << ". ";
//		} cout << endl;
//	}
//}
//double Magnet(vector < vector <double> >& v, int size)
//{
//	double m = 0;
//	for (int i = 0; i < size; i++){
//		for (int j = 0; j < size; j++){
//			m = m + v[i][j];
//		}
//	}
//	m = abs(m / (double(size) * size)); //absolute value of average spin
//	return m;
//}
//double Energy(vector < vector <double> >& v, int size)
//{
//	double e = 0;
//	for (int i = 0; i < size; i++){
//		for (int j = 0; j < size; j++){
//			e = e - v[i][j] * (v[(i + 1) % size][j] + v[i][(j + 1) % size]);
//		}
//	}
//	return e;
//}
//double delU(vector < vector <double> >& v, int size, int x, int y)
//{
//	double top = v[(x + (size-1)) % size][y % size];
//	double left = v[x % size][(y + (size - 1)) % size];
//	double bottom = v[(x + 1) % size][y % size];
//	double right = v[x % size][(y + 1) % size];
//	double E = 2 * v[x][y] * (top + left + bottom + right);
//	return E;
//}
//double exp_delU(double E, double* expE)
//{
//	double result;
//	if (E == 8) result = *(expE);
//	else if (E == 4) result = *(expE + 1);
//	else result = 1;
//	return result;
//}
//void MC_1step(vector < vector <double> >& v, int size, double* expE)
//{
//	for (int i = 0; i < size; i++) {
//		int j = i % 2;
//		for (int k = 0; k<int(size / 2); k++){
//			double Ediff = delU(v, size, i, j);
//			if (Ediff <= 0) v[i][j] = -v[i][j];
//			else {
//				if (dis(gen)<=exp_delU(Ediff, expE)) v[i][j]=-v[i][j];
//			}
//			j = j + 2;
//		}
//	}
//	for (int i = 1; i < size+1; i++) {
//		int j = i % 2;
//		for (int k = 0; k<int(size / 2); k++)
//		{
//			double Ediff = delU(v, size, i-1, j);
//			if (Ediff <= 0) v[i-1][j] = -v[i-1][j];
//			else {
//				if (dis(gen) < exp_delU(Ediff, expE)) v[i-1][j] = -v[i-1][j];
//			}
//			j = j + 2;
//		}
//	}
//
//}
//void MC_1cycle(int size, double T, double& mag, double& ene, double& mag_sus, double& sp_heat)
//{
//	int step1 = 1000, step2 = 1000, step3 = 1000;
//	double expE[2] = { exp(-8 / T), exp(-4 / T) };
//	vector < vector <double> > array(size, vector<double>(size, 0));
//
//	initialize(array, size);
//
//	/*int num = int(step1 / 10);
//	double TT = T + 2;
//	double expE1[5] = { exp(-8 / TT), exp(-4 / TT)};
//	for (int k = 0; k < step1; k++){
//		MC_1step(array, size, &expE1[0]);
//		if (k % num == 0) {
//			TT = TT - double(k) / (double(num) * 5);
//			double expE1[5] = { exp(-8 / TT), exp(-4 / TT) };
//		}
//	}*/
//	for (int k = 0; k < step2; k++) { MC_1step(array, size, &expE[0]); }
//
//	vector<double> magnet(step3, 0);
//	vector<double> energy(step3, 0);
//	for (int k = 0; k < step3; k++){
//		MC_1step(array, size, &expE[0]);
//		magnet.at(k) = Magnet(array, size);
//		energy.at(k) = Energy(array, size);
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
//	mag = Mag/step3;
//	ene = Ene/step3;
//	mag_sus = pow(size, 2) * (Mag2 / step3 - pow(Mag / step3, 2)) / T;
//	sp_heat = (Ene2 / step3 - pow(Ene / step3, 2)) / pow(size * T, 2);
//	//return Mag, Ene, mag_sus, sp_heat;
//}
//void MC_1cycle_graphing(int size, double T)
//{
//	int step1 = 1000, step2 = 700, step3 = 800;
//	double expE[2] = { exp(-8 / T), exp(-4 / T) };	
//	vector < vector <double> > array(size, vector<double>(size, 0));
//
//	initialize(array, size);
//	cout << "Initial state: " << endl;
//	color(array, size);
//	cout << endl;
//
//	/*int num = int(step1 / 10);
//	double TT = T + 2;
//	double expE1[5] = { exp(-8 / TT), exp(-4 / TT) };
//	for (int k = 0; k < step1; k++){
//		MC_1step(array, size, &expE1[0]);
//		if (k % num == 0) { 
//			TT = TT - double(k) / (double(num) * 5);
//			double expE1[5] = { exp(-8 / TT), exp(-4 / TT) };
//		}
//	}*/
//	for (int k = 0; k < step2; k++) { MC_1step(array, size, &expE[0]); }
//
//	vector<double> magnet(step3, 0);
//	vector<double> energy(step3, 0);
//	for (int k = 0; k < step3; k++){
//		MC_1step(array, size, &expE[0]);
//		magnet.at(k) = Magnet(array, size);
//		energy.at(k) = Energy(array, size);
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
//	cout << "Final state with temperature " << T <<" :"<< endl;
//	color(array, size);
//	cout << endl;
//	cout << "Magnetization: " << Mag/step3 << endl;
//	cout << "Energy (H): " << Ene/step3 << endl;
//	cout << "Magnetic susceptibility: " << pow(size, 2) * (Mag2 / step3 - pow(Mag / step3, 2)) / T << endl;
//	cout << "Specific heat: " << (Ene2 / step3 - pow(Ene / step3, 2)) / pow(size * T, 2) << endl;
//}
//
//int main()
//{
//	clock_t start = clock();
//
//	double Mag = 0, Ene = 0, mag_sus = 0, sp_heat = 0;
//	ofstream File;
//	File.open("ising.txt");
//
//	for (int k = 1; k < 51; k++) {
//		for (int h = 0; h < 10; h++) {
//			MC_1cycle(10, 0.1 * k, Mag, Ene, mag_sus, sp_heat);
//			//cout << 0.1 * k << " " << Mag << " " << Ene << " " << mag_sus << " " << sp_heat << " " << endl;
//		}
//		cout << 0.1 * k << " finished" << endl;
//	}
//
//	File.close();
//	//MC_1cycle_graphing(20, 2.4);
//
//	cout << endl << "total time: " << (double(clock()) - double(start)) / 1000 << " sec" << endl;
//	return 0;
//}