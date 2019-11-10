// CompAero.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES 
#include <cmath>
#include <vector>

//int main()
//{
//    std::cout << "Hello World!\n";
//}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

using namespace std;

int main() {
	ofstream myfile;
	myfile.open("coordinates.plt");
	// float x, y, rx, ry_in, ry, theta, delta_theta, imax, jmax;
	int j, i, jmax, imax;
	int imax_in;
	// ask for imax
	while (true) {
		cout << "Choose imax (imax must be divisible by 8): ";
		cin >> imax_in;
		if (imax_in % 8 == 0) {
			imax = imax_in;
			break;
		}
		else {
			cout << "imax in not divisible by 8. Please choose a different imax." << endl;
		}
	}
	const int imax_array = imax_in + 1;

	// ask for jmax
	cout << "Choose jmax: ";
	cin >> jmax;

	double rx, ry_in, ry, theta, delta_theta, rin, r, rloop;
	// x[imax_array], y[imax_array];
	rx = 0.5;
	cout << "Choose thickness as a percentage of chord (e.g. 1%) : ";
	cin >> ry_in;
	cout << "Choose Square Radius (Nfar as a multiplication of chord) :";
	cin >> rin;
	r = rin*2*rx;
	rloop = rin * rx;
	ry = (ry_in / 100) * (2*rx);

	// Calculate angle increment if the ellipse
	delta_theta = (2 * 3.14159265358979323846) / imax;

	vector<double> x, y;

	theta = 0;

	for (j = 0; j <= imax_in; j++) {
		x.push_back(rx*cos(theta));
		y.push_back(ry*sin(theta));
		theta = theta + delta_theta;
		cout << j << ": " << x[j] << ", " << y[j] <<", " << "theta: " << theta << ", \n";
		myfile << j << ": " << x[j] << ", " << y[j] << ", \n";
	}

	myfile.close();

	ofstream myfile1;
	myfile1.open("square.plt");

	//// find the circumference of the square box and the i increment
	double circumf, sqr_increment;
	circumf = 8*r;
	sqr_increment = circumf / imax_in;
	// check delta i for clockwise steps
	cout << imax_in << "i increment is: " << sqr_increment << '\n';

	vector<double> sqrx, sqry;
	for (i = 0; i <= (imax / 8); i++) {
		sqrx.push_back(r);
		sqry.push_back(0 - sqr_increment * i);
		cout << i << ' ' << sqrx[i] << ',' << sqry[i] << ", 1 \n"; 
	}
	for (i = (imax / 8)+1; i <= ((imax / 8) * 3); i++) {
		sqry.push_back(-r);
		sqrx.push_back((2*r) - (sqr_increment*i));
		cout << i << ' ' << sqrx[i] << ',' << sqry[i] << ", 2 \n";
	}
	for (i = ((imax / 8) * 3)+1; i <= ((imax / 8) * 5); i++) {
		sqrx.push_back(-r);
		sqry.push_back((sqr_increment * i) - 4 * r);
		cout << i << ' ' << sqrx[i] << ',' << sqry[i] <<", 3 \n";
	}
	for (i = ((imax / 8) * 5)+1; i <= ((imax / 8) * 7); i++) {
		sqry.push_back(r);
		sqrx.push_back((sqr_increment * i) - 6*r);
		cout << i << ' ' << sqrx[i] << ',' << sqry[i] <<", 4 \n";
	}
	for (i = ((imax / 8) * 7)+1; i <= imax; i++) {
		sqrx.push_back(r);
		sqry.push_back(8 * r - (sqr_increment * i));
		cout << i << ' ' << sqrx[i] << ',' << sqry[i] <<", 5 \n";
	}

	//cout << "zeroth element" << sqrx[0] << ", " << sqry[0] << ", \n";
	for (i = 0; i <= imax_in; i++) {
		myfile1 << i << sqrx[i] << ", " << sqry[i] << ",\n";
	}
	myfile1.close();
	//return 0;
}



