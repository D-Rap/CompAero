// CompAero.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
//#define _USE_MATH_DEFINES 
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
	double rx, ry_in, ry, theta, delta_theta, rin, r, rloop, i, j, jmax, imax;
	// int i, j, jmax, imax, imax_in;
	int imax_in;
	rx = 0.5;
	cout << "Choose thickness as a percentage of chord (e.g. 1%) : ";
	cin >> ry_in;
	cout << "Choose Square Radius (Nfar as a multiplication of chord) :";
	cin >> rin;
	r = rin * 2 * rx;
	rloop = rin * rx;
	ry = (ry_in / 100) * (2 * rx);

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

	// Calculate angle increment if the ellipse
	delta_theta = (2 * 3.14159265358979323846) / imax;

	// ask for jmax
	cout << "Choose jmax: ";
	cin >> jmax;

	vector<double> x, z;

	theta = 0;

	for (i = 0; i <= imax_in; i++) {
		x.push_back(rx*cos(theta));
		z.push_back(ry*sin(theta));
		theta = theta - delta_theta;
		cout << i << ": " << x[i] << ", " << z[i] <<", " << "theta: " << theta << ", \n";
		myfile << x[i] << ", " << z[i] << ", \n";
	}

	/*myfile.close();*/

	//ofstream myfile1;
	//myfile1.open("square.plt");

	// find the circumference of the square box and the i increment
	double circumf, sqr_increment;
	circumf = 8*r;
	sqr_increment = circumf / imax_in;
	// check delta i for clockwise steps
	cout << imax_in << ": i increment is: " << sqr_increment << '\n';

	vector<double> sqrx, sqrz;
	for (i = 0; i <= (imax / 8); i++) {
		sqrx.push_back(r);
		sqrz.push_back(0 - sqr_increment * i);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 1 \n"; 
	}
	for (i = (imax / 8)+1; i <= ((imax / 8) * 3); i++) {
		sqrz.push_back(-r);
		sqrx.push_back((2*r) - (sqr_increment*i));
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 2 \n";
	}
	for (i = ((imax / 8) * 3)+1; i <= ((imax / 8) * 5); i++) {
		sqrx.push_back(-r);
		sqrz.push_back((sqr_increment * i) - 4 * r);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] <<", 3 \n";
	}
	for (i = ((imax / 8) * 5)+1; i <= ((imax / 8) * 7); i++) {
		sqrz.push_back(r);
		sqrx.push_back((sqr_increment * i) - 6*r);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] <<", 4 \n";
	}
	for (i = ((imax / 8) * 7)+1; i <= imax; i++) {
		sqrx.push_back(r);
		sqrz.push_back(8 * r - (sqr_increment * i));
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] <<", 5 \n";
	}

	// cout << "zeroth element" << sqrx[0] << ", " << sqry[0] << ", \n";
	for (i = 0; i <= imax_in; i++) {
		myfile << sqrx[i] << ", " << sqrz[i] << endl;
	}
	
	vector<double> deltax, deltaz; 

	//for (i = 0; i <= imax_in; i++) {
	//	// Going outwards
	//	deltax.push_back(sqrx[i] - x[i]);
	//	deltaz.push_back(sqrz[i] - z[i]);

	//	vector<double> xj, zj;

	//	for (j = 0; j <= jmax - 1; j++) {
	//		xj.push_back(x[i] + (pow((j /(jmax - 1)), 2.0) * deltax[i]));
	//		zj.push_back(z[i] + (pow((j / (jmax - 1)), 2.0) * deltaz[i]));
	//		cout << "outwards " << xj[j] << ", " << zj[j] << endl;
	//		myfile << xj[j] << ", " << zj[j] << endl;
	//	}

	//	xj.clear();
	//	zj.clear();
	//}

	for (i = 0; i <= imax_in; i++) {
		// Going outwards
		deltax.push_back(sqrx[i] - x[i]);
		deltaz.push_back(sqrz[i] - z[i]);

		vector<double> xj, zj;

		for (j = 0; j <= jmax - 1; j++) {
			xj.push_back(x[i] + (pow((j / (jmax - 1)), 2.0) * deltax[i]));
			zj.push_back(z[i] + (pow((j / (jmax - 1)), 2.0) * deltaz[i]));
			cout << "outwards " << xj[j] << ", " << zj[j] << endl;
			myfile << xj[j] << ", " << zj[j] << endl;
		}

		xj.clear();
		zj.clear();
	}

	deltax.clear();
	deltaz.clear();

	for (j = 0; j <= jmax - 1 ; j++) {
		// Going clockwise
		vector<double> xj, zj;
		double m, mm;
		m = j / (jmax - 1);
		mm = pow(m, 2.0);

		for (i = 0; i <= imax_in; i++) {
			deltax.push_back(mm*(sqrx[i] - x[i]));
			deltaz.push_back(mm*(sqrz[i] - z[i]));

			xj.push_back(x[i] + (deltax[i]));
			zj.push_back(z[i] + (deltaz[i]));
			cout << xj[i] << ", " << zj[i] << endl;
			myfile << xj[i] << ", " << zj[i] << endl;
		}

		deltax.clear();
		deltaz.clear();
		xj.clear();
		zj.clear();
		
	}

	myfile.close();
	return 0;
}



