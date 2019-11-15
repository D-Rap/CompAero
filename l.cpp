// CompAero.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
//#define _USE_MATH_DEFINES 
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;


// gradient function
double delta(double x2, double x1) {
	double grad;
	grad = (x2 - x1);
	return grad;
}

double mod(double x, double z) {
	double length;
	length = pow((pow(x, 2) + pow(z, 2)), 0.5);
	return length;
}

// smoothing function
double phi0func(double zeta, double alpha) {
	double phi0;
	phi0 = 1 - ((exp(alpha * zeta) - 1 - (alpha * zeta)) / (exp(alpha) - 1 - alpha));
	return phi0;
}
double phi1func(double zeta, double alpha) {
	double phi1;
	phi1 = zeta - ((exp(alpha * zeta) - 1 - (alpha * zeta)) / (exp(alpha) - 1 - alpha));
	return phi1;
}
double phi2func(double zeta, double alpha) {
	double phi2;
	phi2 = ((exp(alpha * zeta) - 1 - (alpha * zeta)) / (exp(alpha) - 1 - alpha));
	return phi2;
}
double lap(double p2, double p1, double p0) {
	double lp;
	lp = (0.5*(p2-p1))+(0.5*(p0-p1));
	return lp;
}

int main() {
	ofstream myfile;
	myfile.open("coordinates.plt");
	double rx, rz_in, rz, theta, delta_theta, rin, r, rloop, i, j, jmax, imax, AR, ch;
	int imax_in, imax_quarter, ii; // iii;
	rx = 0.5; // radius in x direction
	cout << "Chord length = 1\n";

	while (true) {
		cout << "Choose thickness as a percentage of chord (e.g. 1%) : ";
		cin >> rz_in;
		if (rz_in != 0 && rz_in > 0) {
			rz_in = rz_in;
			break;
		}
		else if (rz_in == 0) {
			cout << "Please choose a non-zero percentage e.g. 20% :" << endl;
		}
		else if (rz_in < 0) {
			cout << "Thickness percentage must be positive e.g. 20% :" << endl;
		}
	}
	
	while (true) {
		cout << "Choose far field distance e.g. 10 (Nfar as a multiplication of chord) :";
		cin >> rin;
		if (rin != 0 && rin > 0) {
			rin = rin;
			break;
		}
		else if (rin == 0) {
			cout << "Please choose a non-zero far field distance e.g. 10 :" << endl;
		}
		else if (rz_in < 0) {
			cout << "Far field distance must be positive e.g. 10 :" << endl;
		}
	}
	

	r = rin * 2 * rx;
	rloop = rin * rx;
	rz = (rz_in / 100) * (rx);

	// ask for imax
	while (true) {
		cout << "Choose imax (imax is the number of grid in the clockwise direction. Must be divisible by 8 : ";
		cin >> imax_in;
		if (imax_in % 8 == 0 && imax_in != 0) {
			imax = imax_in;
			break;
		}
		else if (imax_in == 0) {
			cout << "Please choose a non-zero, positive imax which is an even number :" << endl;
		}
		else if (imax_in < 0) {
			cout << "imax must be a positive integer. Please choose again. :"<<endl;
		}
		else {
			cout << "imax in not divisible by 8. Please choose a different imax. :" << endl;
		}
	}

	double PI;
	PI = 3.14159265358979323846;

	// ask for jmax
	while (true) {
		cout << "Choose jmax. (jmax is the number of section from the inner boundary surface to the far field distance) : ";
		cin >> jmax;
		if (jmax != 0 && jmax > 0 && jmax <= imax_in/3) {
			jmax = jmax;
			break;
		}
		else if (jmax == 0) {
			cout << "Please choose a non-zero, positive jmax :" << endl;
		}
		else if (jmax < 0) {
			cout << "jmax must be a positive integer. Please choose again. :" << endl;
		}
		else if (jmax > imax_in / 3) {
			cout << "jmax is recommended to be less than a third of imax. (Less than " << imax_in/3 << ") : " << endl;
		}
	}
	

	// ask for Cell height
	while (true) {
		cout << "Choose first cell height at trailing edge: ";
		cin >> ch;
		if (ch < (rin / jmax) && ch != 0) {
			ch = ch;
			break;
		}
		else if (ch < 0) {
			cout << "Cell height must be a positive number. Please choose again : " << endl;
		}
		else if (ch == 0) {
			cout << "Cell height must be a positive non zero number. Please choose again : " << endl;
		}
		else {
			cout << "Cell height is too big. Please choose a smaller number (less than " << rin/jmax << ") : " << endl;
		}
	}

	// ask for Aspect Ratio
	cout << "Choose AR: ";
	cin >> AR;

	// ELLIPSE BODY
	vector<double> x, z;
	
	// Calculate angle increment of the ellipse
	delta_theta = (2 * 3.14159265358979323846) / imax;

	theta = 0;

	for (i = 0; i <= imax_in; i++) {
		x.push_back(rx * cos(theta));
		z.push_back(rz * sin(theta));
		theta = theta - delta_theta;
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << endl;
	}
	
	// SQAURE BOX
	// find the circumference of the square box and the i increment
	double circumf, sqr_increment;// me;
	circumf = 8 * r;
	sqr_increment = circumf / imax_in;

	// check delta i for clockwise steps
	// cout << imax_in << ": i increment is: " << sqr_increment << '\n';

	vector<double> sqrx, sqrz;
	for (i = 0; i <= (imax / 8); i++) {
		sqrx.push_back(r);
		sqrz.push_back(0 - sqr_increment * i);
		// cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 1 \n";
	}
	for (i = (imax / 8) + 1; i <= ((imax / 8) * 3); i++) {
		sqrz.push_back(-r);
		sqrx.push_back((2 * r) - (sqr_increment * i));
		// out << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 2 \n";
	}
	for (i = ((imax / 8) * 3) + 1; i <= ((imax / 8) * 5); i++) {
		sqrx.push_back(-r);
		sqrz.push_back((sqr_increment * i) - 4 * r);
		// cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 3 \n";
	}
	for (i = ((imax / 8) * 5) + 1; i <= ((imax / 8) * 7); i++) {
		sqrz.push_back(r);
		sqrx.push_back((sqr_increment * i) - 6 * r);
		// cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 4 \n";
	}
	for (i = ((imax / 8) * 7) + 1; i <= imax; i++) {
		sqrx.push_back(r);
		sqrz.push_back(8 * r - (sqr_increment * i));
		// cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 5 \n";
	}

	/*for (i = 0; i <= imax_in; i++) {
		myfile << sqrx[i] << " " << 0 << " " << sqrz[i] << endl;
	}*/


// Find normal
double normalgrad;
vector<double> deltax, deltaz, normalx, normalz;

for (i = 0; i <= imax; i++) {
	if (i == 0) {
		deltax.push_back(0);
		deltaz.push_back(-1);
		normalx.push_back(1);
		normalz.push_back(0);
	}
	else if (i > 0 && i <= imax/8 ) {
		deltax.push_back(delta(x[i + 1], x[i - 1]));
		deltaz.push_back(delta(z[i + 1], z[i - 1]));
		normalgrad = -(1 / (deltaz[i] / deltax[i]));
		normalx.push_back(1 / mod(1, normalgrad));
		normalz.push_back((normalgrad / mod(1, normalgrad)));
	}
	else if (i == imax) {
		deltax.push_back(0);
		deltaz.push_back(-1);
		normalx.push_back(1);
		normalz.push_back(0);
	}
	else if (i == imax / 2) {
		deltax.push_back(0);
		deltaz.push_back(1);
		normalx.push_back(-1);
		normalz.push_back(0);
	}
	else if (i == imax / 4) {
		deltax.push_back(-1);
		deltaz.push_back(0);
		normalx.push_back(0);
		normalz.push_back(-1);
	}
	else if (i == 3 * imax / 4) {
		deltax.push_back(1);
		deltaz.push_back(0);
		normalx.push_back(0);
		normalz.push_back(1);
	}
	
	else if (i > 0 && i < imax/4) {
		deltax.push_back(delta(x[i + 1], x[i - 1]));
		deltaz.push_back(delta(z[i + 1], z[i - 1]));
		// normalgrad = (-(1 / (deltaz[i] / deltax[i])))/(mod(x[i],z[i]));
		normalgrad = -(1 / (deltaz[i] / deltax[i]));
		normalx.push_back(-1 / mod(1, normalgrad));
		normalz.push_back((normalgrad/mod(1,normalgrad)));
	}
	else if (i > imax / 4 && i < imax / 2) {
		deltax.push_back(delta(x[i + 1], x[i - 1]));
		deltaz.push_back(delta(z[i + 1], z[i - 1]));
		// normalgrad = (-(1 / (deltaz[i] / deltax[i])))/(mod(x[i],z[i]));
		normalgrad = -(1 / (deltaz[i] / deltax[i]));
		normalx.push_back(-1 / mod(1, normalgrad));
		normalz.push_back(-(normalgrad / mod(1, normalgrad)));
	}
	else if (i > imax/2 && i < 3 * imax / 4) {
		deltax.push_back(delta(x[i + 1], x[i - 1]));
		deltaz.push_back(delta(z[i + 1], z[i - 1]));
		// normalgrad = (-(1 / (deltaz[i] / deltax[i])))/(mod(x[i],z[i]));
		normalgrad = -(1 / (deltaz[i] / deltax[i]));
		normalx.push_back(-1 / mod(1, normalgrad));
		normalz.push_back(-(normalgrad / mod(1, normalgrad)));
	}
	
	else {
		deltax.push_back(delta(x[i + 1], x[i - 1]));
		deltaz.push_back(delta(z[i + 1], z[i - 1]));
		normalgrad = -(1 / (deltaz[i] / deltax[i]));
		normalx.push_back(1/mod(1,normalgrad));
		normalz.push_back(normalgrad/mod(1,normalgrad));
		
	}
}

	// vector<double> deltaxj, deltazj; 
	/*
	for (i = 0; i <= imax_in; i++) {
	// Going outwards

	deltaxj.push_back(-deltaz[i]);
	deltazj.push_back(-deltax[i]);
	
	vector<double> xj, zj;

	for (j = 0; j <= jmax - 1; j++) {
		xj.push_back(x[i] + ((j * deltaxj[i])*r));
		zj.push_back(z[i] + (((j * deltazj[i])*r)));
		cout << "outwards " << xj[j] << ", " << zj[j] << endl;
		myfile << xj[j] << ", " << zj[j] << endl;
	}

	xj.clear();
	zj.clear();
	}
	*/
	
	/*
	for (i = 0; i <= imax_in; i++) {
		// Going outwards
		deltaxj.push_back(sqrx[i] - x[i]);
		deltazj.push_back(sqrz[i] - z[i]);

		vector<double> xj, zj;

		for (j = 0; j <= jmax - 1; j++) {
			xj.push_back(x[i] + (pow((j / (jmax - 1)), 2.0) * deltaxj[i]));
			zj.push_back(z[i] + (pow((j / (jmax - 1)), 2.0) * deltazj[i]));
			cout << "outwards " << xj[j] << ", " << zj[j] << endl;
			myfile << xj[j] << " " << 0 << " " << zj[j] << endl;
		}

		xj.clear();
		zj.clear();
	}
	*/
	// clockwise that works for simplest case
	/*
	vector<double> deltaxi, deltazi;
	
	for (j = 0; j <= jmax - 1 ; j++) {
		// Going clockwise
		vector<double> xj, zj;
		double m, mm;
		m = j / (jmax - 1);
		mm = pow(m, 2.0);

		for (i = 0; i <= imax_in; i++) {
			deltaxi.push_back(mm*(sqrx[i] - x[i]));
			deltazi.push_back(mm*(sqrz[i] - z[i]));

			xj.push_back(x[i] + (deltaxi[i]));
			zj.push_back(z[i] + (deltazi[i]));
			cout << xj[i] << ", " << zj[i] << endl;
			myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
		}

		deltaxi.clear();
		deltazi.clear();
		xj.clear();
		zj.clear();
		
	}
	*/

	// TFI
	vector<double> xj, zj, lapx, lapz;
	double zet, alph, k, lambda;
	for (j = 0; j <= jmax ; j++){
		zet = (j) / (jmax-1);
		alph = 4;
		for (i = 0; i <= imax; i++) {
			if (j == 0) {
				xj.push_back(x[i]);
				zj.push_back(z[i]);
				// myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
			}
			else {
				xj.push_back(phi0func(zet, alph) * x[i] + phi1func(zet, alph) * normalx[i] + phi2func(zet, alph) * sqrx[i]);
				zj.push_back(phi0func(zet, alph) * z[i] + phi1func(zet, alph) * normalz[i] + phi2func(zet, alph) * sqrz[i]);
				// myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
			}
		}
		
		if (j == 0) {
			for (i = 0; i <= imax; i++) {
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				// cout << "jzeroth: " << xj[i] << " " << 0 << " " << zj[i] << endl;
			}
		}
		else if (j < jmax) {
						
			// Laplace Smoothing
			for (k = 0; k <= 2; k++) {
				lambda = 0.5;
				for (i = 0; i <= imax; i++) {
					if (i == 0) {
						lapx.push_back(xj[i] + lambda * lap(xj[1], xj[0], xj[imax-1]));
						lapz.push_back(zj[i] + lambda * lap(zj[1], zj[0], zj[imax-1]));
						replace(xj.begin(), xj.end(), xj[i], lapx[i]);
						replace(zj.begin(), zj.end(), zj[i], lapz[i]);
						// cout << i << ": "<< xj[i] << " " << 0 << " " << zj[i] << endl;
					}
					else if (i == imax) {
						//lapx.push_back(xj[i] + lambda * lap(xj[1], xj[imax], xj[imax-1]));
						//lapz.push_back(zj[i] + lambda * lap(zj[1], zj[imax], zj[imax-1]));
						lapx.push_back(xj[0]);
						lapz.push_back(zj[0]);
						replace(xj.begin(), xj.end(), xj[i], lapx[i]);
						replace(zj.begin(), zj.end(), zj[i], lapz[i]);
						// cout << i << ": " << xj[i] << " " << 0 << " " << zj[i] << endl;
					}
					else if (i < imax) {
						// lapx.push_back(lap(xj[i] + lambda * xj[i + 1], xj[i], xj[i - 1]));
						// lapz.push_back(lap(zj[i] + lambda * zj[i + 1], zj[i], zj[i - 1]));
						lapx.push_back(xj[i] + lambda * lap(xj[i + 1], xj[i], xj[i - 1]));
						lapz.push_back(zj[i] + lambda * lap(zj[i + 1], zj[i], zj[i - 1]));
						replace(xj.begin(), xj.end(), xj[i], lapx[i]);
						replace(zj.begin(), zj.end(), zj[i], lapz[i]);
						// cout << xj[i] << " " << 0 << " " << zj[i] << endl;
					}

				}
				//replace(xj.begin(), xj.end(), xj[imax], xj[0]);
				//replace(zj.begin(), zj.end(), zj[imax], zj[0]);
				lapx.clear();
				lapz.clear();
				
			}
			for (i = 0; i <= imax; i++) {
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				//cout << "after loop " << i << ", " << j << ": " << xj[i] << " " << 0 << " " << zj[i] << endl;
			}
		}
		else {
			for (i = 0; i <= imax; i++) {
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				//cout << "jmax: " << xj[i] << " " << 0 << " " << zj[i] << endl;
			}
		}
		
		xj.clear();
		zj.clear();
	}

	// GOING OUTWARDS but it creates star

	vector<double> xi,zi;
	for (i = 0; i <= imax; i++) {
		for (j = 0; j <= jmax; j++) {
			zet = (j) / (jmax-1);
			alph = 4;
			if (j == 0) {
				xi.push_back(x[i]);
				zi.push_back(z[i]);
				myfile << xi[j] << " " << 0 << " " << zi[j] << endl;
			}
			else {
				xi.push_back(phi0func(zet, alph) * x[i] + phi1func(zet, alph) * normalx[i] + phi2func(zet, alph) * sqrx[i]);
				zi.push_back(phi0func(zet, alph) * z[i] + phi1func(zet, alph) * normalz[i] + phi2func(zet, alph) * sqrz[i]);
				myfile << xi[j] << " " << 0 << " " << zi[j] << endl;
			}
			//cout << i << "out: " << xi[j] << " " << 0 << " " << zi[j] << endl;
		}
		xi.clear();
		zi.clear();
	}
			
	myfile.close();
	return 0;
}



