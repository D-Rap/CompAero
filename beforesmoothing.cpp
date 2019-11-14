// CompAero.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
//#define _USE_MATH_DEFINES 
#include <cmath>
#include <vector>
#include <numeric>

using namespace std;

int main() {
	ofstream myfile;
	myfile.open("coordinates.plt");
	double rx, ry_in, ry, theta, delta_theta, rin, r, rloop, i, j, jmax, imax, AR, ch;// delta;
	// int i, j, jmax, imax, imax_in;
	int imax_in, imax_quarter, ii; // iii;
	rx = 0.5;
	cout << "Choose thickness as a percentage of chord (e.g. 1%) : ";
	cin >> ry_in;
	cout << "Choose Square Radius (Nfar as a multiplication of chord) :";
	cin >> rin;
	r = rin * 2 * rx;
	rloop = rin * rx;
	ry = (ry_in / 100) * (rx);

	// ask for imax
	while (true) {
		cout << "Choose imax (imax is the number of grid clockwise. Must be divisible by 8: ";
		cin >> imax_in;
		if (imax_in % 8 == 0 && imax_in != 0) {
			imax = imax_in;
			imax_quarter = imax_in/4;
			break;
		}
		else if (imax_in == 0) {
			cout << "Please choose a non-zero, positive imax which is an even number" << endl;
		}
		else {
			cout << "imax in not divisible by 8. Please choose a different imax." << endl;
		}
	}

	// Calculate angle increment if the ellipse
	double PI;
	PI = 3.14159265358979323846;
	// delta_theta = (2 * PI); // / imax;
	vector<double> section;
	int sec;
	for (i = imax_quarter; i >= 1; i--) {
		section.push_back(i);
	}
	sec = accumulate(section.begin(), section.end(), 0);
	cout << "numbers of sections: " << sec << endl;

	// delta_theta = (PI / 2) / delta;

	delta_theta = (PI / 2) / sec;



	// ask for jmax
	cout << "Choose jmax: ";
	cin >> jmax;

	// ask for Cell height
	while (true) {
		cout << "Choose first cell height at trailing edge: ";
		cin >> ch;
		if (ch < (rin / jmax)) {
			ch = ch;
			break;
		}
		else {
			cout << "Cell height is too big. Please choose a smaller number: " << endl;
		}
	}

	// ask for Aspect Ratio
	cout << "Choose AR: ";
	cin >> AR;

	// ELLIPSE BODY
	vector<double> x, z;

	theta = 0;

	// Calculate angle increment if the ellipse
	delta_theta = (2 * 3.14159265358979323846) / imax;

	// vector<double> x, z;

	theta = 0;

	for (i = 0; i <= imax_in; i++) {
		x.push_back(rx * cos(theta));
		z.push_back(ry * sin(theta));
		theta = theta - delta_theta;
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << ", \n";
	}
	/*
	double dtheta;
	vector<double> theta_increment, theta_store;
	
	for (i = 0; i <= (imax/4); i++) {
		x.push_back(rx * cos(theta));
		z.push_back(ry * sin(theta));
		theta_store.push_back(theta);
		dtheta = (i + 1) * delta_theta;
		theta = theta - dtheta;
		theta_increment.push_back(dtheta);
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta_store[i] << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << ", \n";
	}
	ii = 1;
	theta = - (PI / 2);
	for (i = (imax / 4) + 1; i <= (imax / 2); i++) {
		dtheta = theta_increment[i - ii - 1];
		theta = theta - dtheta;
		x.push_back(rx * cos(theta));
		z.push_back(ry * sin(theta));
		theta_store.push_back(theta);
		ii = ii + 2;
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta_store[i] << "ii: " << ii << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << ", \n";
	}
	for (i = (imax / 2) + 1; i <= ((3 * imax) / 4); i++) {
		dtheta = (i - (imax / 2)) * delta_theta;
		theta = theta - dtheta;
		x.push_back(rx * cos(theta));
		z.push_back(ry * sin(theta));
		theta_store.push_back(theta);
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta_store[i] << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << ", \n";
	}
	ii = 1;
	for (i = ((3 * imax) / 4) + 1; i <= imax; i++) {
		dtheta = ((i - ((imax / 2))) - ii) * delta_theta;
		theta = theta - dtheta;
		x.push_back(rx * cos(theta));
		z.push_back(ry * sin(theta));
		theta_store.push_back(theta);
		ii = ii + 2;
		cout << i << ": " << x[i] << ", " << z[i] << ", " << "theta: " << theta_store[i] << " ii: " << ii << ", \n";
		myfile << x[i] << " " << 0 << " " << z[i] << ", \n";
	}
	*/
	// SQAURE BOX
	// find the circumference of the square box and the i increment
	double circumf, sqr_increment, me;
	circumf = 8 * r;
	sqr_increment = circumf / imax_in;
	// check delta i for clockwise steps
	cout << imax_in << ": i increment is: " << sqr_increment << '\n';

	vector<double> sqrx, sqrz;
	for (i = 0; i <= (imax / 8); i++) {
		sqrx.push_back(r);
		sqrz.push_back(0 - sqr_increment * i);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 1 \n";
	}
	for (i = (imax / 8) + 1; i <= ((imax / 8) * 3); i++) {
		sqrz.push_back(-r);
		sqrx.push_back((2 * r) - (sqr_increment * i));
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 2 \n";
	}
	for (i = ((imax / 8) * 3) + 1; i <= ((imax / 8) * 5); i++) {
		sqrx.push_back(-r);
		sqrz.push_back((sqr_increment * i) - 4 * r);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 3 \n";
	}
	for (i = ((imax / 8) * 5) + 1; i <= ((imax / 8) * 7); i++) {
		sqrz.push_back(r);
		sqrx.push_back((sqr_increment * i) - 6 * r);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 4 \n";
	}
	for (i = ((imax / 8) * 7) + 1; i <= imax; i++) {
		sqrx.push_back(r);
		sqrz.push_back(8 * r - (sqr_increment * i));
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 5 \n";
	}
	
	// cout << "zeroth element" << sqrx[0] << ", " << sqry[0] << ", \n";
	for (i = 0; i <= imax_in; i++) {
		myfile << sqrx[i] << " " << 0 << " " << sqrz[i] << endl;
	}
	
	/*
	for (i = 0; i <= (imax / 8); i++) {
		me = z[i] / x[i];
		sqrx.push_back(r);
		sqrz.push_back(me*sqrx[i]);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 1 \n";
	}
	for (i = (imax / 8) + 1; i <= ((imax / 8) * 3); i++) {
		me = z[i] / x[i];
		sqrz.push_back(-r);
		sqrx.push_back(me*sqrz[i]);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 2 \n";
	}
	for (i = ((imax / 8) * 3) + 1; i <= ((imax / 8) * 5); i++) {
		me = z[i] / x[i];
		sqrx.push_back(-r);
		sqrz.push_back(me*sqrx[i]);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 3 \n";
	}
	for (i = ((imax / 8) * 5) + 1; i <= ((imax / 8) * 7); i++) {
		me = z[i] / x[i];
		sqrz.push_back(r);
		sqrx.push_back(me*sqrz[i]);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 4 \n";
	}
	for (i = ((imax / 8) * 7) + 1; i <= imax; i++) {
		me = z[i] / x[i];
		sqrx.push_back(r);
		sqrz.push_back(me*sqrx[i]);
		cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 5 \n";
	}
	*/

	// AR consideration
	double deltax1, deltax2, deltaz1, deltaz2;
	vector<double> deltax, deltaz;
	for (i = 0; i <= imax; i++) {
		if (i == 0) {
			deltax1 = x[i] - x[imax-1];
			deltax2 = x[i + 1] - x[i];
			deltax.push_back(AR* (0.5 * (deltax2 + deltax1)));

			deltaz1 = z[i] - z[imax-1];
			deltaz2 = z[i + 1] - z[i];
			deltaz.push_back(AR* (0.5 * (deltaz2 + deltaz1)));

			cout << i << ": " << deltax[i] << ", " << deltaz[i] << "\n";
		}
		else if (i <= imax - 1) {
			deltax1 = x[i] - x[i - 1];
			deltax2 = x[i + 1] - x[i];
			deltax.push_back(AR * (0.5 * (deltax2 + deltax1)));

			deltaz1 = z[i] - z[i - 1];
			deltaz2 = z[i + 1] - z[i];
			deltaz.push_back(AR * (0.5 * (deltaz2 + deltaz1)));

			cout << i << ": " << deltax[i] << ", " << deltaz[i] << "\n";
		}
		else if (i == imax) {
			deltax1 = x[i] - x[i - 1];
			deltax2 = x[1] - x[i];
			deltax.push_back(AR * (0.5 * (deltax2 + deltax1)));

			deltaz1 = z[i] - z[i-1];
			deltaz2 = z[1] - z[i];
			deltaz.push_back(AR * (0.5 * (deltaz2 + deltaz1)));

			cout << i << ": " << deltax[i] << ", " << deltaz[i] << "\n";
		}
	}

	// Find normal
	vector<double> normalgrad;
	for (i = 0; i <= imax; i++) {
		if (deltax[i] != 0) {
			if (deltaz[i] != 0) {
				// grad.push_back(deltaz[i] / deltax[i]);
				normalgrad.push_back(-((deltaz[i]) / (deltax[i])));
				cout << i << ": " << normalgrad[i] << endl;
			}
			else {
				normalgrad.push_back(-deltax[i]);
				cout << i << ": " << normalgrad[i] << endl;
			}
		}
		else {
			normalgrad.push_back(-deltaz[i]);
			cout << i << ": " << normalgrad[i] << endl;
		}
	}
	/*
	double sqrangle;
	// with angle boundaries
	for (i = 0; i <= imax; i++) {
		sqrangle = theta_store[i];
		if (sqrangle >= - PI / 4) {
			// me = normalgrad[i];
			sqrx.push_back(r);
			sqrz.push_back(normalgrad[i]* sqrx[i]);
			cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 1 \n";
			myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
		}
		else if (sqrangle >= -((3 * PI) / 4)) {
			// me = z[i] / x[i];
			sqrz.push_back(-r);
			sqrx.push_back((1/normalgrad[i]) * sqrz[i]);
			cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 2 \n";
			myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
		}
		else if (sqrangle >= -((5 * PI) / 4)) {
			// me = z[i] / x[i];
			sqrx.push_back(-r);
			sqrz.push_back(normalgrad[i] * sqrx[i]);
			cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 3 \n";
			myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
		}
		else if (sqrangle >= -((7*PI)/4)) {
			// me = z[i] / x[i];
			sqrz.push_back(r);
			sqrx.push_back((1/normalgrad[i]) * sqrz[i]);
			cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 4 \n";
			myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
		}
		else if (sqrangle >= -2*PI) {
			// me = z[i] / x[i];
			sqrx.push_back(r);
			sqrz.push_back(normalgrad[i] * sqrx[i]);
			cout << i << ' ' << sqrx[i] << ',' << sqrz[i] << ", 5 \n";
			myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
		}
	}
	*/
	/*
	for (i = 0; i <= imax_in; i++) {
		myfile << sqrx[i] << ", " << sqrz[i] << ", \n";
	}
	*/
	
	// myfile.close();
	//ofstream myfile1;
	//myfile1.open("square.plt");

	
	
	vector<double> deltaxj, deltazj; 
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

	// clockwise by updating radius of ellipse
	double newr;
	vector<double> xj, zj;
	for (j = 0; j <= jmax-1; j++) {
		theta = 0;
		for (i = 0; i <= imax; i++) {
			
			if (j == 0) {
				xj.push_back(x[i]);
				zj.push_back(z[i]);
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				cout << i << ", " << j << ": " << xj[i] << ", " << zj[i] << endl;
			}
			/*
			else if (j == 1) {
				// newr = rx + ch;
				// xj.push_back(newr*cos(theta));
				// zj.push_back((ry_in / 100)* newr*sin(theta));
				// theta = theta - delta_theta;
				xj.push_back(x[i] + (ch * deltaxj[i]));
				zj.push_back(z[i] + (ch * deltazj[i]));
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				cout << i << ", " << j << ": " << xj[i] << ", " << zj[i] << endl;
			}
			*/
			else {
				//newr = rx + (pow((j / (jmax - 1)), 2.0) * r);
				//xj.push_back(newr* cos(theta));
				//zj.push_back((ry_in / 100) * newr* sin(theta));
				//theta = theta - delta_theta;

				// xj.push_back(x[i] + (ch * j/rin * deltaxj[i]));
				// zj.push_back(z[i] + (ch * j /rin * deltazj[i]));
				xj.push_back(x[i] + (pow((j / (jmax - 1)), 2.0) * deltaxj[i]));
				zj.push_back(z[i] + (pow((j / (jmax - 1)), 2.0) * deltazj[i]));
				myfile << xj[i] << " " << 0 << " " << zj[i] << endl;
				cout << i << ", " << j << ": " << xj[i] << ", " << zj[i] << endl;
			}
		}
		xj.clear();
		zj.clear();
	}

	myfile.close();
	return 0;
}



