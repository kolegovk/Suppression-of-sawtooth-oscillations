/*
CaphydropFDM (capillary hydrodynamics in droplet via Finite Difference Method)
(c) 2023 Konstantin S. Kolegov
https://orcid.org/0000-0002-9742-1308 
email: k.kolegov87@gmail.com

This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

The program calculates the capillary hydrodynamics in a liquid droplet evaporating
on a substrate and containing a suspended or dissolved substance
based on the lubrication approximation. The finite difference method is
used here. The problem of sawtooth oscillations has been solved. 

For more information see the paper:
Kolegov, K. S. (2023). Suppression of sawtooth oscillations when using a finite-difference scheme 
for mass transfer simulation via the lubrication approximation in a droplet evaporated on a substrate.
https://doi.org/10.48550/arXiv.2301.06983
(in Russian language)

Version 7
Checking of the  mass conservation law has been added.
*/
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;
using counter = unsigned long long;

//////////////////////
//Geometric parameters
//////////////////////

const double h0 = 0.1e-3;//droplet height [m]
const double h0_ = 1.0;//dimensionless value
const double R = 1.0e-3;//radius of the droplet base [m]
const double R_ = 1.0;//dimensionless value
const double hf = 0.01 * h0;//height of the liquid layer near the contact line [m]
const double hf_ = hf / h0;//dimensionless value
const double eps = h0 / R;//aspect ratio
const double theta = 2.0 * eps;//contact angle, radians

////////////
// Mesh
////////////

// Mesh resolution (N+1) x (M+1) nodes
const int N = 75;// horizontal direction
const int M = 75;// vertical direction
const int Np1 = N + 1;
const int Nm1 = N - 1;
const int Nm2 = N - 2;
const int Mp1 = M + 1;
//space step counters
int n = 0, m = 0;
// spatial steps
double dr = R / N;
double dr_ = dr / R;//dimensionless value
double dz = h0 / M;
double dz_ = dz / h0;//dimensionless value
//time step, seconds
double dt = 1.0e-6; //1.0e-4;  5.0e-6;
//vectors of spatial coordinates
double r[Np1];
double z[Mp1];

//////////////////////
// Physical parameters
///////////////////////

const double C0 = 0.05;//initial mass fraction of the substance
const double Cg = 0.7;// gelation concentration
const double Cg2 = Cg * Cg;
const double S = 1.692; //parameter in Mooney’s formula
const double K = 1.236; //parameter in Mooney’s formula
const double D = 1.0e-10;//diffusion coefficient of a dissolved or suspended substance [m2/s]
const double eta = 0.001;//viscosity of the liquid [Pa*s]
const double eta_Cg = exp(S * Cg / (1.0 - K * Cg));//dimensionless viscosity at critical concentration
const double rho = 1000.0;//liquid density [kg/m3]
const double sigma = 0.072;//surface tension [N/m]

const double uc = eta / (rho * h0);//characteristic velocity [m/s]
const double tc = R / uc;//characteristic time [s]

double tf = 450.0;//full evaporation time [s]

const double Dv = 2.4e-5;//vapor diffusion coefficient [m2/s]
const double Cv = 2.32e-2;//saturation vapor density [kg/m3]
const double H = 0.4;//relative humidity
const double J0 = (Dv * Cv * (1 - H) / R) * \
(0.27 * theta * theta + 1.3) * \
(0.6381 - 0.2239 * (theta - M_PI / 4.0) * (theta - M_PI / 4.0));//vapour flux density [kg/(s*m2)]
const double Jc = eps * rho * uc;//characteristic vapour flux density
const double J0_ = J0 / Jc;//dimensionless vapour flux density
const double width = 30.0; //deposit width adjustment for initial condition
const double kappa = 1.0;//emperical parameter for J

//////////////////////////////
// Dimensionless parameters
//////////////////////////////

const double Ca = eta * uc / (sigma * eps * eps * eps); // Capillary number
const double Pe = R * uc / D; // Peclet number

/////////////////////////////////////////
//Parameters for the numerical method
/////////////////////////////////////////
double ErrorLimit = 1.0e-17;
double maxError = 0.0;
const double Relaxation = 1.0e-3; // Relaxation parameter
const int iterNum = 101;//number of iterations

///////////////////////////////////
//Dimensionless functions
//////////////////////////////////

double h[Np1];//current time layer of h
double hTMP[Np1];//temporary value of h
double hn[Np1];//new time layer of h
double P[Np1];//capillary pressure
double u_a[Np1];//flow velocity averaged over the thickness of the liquid layer
double C[Np1];//mass fraction
double Cn[Np1];//mass fraction on next time step
double J[Np1];//evaporation flux density
double U[Mp1][Np1];// velocity vector field of fluid flow: horizontal component
double W[Mp1][Np1];// velocity vector field of fluid flow: vertical component
double eta_[Np1];//viscosity
double Ha[Np1]; //Analytical form of Heaviside function

/////////////////////////////////////
//Working with strings and files
/////////////////////////////////////
std::ofstream fout; // File access
std::string str = ""; // Message output

inline void tridiag(double a[], double b[], double c[], double d[], double x[]);
inline void SaveDataToFile(double);

int main(int argc, char* argv[]) {
	cout << "CaphydropFDM (capillary hydrodynamics in droplet via Finite Difference Method)." << endl;
	cout << "(c) 2023 Konstantin S. Kolegov" << endl;
	cout << "email: k.kolegov87@gmail.com" << endl;
	cout << endl;
	cout << "This program is free software (GNU General Public License v3.0)." << endl;
	cout << "It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY!" << endl;
	cout << endl;
	cout << "The program calculates the capillary hydrodynamics in a liquid droplet evaporating \
on a substrate and containing a suspended or dissolved substance \
based on the lubrication approximation.The finite difference method is \
used here.The problem of sawtooth oscillations has been solved." << endl;
	cout << endl;
	cout << "Kolegov, K. S. (2023). Suppression of sawtooth oscillations when using a finite-difference scheme \
for mass transfer simulation via the lubrication approximation in a droplet evaporated on a substrate. \
https://doi.org/10.48550/arXiv.2301.06983" << endl;
	cout << endl;
	cout << "Start of a numerical calculation..." << endl;
	//start timer
	clock_t t = clock();

	///////////////////////////////////////////////
	//Reading the parameters from the command line
	///////////////////////////////////////////////
	//Running the program with parameters on the command line: program.exe tf dt
	for (int i = 0; i < argc; i++) {
		if (i == 1) {
			std::string str1(argv[i]);
			tf = std::stod(str1);
		}
		if (i == 2) {
			std::string str2(argv[i]);
			dt = std::stod(str2);
		}
	}


	//Time 
	const double tf_ = tf/ tc;
	const double dt_ = dt/ tc;
	const counter NumTimeSteps = (counter)trunc(tf_ / dt_);
	cout << "Number of time steps is needed= " << NumTimeSteps << endl;
	double currentModelTime = 0.0;//time of process, seconds
	const int numberOfPeriods = 45;//number of time periods
	counter periodSize = NumTimeSteps / (counter)numberOfPeriods;//Time period for adding plot data
	int periodCount = 0; // Counter for adding data to the plot/ file
	bool finishFlag = false;//flag of the calculation end 

	//auxiliary variables
	const double halfDt_ = 0.5 * dt_;
	const double dr_Pe = Pe * dr_;
	const double dr_2 = 2.0 * dr_;
	const double dr_2Ca = dr_2 * Ca;
	const double dr_dr_2Ca = dr_2Ca * dr_;
	const double dr_6 = 6.0 * dr_;

	//Values of spatial grid nodes
	for (m = 0; m <= M; m++) z[m] = (double)m * dz_;
	for (n = 0; n <= N; n++) r[n] = (double)n * dr_;

	/////////////////////////////////////////////
	//Initial conditions
	/////////////////////////////////////////////

	for (n = 0; n <= N; n++) {
		hTMP[n] = hn[n] = h[n] = hf_ + h0_ * (1 - r[n] * r[n] / (R_ * R_));//droplet shape
		u_a[n] = 0.0;//the flow is missing at the initial moment of time
		C[n] = Cn[n] = Cg * (2.0 - C0 + 2.0 * (C0 - 1.0) / (1.0 + exp(width * (r[n] - 1.0))));//mass fraction
		J[n] = J0_* (1.0 - C[n] * C[n] / Cg2) / (kappa + h[n]);//evaporation flux density
		eta_[n] = exp(S * C[n] / (1.0 - K * C[n]));//viscosity
		Ha[n] = 1.0 / (1.0 + exp(-4000.0 * (Cg - Cn[n] - 0.005))); //Analytical form of Heaviside function
		for (m = 0; m <= M; m++) {
			//The flow is missing at the initial moment of time.
			U[m][n] = 0.0;
			W[m][n] = 0.0;
		}
	}

	for (n = 1; n < N; n++) //central finite difference scheme
		P[n] = ((r[n] + r[n + 1]) * (hn[n] - hn[n + 1]) / dr_ \
			- (r[n] + r[n - 1]) * (hn[n - 1] - hn[n]) / dr_) / (r[n] * dr_2Ca);

	//boundary condition dP/dr=0
	P[0] = (4.0 * P[1] - P[2]) / 3.0;//second-order approximation
	P[N] = 4.0 * P[Nm1] - 3.0 * P[Nm2];//second-order approximation

	//Checking of the  mass conservation law.
	double mass = 0.0;
	for (int n = 0; n <= N; n++)
		mass += hn[n] * Cn[n];

	cout << "Parameters:" << endl;
	cout << "Ca=" << Ca << endl;
	cout << "Pe=" << Pe << endl;
	cout << "dt_=" << dt_ << endl;
	cout << "dr_=" << dr_ << endl;
	cout << "tc=" << tc << " s" << endl;
	cout << "tf=" << tf << " s" << endl;
	cout << "J0=" << J0 * 1.0e3 << " g/(s*m^2)" << endl;
	cout << "N=" << N << endl;
	cout << "M=" << M << endl;
	cout << "The mass of the substance at the beginning: " << mass << endl;

	///////////////////////////////////////////////////////
	//Write initial data to file
	///////////////////////////////////////////////////////

	SaveDataToFile(currentModelTime);

	///////////////////////////////////////////////
	//Time periods
	///////////////////////////////////////////////
	counter k = 1;//time step counter
	while (periodCount < numberOfPeriods) {
		//////////////
		//Time steps
		/////////////
		for (k = 1; k <= periodSize; k++) {

			currentModelTime = (double)periodCount * dt_ * tc * (double)periodSize + \
				(double)k * dt_ * tc;
	
			/////////////////////////////////////////////////////////////
			//Calculation of the droplet shape as a result of evaporation
			/////////////////////////////////////////////////////////////
			double sumC = 0.0;
			for (n = 0; n <= N; n++) {
				Ha[n] = 1.0 / (1.0 + exp(-4000.0 * (Cg - Cn[n] - 0.005)));
				J[n] = J0_ * (1.0 - C[n] * C[n] / Cg2) / (kappa + h[n]);
				sumC += C[n];
				hTMP[n] = hn[n] = h[n] - J[n] * dt_;
				if (hTMP[n] < 0) {
					finishFlag = true;
					str += "Negative h has been detected!\n";
					cout << "Negative h has been detected!" << endl;
				}
				// Calculate the viscosity.
				if (C[n] <= Cg)//Dependence of dimensionless viscosity on the mass fraction.
					eta_[n] = exp(S * C[n] / (1.0 - K * C[n]));// Mooney’s formula.
				else
					eta_[n] = eta_Cg;
			}

			hTMP[0] = hn[0] = (4.0 * hTMP[1] - hTMP[2]) / 3.0;//Second-order approximation
			hTMP[N] = hn[N] = hf_;//The edge of the drop

			if (sumC / Np1 > 1.0) {
				finishFlag = true;
				str += "Sum of C overflow has been detected!\n";
				cout << "Sum of C overflow has been detected!" << endl;
			}

			//If a dry area appears on the substrate or something wrong, the calculation should be finished.
			if (finishFlag) break;


			// P and u_a should be calculated after droplet shape changing 
			//as a result of evaporation. Now use an iterative method for calculations P, u_a, and h.
			int iterCount = 1;
			maxError = ErrorLimit + 1.0;
			while ((iterCount < iterNum) && (maxError > ErrorLimit)) {
				maxError = 0.0;

				///////////////////////////////////
				//Calculate the capillary pressure
				///////////////////////////////////

				double pressurePreviousIteration = 0.0;
				for (n = 1; n < N; n++) {
					pressurePreviousIteration = P[n];
					//Central finite difference scheme
					P[n] = ((r[n] + r[n + 1]) * (hn[n] - hn[n + 1]) \
						- (r[n] + r[n - 1]) * (hn[n - 1] - hn[n])) / (r[n] * dr_dr_2Ca);
					//Explicit Relaxation https://www.youtube.com/watch?v=GSsv2ncNJN8
					P[n] = pressurePreviousIteration + Relaxation * (P[n] - pressurePreviousIteration);
				}

				//Boundary condition dP/dr=0
				P[0] = (4.0 * P[1] - P[2]) / 3.0;//Second-order approximation
				P[N] = 4.0 * P[Nm1] - 3.0 * P[Nm2];//Second-order approximation

				/////////////////////////////////////////
				//Calculate the depth-averaged velocity 
				/////////////////////////////////////////

				for (n = 1; n < N; n++) //Central finite difference scheme
					u_a[n] = - Ha[n] * h[n] * h[n] * (P[n + 1] - P[n - 1]) / (dr_6 * eta_[n]);

				//Linear interpolation of near-border nodes to remove oscillations
				u_a[1] = 0.5 * u_a[2];
				u_a[Nm1] = 0.5 * u_a[Nm2];

				//The boundary conditions have already been satisfied.
				//u_a[0] = 0.0;
				//u_a[N] = 0.0;

				/////////////////////////////////////////////////////////////////
				/*Calculate the droplet profile as a result of
				convective mass transfer.*/
				/////////////////////////////////////////////////////////////////
				
				//Solve the equation based on the conservation law of the solution's mass.
				double A[Np1], B[Np1], C[Np1], D[Np1];
				for (n = 0; n <= N; n++) 
					A[n] = B[n] = C[n] = D[n] = 0.0;
					
				//Boundary condition (B.C.)
				B[0] = 1.0; C[0] = -1.0; D[0] = 0.0;
				A[N] = 0.0; B[N] = 1.0; D[N] = hf_;
				//Matrix coefficients
				for (n = 1; n < N; n++) {
					A[n] = - halfDt_ * r[n - 1] * u_a[n - 1] / dr_ / r[n];
					B[n] = 1.0;
					C[n] = halfDt_ * r[n + 1] * u_a[n + 1] / dr_ / r[n];
					D[n] = hTMP[n];
					if (B[n] < abs(A[n]) + abs(C[n]))//abs(B[n]) in general case.
						cout << "The condition of diagonal predominance is violated." << endl;
				}
				
				tridiag(A, B, C, D, hn);

				double newError = 0.0;
				for (int n = 1; n < N; n++) {
					double rhu = r[n + 1] * hn[n + 1] * u_a[n + 1] - r[n - 1] * hn[n - 1] * u_a[n - 1];
					double dtConvectionTerm = halfDt_ * rhu / (dr_ * r[n]);
					newError = abs(hn[n] - hTMP[n] + dtConvectionTerm);
					if (newError > maxError) maxError = newError;
				}

				//if ((k == 1)&&(iterCount == 2)) {//debugging
				//	cout << "maxError=" << maxError << endl;
				//	finishFlag = true;
				//	break;
				//}	

				iterCount++;
			}
			
			//Calculation of the mass fraction
			for (n = 1; n < N; n++) {
				//Convection
				Cn[n] = C[n] - halfDt_ * (abs(u_a[n]) + u_a[n]) * (C[n] - C[n - 1]) / dr_ + \
					halfDt_ * (abs(u_a[n]) - u_a[n]) * (C[n + 1] - C[n]) / dr_;
				//Diffusion
				Cn[n] += /* (Cn[n] < Cg) * */ Ha[n] * (halfDt_ / (r[n] * hn[n] * dr_Pe)) * \
					((r[n + 1] * hn[n + 1] + r[n] * hn[n]) * (C[n + 1] - C[n]) / dr_ - \
						(r[n] * hn[n] + r[n - 1] * hn[n - 1]) * (C[n] - C[n - 1]) / dr_);
				//Evaporation
				Cn[n] += J[n] * C[n] * dt_ / hn[n];
			}
			/*Boundary conditions.
			In the center due to axial symmetry.
			dC/dr=0.*/
			Cn[0] = (4.0 * Cn[1] - Cn[2]) / 3.0; // second-order approximation

			/*Near the edge of the drop due to the lack of mass transfer of
			the substance outside the liquid. dC/dr=0.*/
			//Cn[N] = 4.0 * Cn[Nm1] - 3.0 * Cn[Nm2];
			
			//The maximum concentration value at the edge for pinning the contact line.
			Cn[N] = Cg;

			/////////////////////////////////////////////////////////////////
			//Preparing for the next time step
			/////////////////////////////////////////////////////////////////
			for (n = 0; n <= N; n++) {
				h[n] = hn[n];
				C[n] = Cn[n];
			}
		}

		//////////////////////////////////////////////
		//Solving the problem of sawtooth oscillations
		//////////////////////////////////////////////

		//Restore the pressure from the depth-averaged radial velocity to eliminate the sawtooth oscillation.
		double sum = 0.0;// because u_a[0] is 0
		P[0] = 0.5 * (P[1] + P[2]);
		for (n = 1; n < N; n++) {
			sum += u_a[n] * eta_[n]  / (h[n] * h[n]);
			P[n] = -3.0 * sum * dr_ + P[0];
		}
		P[N] = 4.0 * P[Nm1] - 3.0 * P[Nm2];//Second-order approximation

		//Filter for velocity and droplet shape
		double filterVel[Np1], filterShape[Np1];
		for (int n = 1; n < N; n++) {
			filterVel[n] = 0.25 * u_a[n + 1] + 0.5 * u_a[n] + 0.25 * u_a[n - 1];
			filterShape[n] = 0.25 * h[n + 1] + 0.5 * h[n] + 0.25 * h[n - 1];
		}

		for (int n = 1; n < N; n++) {
			u_a[n] = filterVel[n];
			hn[n] = h[n] = filterShape[n];
		}
		hn[0] = h[0] = (4.0 * h[1] - h[2]) / 3.0;

		
		/////////////////////////////////////////////////////////////////
		//Two-dimensional flow.
		/////////////////////////////////////////////////////////////////
		
		for (n = 1; n < N; n++)
			for (m = 1; m <= M; m++)
				if (z[m] <= hn[n]) {
					U[m][n] = Ha[n] * (0.5 * z[m] * z[m] - hn[n] * z[m]) * \
						(P[n + 1] - P[n - 1]) / dr_2/ eta_[n];
					double xi_p = z[m] * z[m] * z[m] / 6.0 - 0.25 * (hn[n + 1] + hn[n]) * z[m] * z[m];
					double xi_m = z[m] * z[m] * z[m] / 6.0 - 0.25 * (hn[n] + hn[n - 1]) * z[m] * z[m];
					double r_p = 0.25 * (r[n] + r[n + 1]) * (Ha[n] + Ha[n+1]);
					double r_m = 0.25 * (r[n] + r[n - 1]) * (Ha[n] + Ha[n-1]);
					double eta_p = 0.5 * (eta_[n] + eta_[n + 1]);
					double eta_m = 0.5 * (eta_[n] + eta_[n - 1]);
					W[m][n] = r_p * xi_p * (P[n + 1] - P[n]) / (dr_ * eta_p) - \
						r_m * xi_m * (P[n] - P[n - 1]) / (dr_ * eta_m);
					W[m][n] /= (r[n] * dr_);
					W[m][n] *= -1.0;
				}
				else {
					W[m][n] = 0.0;
					U[m][n] = 0.0;
				}
		//B.C. dW/dr=0 for r=0 (axial symmetry) 
		for (m = 1; m <= M; m++)
			if (z[m] <= hn[0]) //second-order approximation
				W[m][0] = (4.0 * W[m][1] - W[m][2]) / 3.0;
			else
				W[m][0] = 0.0;

		///////////////////////////////////////
		//Writing data to files
		///////////////////////////////////////

		SaveDataToFile(currentModelTime);
		cout << "maxError=" << maxError << endl;

		//////////////////////////////////////////////
		//Checking the completion of the calculation
		//////////////////////////////////////////////

		if (finishFlag) {
			cout << "The program is completed earlier!" << endl;
			cout << "Model time = " << currentModelTime << " seconds." << endl;
			str += "The program is completed earlier!\n";
			str += "Model time = " + std::to_string(t) + " seconds.\n";
			break;
		}
			
		periodCount++;
	}
	///////////////////////////////////////////////
	//Finish of the time periods
	///////////////////////////////////////////////

	//Check mass
	double massFinal = 0.0;
	for (int n = 0; n <= N; n++)
		massFinal += hn[n] * Cn[n];

	cout << "The mass of the substance at the end: " << massFinal << endl;

	int time = (clock() - t) / CLOCKS_PER_SEC;

	///////////////////////////////////////////////
	//Writing parameters to files.
	///////////////////////////////////////////////
	fout.open("parameters.txt");
	if (fout.is_open()) {
		fout << "Ca=" << Ca << endl;
		fout << "Pe=" << Pe << endl;
		fout << "dt_=" << dt_ << endl;
		fout << "dr_=" << dr_ << endl;
		fout << "tc=" << tc << " s" << endl;
		fout << "tf=" << tf << " s" << endl;
		fout << "J0=" << J0 * 1.0e3 << " g/(s*m^2)" << endl;
		fout << "N=" << N << endl;
		fout << "M=" << M << endl;
		fout << "The mass of the substance at the beginning: " << mass << endl;
		fout << "The mass of the substance at the end: " << massFinal << endl;
		fout << str;
		fout << "Execution time: " << time << " s" << endl;
	}
	else
		cout << "parameters.txt is not created!" << endl;
	fout.close();

	cout << "Execution time:" << time << " s" << endl;
	cout << "The program has been completed. Press Enter to exit...";

	std::cin.get();

	return 0;
}



//Tridiagonal matrix algorithm, also known as the Thomas algorithm.
inline void tridiag(double a[], double b[], double c[], double d[], double x[])
{
	//Here, a, b, and c are the tridiagonal matrix coefficient in A, 
	//A is the two-dimensional matrix,
	//d is the vector of right parts, x is the unknown vector.
	//The system of equations in matrix form is A*x = d.

	double coef[2][Np1];//new coefficients
	if (b[0] * b[0] <= 1.0e-18) {
		cout << "Error in tridiag: b[0] = 0." << endl;
		system("pause");
		exit(1);
	}
	//Forward substitution
	coef[0][0] = -c[0] / b[0];
	coef[1][0] = d[0] / b[0];
	double tmp = 0.0;
	for (int n = 1; n < Np1; n++) {
		tmp = b[n] + a[n] * coef[0][n - 1];
		if (tmp * tmp <= 1.0e-18) {
			cout << "Error in tridiag: division by zero." << endl;
			system("pause");
			exit(1);
		}
		coef[0][n] = -c[n] / tmp;
		coef[1][n] = (d[n] - a[n] * coef[1][n - 1]) / tmp;
	}
	//Back substitution
	x[N] = coef[1][N];
	for (int n = Nm1; n >= 0; n--)
		x[n] = coef[0][n] * x[n + 1] + coef[1][n];
}


inline void SaveDataToFile(double Time) {
	std::string fname = "time_";
	std::string strTime = std::to_string((int)Time);
	fname += strTime + "s.txt";
	fout.open(fname);
	if (fout.is_open()) {
		std::string outStr = "";
		//Here,the text markup fits for OriginPro 8.1.
		// \i() - italics, \g(m) - greek letter, e.g. \g(m) means $\mu$ (LaTeX notation), 
		// \-() - subscript, and \+() - superscript.
		// You can change the text markup for any other data processing and visualization program.
		outStr += "\\i(r)\t\\i(h)\t\\i(P)\t<\\i(u)>\t\\i(J)\t\\g(eta)\t\\i(C)\n";
		outStr += "mm\tmm\tPa\t\\g(m)m/s\tg/(sm\\+(2))\tsPa\n";
		outStr += "\\i(t) = " + strTime + " s\t" + "\\i(t) = " + strTime + " s\t" + "\\i(t) = " + strTime + " s\t" + \
			"\\i(t) = " + strTime + " s\t" + "\\i(t) = " + strTime + " s\t" + "\\i(t) = " + strTime + " s\t" + \
			"\\i(t) = " + strTime + " s\n";
		for (n = 0; n <= N; n++)
			outStr += std::to_string(r[n] * R * 1.0e3) + "\t" + std::to_string(hn[n] * h0 * 1.0e3) + "\t" + \
			std::to_string(P[n] * eta * uc * R / (h0 * h0)) + \
			"\t" + std::to_string(u_a[n] * uc * 1.0e6) + "\t" + std::to_string(J[n] * Jc * 1.0e3) + "\t" + \
			 std::to_string(eta_[n] * eta) + "\t" + std::to_string(Cn[n]) + "\n";
		fout << outStr;
		fout.close();
		cout << "File " << fname << " created." << endl;
	}
	else
		cout << fname << " is not created!" << endl;

	std::string fnameU = "U_" + fname;
	fout.open(fnameU);
	if (fout.is_open()) {
		std::string outStr = "";
		for (m = 0; m <= M; m++) {
			for (n = 0; n <= N; n++)//The velocity recorded in microns per second.
				outStr += std::to_string(U[m][n] * uc * 1.0e6) + "\t";
			outStr += "\n";
		}
		fout << outStr;
		fout.close();
	}
	else
		cout << fnameU << " is not created!" << endl;

	std::string fnameW = "W_" + fname;
	fout.open(fnameW);
	if (fout.is_open()) {
		std::string outStr = "";
		for (m = 0; m <= M; m++) {
			for (n = 0; n <= N; n++)//The velocity recorded in microns per second.
				outStr += std::to_string(W[m][n] * uc * eps * 1.0e6) + "\t";
			outStr += "\n";
		}
		fout << outStr;
		fout.close();
	}
	else
		cout << fnameW << " is not created!" << endl;
}
