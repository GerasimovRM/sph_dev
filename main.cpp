#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>

using namespace std;

void read_data_xyz(
					ofstream& out_file,
					double* x,
					double* y,
					double* z,
					double* Fx,
					double* Fy,
					double* Fz,
					double* density,
					double* energy,
					int N) {
						
	int i;
	
	out_file << N << "\n";
	out_file << "Properties=i:I:1:pos:R:3:Force:R:3:Density:R:1:Energy:R:1\n";

	for (i = 0; i < N; ++i) {
		out_file << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t";
		out_file << Fx[i] << "\t" << Fy[i] << "\t" << Fz[i] << "\t";
		out_file << density[i] << "\t";
		out_file << energy[i] << endl;
	}
}

void rho_sum(
				double* x,
				double* y,
				double* z,
				double*& density,
				double* mass,
				double h,
				int N) {
	double del_x, del_y, del_z;
	double distance2, wfd;
	double ih = 1. / h;
	double ihsq = ih * ih;
	double h2 = h * h;
	int i, j;
	
	for (i = 0; i < N; ++i)
	{
		wfd = 2.1541870227086614782 / (h * h * h);
		density[i] = mass[i] * wfd;


		for (j = 0; j < N; ++j)
		{
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 <  h2)
			{			
				wfd = 1. - distance2 * ihsq;
				wfd = wfd * wfd;
				wfd = wfd * wfd;
				wfd = 2.1541870227086614782e0 * wfd * ihsq * ih;		

				density[i] += mass[j] * wfd;
				if (i == 0)
				{
					printf("%f\n", wfd);
					printf("rho[0]:=%f\n", density[i]);
				}
			}
			
		}
	}
}

void taitwater(
					double*& x,
					double*& y,
					double*& z,
					double*& Vx_est,
					double*& Vy_est,
					double*& Vz_est,
					double*& Fx,
					double*& Fy,
					double*& Fz,
					double*& mass,
					double*& density,
					double*& energy,
					double*& d_density,
					double*& d_energy,
					double h,
					double density_0,
					double sound_velocity,
					double viscosity,
					int N) {
		
		
	int i;
	int j;
	
	double distance;
	double distance2;
	
	double del_x;
	double del_y;
	double del_z;
	
	double del_Vx;
	double del_Vy;
	double del_Vz;
	
	double h2 = h * h;
	double ih = 1. / h;
	double ihsq = ih * ih;
	
	double wfd;
	double tmp, fi, fj;
	double fpair;
	double delVdelR;
	double mu;
	double fvisc;
	double deltaE;

	fill(Fx, Fx + N, 0.);
	fill(Fy, Fy + N, 0.);
	//fill(Fz, Fz + N, -9.81);
	fill(Fz, Fz + N, 0.);
	fill(d_density, d_density + N, 0.);
	fill(d_energy, d_energy + N, 0.);
	printf("taitwater!\n");
	for (i = 0; i < N; ++i) {
		tmp = density[i] / density_0;
		fi = tmp * tmp * tmp;
		fi = (sound_velocity * sound_velocity * density_0 / 7) * (fi * fi * tmp - 1) / (density[i] * density[i]);

		for (j = 0; j < N; ++j) {
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 < h2) {				
				distance = sqrt(distance2);
				wfd = h - distance;
				wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
				

				//wfd = 2.08890862808113 * (-12 * distance + 24 * distance2 - 12 * distance * distance2) / (h * h * h * h);
				/*
				if (distance <= 1) {
					tmp = 2.25 * distance2 - 3 * distance;
				}
				else if (distance <= 2) {
					tmp = 2 - distance;
					tmp = -0.75 * tmp * tmp;
				}
				else {
					continue;
				}
				wfd = 0.3183098861837906 * tmp / (h * h * h);
				*/

				tmp = density[j] / density_0;
				fj = tmp * tmp * tmp;
				fj = (sound_velocity * sound_velocity * density_0 / 7) * (fj * fj * tmp - 1) / (density[j] * density[j]);

				del_Vx = Vx_est[i] - Vx_est[j];
				del_Vy = Vy_est[i] - Vy_est[j];
				del_Vz = Vz_est[i] - Vz_est[j];

				delVdelR = del_x * del_Vx + del_y * del_Vy + del_z * del_Vz;
				if (delVdelR < 0.) {
					mu = h * delVdelR / (distance2 + 0.01 * h * h);
					fvisc = -viscosity * (sound_velocity + sound_velocity) * mu / (density[i] + density[j]);
				}
				else {
					fvisc = 0;
				}

				fpair = -mass[i] * mass[j] * (fi + fj + fvisc) * wfd;
				deltaE = -0.5 * fpair * delVdelR;

				Fx[i] += del_x * fpair;
				Fy[i] += del_y * fpair;
				Fz[i] += del_z * fpair;

				d_density[i] += mass[j] * delVdelR * wfd;
				d_energy[i] += deltaE;
				
				if (i == 0)
				{
					printf("%e %e %e %e %e %e\n", x[j], y[j], z[j], fpair, delVdelR, fvisc);
				}
			}
		}
		if (i == 0)
			printf("%e %e %e\n", Fx[i], Fy[i], Fz[i]);
	}
}

void taitwater_morris(
						double*& x,
						double*& y,
						double*& z,
						double*& Vx_est,
						double*& Vy_est,
						double*& Vz_est,
						double*& Fx,
						double*& Fy,
						double*& Fz,
						double*& mass,
						double*& density,
						double*& energy,
						double*& d_density,
						double*& d_energy,
						double h,
						double density_0,
						double sound_velocity,
						double viscosity,
						int N) {
		
		
	int i, j;
	
	double distance;
	double distance2;
	
	double del_x;
	double del_y;
	double del_z;
	
	double del_Vx;
	double del_Vy;
	double del_Vz;
	
	double h2 = h * h;
	double ih = 1. / h;
	double ihsq = ih * ih;
	
	double wfd;
	double tmp, fi, fj;
	double fpair;
	double delVdelR;
	double mu;
	double fvisc;
	double deltaE;

	fill(Fx, Fx + N, 0.);
	fill(Fy, Fy + N, 0.);
	//fill(Fz, Fz + N, -9.81);
	fill(Fz, Fz + N, 0.);
	fill(d_density, d_density + N, 0.);
	fill(d_energy, d_energy + N, 0.);
	printf("taitwater_morris!\n");
	for (i = 0; i < N; ++i) {
		tmp = density[i] / density_0;
		fi = tmp * tmp * tmp;
		fi = (sound_velocity * sound_velocity * density_0 / 7) * (fi * fi * tmp - 1) / (density[i] * density[i]);

		for (j = 0; j < N; ++j) {
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 < h2) {				
				distance = sqrt(distance2);
				wfd = h - distance;
				wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
				

				//wfd = 2.08890862808113 * (-12 * distance + 24 * distance2 - 12 * distance * distance2) / (h * h * h * h);
				/*
				if (distance <= 1) {
					tmp = 2.25 * distance2 - 3 * distance;
				}
				else if (distance <= 2) {
					tmp = 2 - distance;
					tmp = -0.75 * tmp * tmp;
				}
				else {
					continue;
				}
				wfd = 0.3183098861837906 * tmp / (h * h * h);
				*/

				tmp = density[j] / density_0;
				fj = tmp * tmp * tmp;
				fj = (sound_velocity * sound_velocity * density_0 / 7) * (fj * fj * tmp - 1) / (density[j] * density[j]);

				del_Vx = Vx_est[i] - Vx_est[j];
				del_Vy = Vy_est[i] - Vy_est[j];
				del_Vz = Vz_est[i] - Vz_est[j];

				delVdelR = del_x * del_Vx + del_y * del_Vy + del_z * del_Vz;
				
				fvisc = 2 * viscosity / (density[i] * density[j]);
				fvisc *= mass[i] * mass[j] * wfd;

				fpair = -mass[i] * mass[j] * (fi + fj) * wfd;
				deltaE = -0.5 * (fpair * delVdelR + fvisc * (del_Vx * del_Vx + del_Vy * del_Vy + del_Vz * del_Vz));

				Fx[i] += del_x * fpair + del_Vx * fvisc;
				Fy[i] += del_y * fpair + del_Vy * fvisc;
				Fz[i] += del_z * fpair + del_Vz * fvisc;

				d_density[i] += mass[j] * delVdelR * wfd;
				d_energy[i] += deltaE;
				
				if (i == 0)
				{
					printf("%e %e %e %e %e %e\n", x[j], y[j], z[j], fpair, delVdelR, fvisc);
				}
			}
		}
		if (i == 0)
			printf("%e %e %e\n", Fx[i], Fy[i], Fz[i]);
	}
}

void heatconduction(
						double*& x,
						double*& y,
						double*& z,
						double*& mass,
						double*& density,
						double*& energy,
						double*& d_energy,
						double h,
						double alpha,
						int N) {
							
	int i, j;
	
	double distance;
	double distance2;
	
	double del_x;
	double del_y;
	double del_z;

	double h2 = h * h;
	double ih = 1. / h;
	double ihsq = ih * ih;
	
	double wfd;

	double deltaE;
	
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 < h2) {				
				distance = sqrt(distance2);
				wfd = h - distance;
				wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
			
				deltaE = 2. * mass[i] * mass[j] / (mass[i] + mass[j]);
				deltaE *= (density[i] + density[j]) / (density[i] * density[j]);
				deltaE *= alpha * (energy[i] - energy[j]) * wfd;
				
				d_energy[i] += deltaE;
			}
		}
	}
	
}

int main() {
	// Parameters
	/*=====================================================================================================*/
	int count_of_steps = 500;
	
	double mesh_size = 0.005; // m
	double h = 1.5 * mesh_size; // m
	double h2 = h * h;
	double sound_velocity = 10; // m/s
	double density_0 = 1000.; // kg/m^3
	double heat_capacity_volume = 700.; // [Cv]: Dj / (K * m^3)
	double viscosity = 1.;
	//double temperature_0 = 300.; // K
	//double temperature_melting = 1688.; // K
	//double lambda = 16.76;
	double volume_particle = mesh_size * mesh_size * mesh_size; // m^3
	double mass_particle = volume_particle * density_0; // kg
	//printf("%f", mas_particle);
	//double d_coeff = lambda / (heat_capacity_volume * density_0);
	double d_coeff = 1e-4;
	double dt = 0.01 * h / sound_velocity;
	printf("%e %e\n", mass_particle, dt);

	ofstream out_file("output.xyz");
	/*=====================================================================================================*/

	// Data variables
	/*=====================================================================================================*/
	int a = 10;
	int b = 10;
	int c = 10;
	int N = a * b * c;

	double* x = new double[N];
	double* y = new double[N];
	double* z = new double[N];

	double* Vx = new double[N];
	double* Vy = new double[N];
	double* Vz = new double[N];

	double* Vx_est = new double[N];
	double* Vy_est = new double[N];
	double* Vz_est = new double[N];

	double* Fx = new double[N];
	double* Fy = new double[N];
	double* Fz = new double[N];

	double* energy = new double[N];
	double* density = new double[N];
	double* pressure = new double[N];

	double* mass = new double[N];

	double* d_energy = new double[N];
	double* d_density = new double[N];
	/*=====================================================================================================*/

	/*Buffers*/
	/*=====================================================================================================*/
	int time, i, j;
	double distance, distance2;
	double del_x, del_y, del_z;
	double del_Vx, del_Vy, del_Vz;
	double ih = 1. / h, ihsq = ih * ih, wfd;
	double tmp, fi, fj;
	double fpair;
	double delVdelR;
	double mu;
	double fvisc;
	/*=====================================================================================================*/

	// Data initialize
	/*=====================================================================================================*/
	int numb;
	for (int i = 0; i < a; ++i)
		for (int j = 0; j < b; ++j)
			for (int k = 0; k < c; ++k)
			{
				numb = i * (b * c) + j * c + k;

				x[numb] = mesh_size * i;
				y[numb] = mesh_size * j;
				z[numb] = mesh_size * k;
			}

	fill(Vx, Vx + N, 0.);
	fill(Vy, Vy + N, 0.);
	fill(Vz, Vz + N, 0.);

	fill(Vx_est, Vx_est + N, 0.);
	fill(Vy_est, Vy_est + N, 0.);
	fill(Vz_est, Vz_est + N, 0.);

	fill(mass, mass + N, mass_particle);

	fill(pressure, pressure + N, 0.);
	fill(energy, energy + N, 0.);
	fill(density, density + N, 0.);
	
	rho_sum(x, y, z, density, mass, h, N);
	
	
	
	
	
	/*=====================================================================================================*/


	// Time integrartion
	/*=====================================================================================================*/
	
	/* Step 0*/
	taitwater_morris(x, y, z, Vx_est, Vy_est, Vz_est, Fx, Fy, Fz, mass, density, energy, d_density, d_energy, h, density_0, sound_velocity, viscosity, N);
		/*Step 1: Initial integration*/
	for (time = 0; time < count_of_steps; ++time)
	{
		cout << time << endl;
		
		/*Write file on 0 step*/
		/*
		if (time == 0)
			read_data_xyz(out_file, x, y, z, Fx, Fy, Fz, density, energy, N);
		*/
		
		/*Step 1: Initial integration*/
		for (i = 0; i < N; ++i)
		{
			energy[i] += d_energy[i] * dt / 2.; // half-step update of particle internal energy
			density[i] += d_density[i] * dt / 2.; // ... and density

			// extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
			Vx_est[i] = Vx[i] + Fx[i] * dt / mass[i];
			Vy_est[i] = Vy[i] + Fy[i] * dt / mass[i];
			Vz_est[i] = Vz[i] + Fz[i] * dt / mass[i];

			Vx[i] += Fx[i] * dt / (2 * mass[i]);
			Vy[i] += Fy[i] * dt / (2 * mass[i]);
			Vz[i] += Fz[i] * dt / (2 * mass[i]);

			x[i] += Vx[i] * dt;
			y[i] += Vy[i] * dt;
			z[i] += Vz[i] * dt;	
		}

		/*Step 2: Calculate forces, d_density, d_energy*/
		taitwater_morris(x, y, z, Vx_est, Vy_est, Vz_est, Fx, Fy, Fz, mass, density, energy, d_density, d_energy, h, density_0, sound_velocity, viscosity, N);
		heatconduction(x, y, z, mass, density, energy, d_energy, d_coeff, h, N);
		
		/*Step 3: Final integration*/
		for (i = 0; i < N; ++i) 
		{
				Vx[i] += Fx[i] * dt / (2 * mass[i]);
				Vy[i] += Fy[i] * dt / (2 * mass[i]);
				Vz[i] += Fz[i] * dt / (2 * mass[i]);

				energy[i] += d_energy[i] * dt / 2;
				density[i] += d_density[i] * dt / 2;
		}
		/*Read data*/
		if (time % 10 == 0)
			read_data_xyz(out_file, x, y, z, Fx, Fy, Fz, density, energy, N);
	}
	/*=====================================================================================================*/
	
	
	// It's a DESTRUCTION TIME
	/*=====================================================================================================*/
	out_file.close();
	delete[]x;
	delete[]y;
	delete[]z;

	delete[]Vx;
	delete[]Vy;
	delete[]Vz;

	delete[]Vx_est;
	delete[]Vy_est;
	delete[]Vz_est;

	delete[]Fx;
	delete[]Fy;
	delete[]Fz;
	
	delete[]energy;
	delete[]density;
	delete[]pressure;
	
	delete[]mass;
	
	delete[]d_energy;
	delete[]d_density;
	/*=====================================================================================================*/

	
	return 0;
}