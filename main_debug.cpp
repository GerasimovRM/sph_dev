#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>

using namespace std;


int main()
{
	// Parameters
	/*=====================================================================================================*/
	int count_of_steps = 10;
	
	double mesh_size = 0.005; // m
	double h = mesh_size * 1.5; // m
	double h2 = h * h;
	double sound_velocity = 10; // m/s
	double density_0 = 1000.; // kg/m^3
	double heat_capacity_volume = 700.; // [Cv]: Dj / (K * m^3)
	double viscosity = 1.;
	//double temperature_0 = 300.; // K
	//double temperature_melting = 1688.; // K
	//double lambda = 16.76;
	double volume_particle = mesh_size * mesh_size * mesh_size; // m^3
	double mass_particle = volume_particle * density_0 / 10; // kg

	//double d_coeff = lambda / (heat_capacity_volume * density_0);
	double dt = 0.01 * h / sound_velocity;
	printf("%e %e\n", mass_particle, dt);

	ofstream out_file("output.xyz");
	/*=====================================================================================================*/

	// Data variables
	/*=====================================================================================================*/
	int a = 2;
	int b = 2;
	int c = 2;
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

	fill(Fx, Fx + N, 0.);
	fill(Fy, Fy + N, 0.);
	//fill(Fz, Fz + N, -9.81);
	fill(Fz, Fz + N, 0.);

	fill(energy, energy + N, 0.);
	fill(pressure, pressure + N, 0.);
	fill(density, density + N, 0.);

	fill(mass, mass + N, mass_particle);

	fill(d_density, d_density + N, 0.);
	fill(d_energy, d_energy + N, 0.);

	/*Density*/
	for (i = 0; i < N; ++i)
	{
		wfd = 2.1541870227086614782 / (h * h * h);
		density[i] = mass[i] * wfd;
		
		//if (i == 0)
		//{
		//	printf("rho[0]:=%f\n", density[i]);
		//	printf("%f %f %f\n", x[i], y[i], z[i]);
		//}
		
		for (j = 0; j < N; ++j)
		{
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 < h2)
			{
				wfd = 1. - distance2 * ihsq;
				
				wfd = wfd * wfd;
				wfd = wfd * wfd;

				wfd = 2.1541870227086614782e0 * wfd * ihsq * ih;

				density[i] += mass[j] * wfd;
				
				//if (i == 0)
				//{
				//	printf("%f\n", wfd);
				//	printf("rho[0]:=%f\n", density[i]);
				//}
			}
			//if (i == 40)
				//cout << density[i] << endl;
		}
		
		
	}
	//system("pause");
	/*=====================================================================================================*/

	// Time integrartion
	/*=====================================================================================================*/
	
	/* Step 0*/
	printf("taitwater!\n");
	for (i = 0; i < N; ++i)
	{
		tmp = density[i] / density_0;
		fi = tmp * tmp * tmp;
		fi = (sound_velocity * sound_velocity * density_0 / 7) * (fi * fi * tmp - 1) / (density[i] * density[i]);

		for (j = 0; j < N; ++j)
		{
			if (i == j)
				continue;
			del_x = x[i] - x[j];
			del_y = y[i] - y[j];
			del_z = z[i] - z[j];

			distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
			if (distance2 < h2)
			{
				distance = sqrt(distance2);

				wfd = h - distance;
				wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;

				tmp = density[j] / density_0;
				fj = tmp * tmp * tmp;
				fj = (sound_velocity * sound_velocity * density_0 / 7) * (fj * fj * tmp - 1) / (density[j] * density[j]);

				del_Vx = Vx_est[i] - Vx_est[j];
				del_Vy = Vy_est[i] - Vy_est[j];
				del_Vz = Vz_est[i] - Vz_est[j];

				delVdelR = del_x * del_Vx + del_y * del_Vy + del_z * del_Vz;
				if (delVdelR < 0)
				{
					mu = h * del_x / (distance2 + 0.01 * h * h);
					fvisc = -viscosity * (sound_velocity + sound_velocity) * mu / (density[i] + density[j]);
				}
				else
				{
					fvisc = 0.;
				}

				fpair = -mass[i] * mass[j] * (fi + fj + fvisc) * wfd;
				//printf("%e\n", fpair);

				Fx[i] += del_x * fpair;
				Fy[i] += del_y * fpair;
				Fz[i] += del_z * fpair;

				d_density[i] += mass[j] * delVdelR * wfd;
				
				if (i == 0)
				{
					printf("%e %e %e %e %e %e\n", x[j], y[j], z[j], fpair, delVdelR, fvisc);
				}
			}
			
			
		}
		if (i == 0)
			printf("%e %e %e\n", Fx[i], Fy[i], Fz[i]);
		
	}

		/*Step 1: Initial integration*/
	for (time = 0; time < count_of_steps; ++time)
	{
		//cout << time << endl;
		
		/*Write file on 0 step*/
		if (time == 0)
		{
			out_file << N << "\n";
			out_file << "Properties=i:I:1:pos:R:3:Force:R:3:Density:R:1\n";

			for (int i = 0; i < N; ++i) {
				out_file << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t";
				out_file << Fx[i] << "\t" << Fy[i] << "\t" << Fz[i] << "\t";
				out_file << density[i] << endl;
			}
		}
		printf("inital integrate!\n");
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
			
			if (i == 0)
			{
				printf("%f %f %f %f\n", Fx[i], Fy[i], Fz[i], d_density[i]);
			}

			//if (z[i] < 0.)
				//z[i] = 0.;
		}

		/*Step 2: Calculate forces, d_density, d_energy*/
		fill(Fx, Fx + N, 0.);
		fill(Fy, Fy + N, 0.);
		//fill(Fz, Fz + N, -9.81);
		fill(Fz, Fz + N, 0.);

		fill(d_density, d_density + N, 0.);
		fill(d_energy, d_energy + N, 0.);
		printf("taitwater!\n");
		for (i = 0; i < N; ++i)
		{
			tmp = density[i] / density_0;
			fi = tmp * tmp * tmp;
			fi = (sound_velocity * sound_velocity * density_0 / 7) * (fi * fi * tmp - 1) / (density[i] * density[i]);

			for (j = 0; j < N; ++j)
			{
				if (i == j)
					continue;
				del_x = x[i] - x[j];
				del_y = y[i] - y[j];
				del_z = z[i] - z[j];

				distance2 = del_x * del_x + del_y * del_y + del_z * del_z;
				if (distance2 < h2)
				{
					distance = sqrt(distance2);
	
					wfd = h - distance;
					wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;

					tmp = density[j] / density_0;
					fj = tmp * tmp * tmp;
					fj = (sound_velocity * sound_velocity * density_0 / 7) * (fj * fj * tmp - 1) / (density[j] * density[j]);

					del_Vx = Vx_est[i] - Vx_est[j];
					del_Vy = Vy_est[i] - Vy_est[j];
					del_Vz = Vz_est[i] - Vz_est[j];

					delVdelR = del_x * del_Vx + del_y * del_Vy + del_z * del_Vz;
					if (delVdelR < 0.)
					{
						mu = h * delVdelR / (distance2 + 0.01 * h * h);
						fvisc = -viscosity * (sound_velocity + sound_velocity) * mu / (density[i] + density[j]);
					}
					else
					{
						fvisc = 0;
					}

					fpair = -mass[i] * mass[j] * (fi + fj + fvisc) * wfd;

					Fx[i] += del_x * fpair;
					Fy[i] += del_y * fpair;
					Fz[i] += del_z * fpair;

					d_density[i] += mass[j] * delVdelR * wfd;
					if (i == 0)
					{
						printf("%e %e %e %e %e %e\n", x[j], y[j], z[j], fpair, delVdelR, fvisc);
					}
				}
				
			}
			if (i == 0)
				printf("%e %e %e\n", Fx[i], Fy[i], Fz[i]);
		}

		/*Step 3: Final integration*/
		printf("final integrate!\n");
		for (i = 0; i < N; ++i) 
		{
			Vx[i] += Fx[i] * dt / (2 * mass[i]);
			Vy[i] += Fy[i] * dt / (2 * mass[i]);
			Vz[i] += Fz[i] * dt / (2 * mass[i]);

			energy[i] += d_energy[i] * dt / 2;
			density[i] += d_density[i] * dt / 2;
			
			if (i == 0)
			{
				printf("%f %f %f %f\n", Fx[i], Fy[i], Fz[i], d_density[i]);
			}
		}
		/*Read data*/
		if (time % 1 == 0)
		{
			out_file << N << "\n";
			out_file << "Properties=i:I:1:pos:R:3:Force:R:3:Density:R:1\n";

			for (i = 0; i < N; ++i) {
				out_file << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t";
				out_file << Fx[i] << "\t" << Fy[i] << "\t" << Fz[i] << "\t";
				out_file << density[i] << endl;
			}
		}
	}
	/*=====================================================================================================*/

	out_file.close();
	//system("pause");
	return 0;
}