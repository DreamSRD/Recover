#include <stdio.h>
#include "math.h"
#include <complex.h>
#include <fftw3.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define M_PI 3.1415926535
#define N 8192

int k0 = 100;

double g = 9.81;



int Number_of_files = 10000;
// int Number_of_files = 1000000;
int number_of_numbers;
char Name[100];

double a = 0;
double b = 10000;
double Dk;

FILE *pfile;

fftw_complex RkCur[N];
fftw_complex VkCur[N];
fftw_complex VuCur[N];

fftw_complex dPhiu[N];
fftw_complex dPhik[N];
fftw_complex Phiu[N];
fftw_complex Phik[N];

fftw_complex dPsiu[N];
fftw_complex dPsik[N];
fftw_complex Psiu[N];
fftw_complex Psik[N];

fftw_complex dZk[N];
fftw_complex dZu[N];
fftw_complex dXTk[N];  
fftw_complex dXTu[N]; 
fftw_complex XTk[N]; 
fftw_complex XTu[N]; 

fftw_complex Yk[N]; 
fftw_complex Yu[N]; 
fftw_complex YuCur[N];
fftw_complex dYu[N];
fftw_complex dYk[N];

fftw_plan dPhiu_2_dPhik;
fftw_plan Phik_2_Phiu_;

fftw_plan dPsiu_2_dPsik;
fftw_plan Psik_2_Psiu_;

fftw_plan etax_2_etak;
fftw_plan psix_2_psik;

fftw_plan dZk_2_dZu_;
fftw_plan dXTk_2_dXTu_;
fftw_plan Yk_2_Yu_;
fftw_plan VkCur_2_VuCur_;

fftw_plan Etax_to_Etak_half;

fftw_complex ck[N];
fftw_complex cx[N];
fftw_complex Cx[N];
fftw_complex Ck[N];

fftw_complex eta_x[N];
fftw_complex eta_half_x[N/2];
fftw_complex psi_x[N];

fftw_complex eta_k[N];
fftw_complex eta_half_k[N/2];
fftw_complex psi_k[N];

fftw_plan ck_2_cx_, Cx_2_Ck;

void CheckFile(FILE* file_name)
{
	if (file_name == NULL)
	{
		printf("Cannot open or find the file \n");
		exit(0);
	}
}

void Normalize(fftw_complex *Arr_x)
{
	for (int j = 0; j < N; j++)
	{
		Arr_x[j] = Arr_x[j]/N;
	}

	Arr_x[N/2] = 0.0 + I * 0.0;
	Arr_x[0] = 0.0 + I * 0.0;
}

void fftw_plans_create()
{
	etax_2_etak = fftw_plan_dft_1d(N, &eta_x[0], &eta_k[0], FFTW_FORWARD, FFTW_ESTIMATE);
	Etax_to_Etak_half = fftw_plan_dft_1d(N/2, &eta_half_x[0], &eta_half_k[0], FFTW_FORWARD, FFTW_ESTIMATE);

	psix_2_psik = fftw_plan_dft_1d(N, &psi_x[0], &psi_k[0], FFTW_FORWARD, FFTW_ESTIMATE);

	VkCur_2_VuCur_ = fftw_plan_dft_1d(N, &VkCur[0], &VuCur[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dZk_2_dZu_ = fftw_plan_dft_1d(N, &dZk[0], &dZu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dXTk_2_dXTu_ = fftw_plan_dft_1d(N, &dXTk[0], &dXTu[0], FFTW_BACKWARD, FFTW_ESTIMATE );

	Yk_2_Yu_ = fftw_plan_dft_1d(N, &Yk[0], &Yu[0], FFTW_BACKWARD, FFTW_ESTIMATE );
	dPsiu_2_dPsik = fftw_plan_dft_1d(N, &dPsiu[0], &dPsik[0], FFTW_FORWARD, FFTW_ESTIMATE );
	Psik_2_Psiu_ = fftw_plan_dft_1d(N, &Psik[0], &Psiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );

	dPhiu_2_dPhik = fftw_plan_dft_1d(N, &dPhiu[0], &dPhik[0], FFTW_FORWARD, FFTW_ESTIMATE );

	Phik_2_Phiu_ = fftw_plan_dft_1d(N, &Phik[0], &Phiu[0], FFTW_BACKWARD, FFTW_ESTIMATE );

	ck_2_cx_ = fftw_plan_dft_1d(N,&ck[0],&cx[0],FFTW_BACKWARD,FFTW_ESTIMATE);
	Cx_2_Ck = fftw_plan_dft_1d(N,&Cx[0],&Ck[0],FFTW_FORWARD,FFTW_ESTIMATE);
}

void read_RV(int file_number)
{
	int number_of_numbers;
	char Name[40];
	sprintf(Name, "./restart_R/%d_restart", file_number);
	pfile = fopen (Name,"rb");
	CheckFile(pfile);
	number_of_numbers = fread(&RkCur[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);

	sprintf(Name, "./restart_V/%d_restart", file_number);
	pfile = fopen (Name,"rb");
	CheckFile(pfile);
	number_of_numbers = fread(&VkCur[0], sizeof(fftw_complex), N, pfile);
	fclose(pfile);

	printf("File %d has been successfully read \n", file_number);
}

// void read_eta_psi(int file_number)
// {
// 	int number_of_numbers;
// 	char Name[40];
// 	sprintf(Name, "./eta_k_%d", file_number);
// 	pfile = fopen (Name,"rb");
// 	CheckFile(pfile);
// 	number_of_numbers = fread(&Yk[0], sizeof(fftw_complex), N, pfile);
// 	fclose(pfile);

// 	sprintf(Name, "./psi_k_%d", file_number);
// 	pfile = fopen (Name,"rb");
// 	CheckFile(pfile);
// 	number_of_numbers = fread(&Psik[0], sizeof(fftw_complex), N, pfile);
// 	fclose(pfile);

// 	printf("File %d has been successfully read \n", file_number);
// }

void recover_z()
{
	int i, j;
	Dk = 2.0 * M_PI /(b - a);
	
	dZk[0] = 1.0;
	fftw_complex sum;
	for(i=1; i<N/2; i++)
	{
		sum = 0.0 + I*0.0;

		for(j=1; j<i; j++)
		{
			sum = sum + RkCur[N-j]*dZk[N-(i-j)];
		}
		dZk[N-i] = - RkCur[N-i] - sum;
	}	
	for(i=1; i<=N/2; i++)
	{
		dZk[i] = 0.0 + I*0.0;
	}

	fftw_execute(dZk_2_dZu_);

	
	for(i=1; i<N/2; i++){
		dXTk[N-i] = dZk[N-i]/2.0;
		XTk[N-i] = I*dXTk[N-i]/((double)i*Dk);
		Yk[N-i] = dZk[N-i]/2.0/((double)i*Dk);
		dXTk[i] = conj(dXTk[N-i]);
		XTk[i] = conj(XTk[N-i]);
		Yk[i] = conj(Yk[N-i]);
	}
	dXTk[0] = 0.0 + I*0.0;
	dXTk[N/2] = 0.0 + I*0.0;
	XTk[0] = 0.0 + I*0.0;
	XTk[N/2] = 0.0 + I*0.0;
	Yk[0] = 0.0 + I*0.0;
	Yk[N/2] = 0.0 + I*0.0;


	fftw_execute(dXTk_2_dXTu_);

	sum = 0.0 + I*0.0;
	for(i=1; i<N/2; i++){
		sum = sum + (double)i*Dk*Yk[i]*conj(Yk[i]);
	}
	Yk[0] = - 2.0*sum;
	

	fftw_execute(Yk_2_Yu_);

}

double Find_XTu(double u_current)
{
	fftw_complex X_tilde_u;
	X_tilde_u = XTk[0];
	for (int j = 1; j < N/2; j++)
	{
		X_tilde_u = X_tilde_u + (XTk[j] * cexp(I * Dk * j * u_current) + conj(XTk[j] * cexp(I * Dk * j * u_current)));	
	}

	return creal(X_tilde_u);
	
}

double find_yu(double u_current)
{
	fftw_complex Yu_current;

	Yu_current = Yk[0];
	for (int j = 1; j < N/2; j++)
	{
		Yu_current = Yu_current + (Yk[j] * cexp(I * Dk * j * u_current) + conj(Yk[j] * cexp(I * Dk * j * u_current)));	
	}

	return creal(Yu_current);
	
}

double Find_Psi_u(double u_current)
{
	fftw_complex Psiu_current;

	Psiu_current = Psik[0];
	for (int j = 1; j < N/2; j++)
	{
		Psiu_current = Psiu_current + (Psik[j] * cexp(I * Dk * j * u_current) + conj(Psik[j] * cexp(I * Dk * j * u_current)));	
	}

	return creal(Psiu_current);
}

void Recover_eta_psi_uniform_x(int file_number)
{
	double dx;
	dx = (b - a)/N;
	double ERR = 1.e6;
	double MIN_ERR = 1e-10;

	double dU[N];
	fftw_complex psi_u[N];
	fftw_complex y_u[N];

	char Name[100];

	FILE *pfile;

	for (int j = 0; j < N; j++)
	{
		dU[j] = j * dx;
		eta_x[j] = 0.0 + I * 0.0;
		psi_x[j] = 0.0 + I * 0.0;
	}

	int Iter = 0;

	while (ERR > MIN_ERR)
	{
		ERR = 0.0;
		for (int j = 0; j < N; j++)
		{
			dU[j] = j * dx - Find_XTu(dU[j]); // iterations to find dU
			//psi_u[j] = Find_Psi_u(dU[j]);
			y_u[j] = find_yu(dU[j]);
			ERR += fabs(creal(y_u[j]) - creal(eta_x[j]));
			eta_x[j] = y_u[j];
			//ERR += fabs(dU[j] + Find_XTu(dU[j]) - j * dx);
		}

		ERR = ERR/N;	

		Iter++;

	}

	for (int j = 0; j < N; j++)
	{
		//eta_x[j] = find_yu(dU[j]);
		psi_x[j] = Find_Psi_u(dU[j]);
	}

	printf("Uniform x-grid has been obtained, ERR = %e\n", ERR);

	fftw_execute(etax_2_etak);
	fftw_execute(psix_2_psik);
	Normalize(&eta_k[0]);
	Normalize(&psi_k[0]);

/*	for(int j = 0; j < N/2; j++)
	{
		eta_half_x[j] = eta_x[j + N/4];
	}

	fftw_execute(Etax_to_Etak_half);

	for(int j = 0; j < N/2; j++)
	{
		eta_half_k[j] = eta_half_k[j]/(N/2);
	}

	eta_half_k[0] = 0.0 + I * 0.0;
	eta_half_k[N/8] = 0.0 + I * 0.0;

	sprintf(Name, "./eta_halfed_k//eta_half_k_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);

	for (int j = N/4; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\n", 2 * (j - N/2), cabs(eta_half_k[j]));
	}

	for (int j = 0; j < N/4; j++)
	{
		fprintf(pfile, "%d\t%e\n", 2 * j, cabs(eta_half_k[j]));
	}

	fclose(pfile);*/


}

void recover_psi()
{
	int u, i;
	fftw_execute(VkCur_2_VuCur_);
	for(u=0; u<N; u++){
		dPhiu[u] = -I*VuCur[u]*dZu[u];
		dPsiu[u] = creal(dPhiu[u]);
	}
	fftw_execute(dPhiu_2_dPhik);
	fftw_execute(dPsiu_2_dPsik);

	Normalize(&dPhik[0]);
	Normalize(&dPsik[0]);

	for(i=1; i<N/2; i++){
		Phik[i] = -I*dPhik[i]/((double)i*Dk);
		Psik[i] = -I*dPsik[i]/((double)i*Dk);
		Phik[N-i] = I*dPhik[N-i]/((double)i*Dk);
		Psik[N-i] = I*dPsik[N-i]/((double)i*Dk);
	}
	Phik[0] = 0.0 + I*0.0;
	Phik[N/2] = 0.0 + I*0.0;
	Psik[0] = 0.0 + I*0.0;
	Psik[N/2] = 0.0 + I*0.0;

	fftw_execute(Phik_2_Phiu_);
	fftw_execute(Psik_2_Psiu_);

}

void Close()
{
	fftw_destroy_plan(etax_2_etak);
	fftw_destroy_plan(Etax_to_Etak_half);

	fftw_destroy_plan(psix_2_psik);

	fftw_destroy_plan(VkCur_2_VuCur_);	
	fftw_destroy_plan(dZk_2_dZu_);
	fftw_destroy_plan(dXTk_2_dXTu_);

	fftw_destroy_plan(Yk_2_Yu_);
	fftw_destroy_plan(dPsiu_2_dPsik);
	fftw_destroy_plan(Psik_2_Psiu_);

	fftw_destroy_plan(ck_2_cx_);
	fftw_destroy_plan(Cx_2_Ck);
}

void Calculate_c(int file_number)
{
	double L = b - a;
	double h = (b - a)/N;
	double Dk = 2.0 * M_PI/L;
	double w0 = sqrt(g * cabs(k0 * Dk));
	double t = file_number * 100.0;


	for (int j = 0; j < N/2; j++)
	{
		//cx[j] = sqrt(2.0)/2.0 * pow(g, 0.25) * pow(k0 * Dk, 0.25) * (eta_x[j] + I * sqrt(k0 * Dk/g) * psi_x[j]);
		ck[j] = sqrt(2.0)/2.0 * pow(g, 0.25) * pow(cabs(j * Dk), 0.25) * (eta_k[j] + I * sqrt(cabs(j * Dk)/g) * psi_k[j]);

		// ck[j] = sqrt(2.0)/2.0 * (sqrt(g/(sqrt(g * fabs(j * Dk)))) * eta_k[j] + sqrt(sqrt(g * fabs(j * Dk))/g) * I * psi_k[j]); //// ak instead of ck !
		// ck[N - j] = sqrt(2.0)/2.0 * (sqrt(g/(sqrt(g * fabs(j * Dk)))) * conj(eta_k[j]) + sqrt(sqrt(g * fabs(j * Dk))/g) * I * conj(psi_k[j]));
	}

	ck[0] = 0.0 + I * 0.0;
	ck[N/2] = 0.0 + I * 0.0;

	for (int j = N/2; j < N; j++)
	{
		ck[j] = 0.0 + I * 0.0;
		//ck[j] = sqrt(2.0)/2.0 * pow(g, 0.25) * pow(cabs((j - N) * Dk), 0.25) * (conj(eta_k[j]) + I * sqrt(cabs((j - N) * Dk)/g) * conj(psi_k[j]));
	}

	sprintf(Name, "./ck//ck_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);

	for (int j = N/2; j < N; j++)
	{
		fprintf(pfile, "%d\t%e\n", j - N, cabs(ck[j]));
	}

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\n", j, cabs(ck[j]));
	}

	fclose(pfile);

	fftw_execute(ck_2_cx_);

	for (int j = 0; j < N; j++)
	{
		Cx[j] = cx[j] * cexp(-I * k0 * Dk * (a + j * h) - I * w0 * t);
	}

	fftw_execute(Cx_2_Ck);

	for (int j = 0; j < N; j++)
	{
		Ck[j] = Ck[j]/N;
	}
}

void Print_data(int file_number)
{
	double h = (b - a)/N;
	double L = b - a;
	double Dk = 2 * M_PI/L;
	double k0_d = k0 * Dk;

	double w0 = sqrt(g * k0_d);

	//double Beta = 2.0 * k0_d * k0_d/sqrt(w0); // for dimless NLS ( coefficients chosen in such way that amplitudes does not change --- C(x,t) = Psi(x',tau))

	double x_scale = 2.0 * k0_d; // dimless x for NLS
	double Psi_scale = k0_d/sqrt(w0); // psi scaling for steepness to be unchanged

	sprintf(Name, "./Eta_Psi_x//Eta_Psi_x_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);

	//fprintf(pfile, "\t%s\t\t\t  %s\t\t\t%s\t\t\t\n", "[x]", "Re[c(x)]", "|c(x)|");

	for (int j = 0; j < N; j++)
	{
		//fprintf(pfile, "%e\t%e\t%e\n", a + j * h, creal(eta_x[j]), creal(psi_x[j]));
		fprintf(pfile, "%e\t%e\n", a + j * h, creal(eta_x[j]));
	}

	fclose(pfile);

/*	sprintf(Name, "./Eta_Psi_k//Eta_Psi_k_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);

	for (int j = N/2; j < N; j++)
	{
		fprintf(pfile, "%d\t%e\t%e\n", j - N, cabs(eta_k[j]), cabs(psi_k[j]));
	}

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\t%e\n", j, cabs(eta_k[j]), cabs(psi_k[j]));
	}

	fclose(pfile);*/

	sprintf(Name, "./qx//qx_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);
	
	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%e\t%e\t%e\n", x_scale * (a + j * h), Psi_scale * creal(Cx[j]), Psi_scale * cimag(Cx[j]));
	}

/*	for (int j = 0; j < N; j++)
	{
		fprintf(pfile, "%e\t%e\n", a + j * h, cabs(Cx[j]));
	}*/

	fclose(pfile);

/*	sprintf(Name, "./Ck//Ck_%d.dat", file_number);
	pfile = fopen (Name, "w");
	CheckFile(pfile);

	for (int j = N/2; j < N; j++)
	{
		fprintf(pfile, "%d\t%e\n", j - N, cabs(Ck[j]));
	}

	for (int j = 0; j < N/2; j++)
	{
		fprintf(pfile, "%d\t%e\n", j, cabs(Ck[j]));
	}

	fclose(pfile);*/
}

void Clean()
{
	system("exec rm -r Eta_Psi_x/*");
	system("exec rm -r Eta_Psi_k/*");
	system("exec rm -r Cx/*");
	system("exec rm -r Ck/*");
	system("exec rm -r ck/*");
}


int main()
{
	fftw_plans_create();

	Clean();

	int k = 0;

	for (int j = 1; j < 201; j++)
	{
		read_RV(j);
		recover_z();
		recover_psi();
		Recover_eta_psi_uniform_x(j);
		Calculate_c(j);
		Print_data(j);
	}

	Close();
	return 0;
}
