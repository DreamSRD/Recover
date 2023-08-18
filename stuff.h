	#include <time.h>
	
	// #define N 4096
	#define N 8192
	// #define N 16384
	// #define N 32768

	
// ================================ position of bi-soliton ====================
	double Tp = 1000.0;
	double x0;
	int freq_ncalc_pos = 100;
	int nstep_calc_pos = 0;
	int flag_gr_3 = 0;
	int ncalc_gr_3 = 0;
	int Period_calc_pos;
	double Max_thr = 0.5;
	double PosCur;
	double PosAvTpNew = 0.0;
	double PosAvTpOld = 0.0;
	int first_time = 1;
	double UAvTp;
	double Vgr = 6.55;

	// double dx_correction = 0.0;

	int kx0 = 100;

// ================================ position of bi-soliton ====================	
// ================================ damping in x-space ========================
	double Dcoeff = 0.005;
	fftw_plan cx_2_ckCur;
// ================================ damping in x-space ========================	
	
	// double EtaMax_duri, EtaMin;
	

	double tau = 0.01;
	double Time = 0.0;
	double FinTime = 110000.0;
	clock_t tCur;
//============================ const ===================	
	// double Ucr = 7.0;
	// double Umax = 0.0;
	// double Dcoeff = 400.0;
	// double alpha_k2 = 0.75;
	


	double print_Time_0;// = 0.0;
	double print_Time_f;// = 835.0;// INITIAL!!!


	// int print_Restart = 2500;
	int print_Restart = 100000;
	int print_Enrg = 100000; // calc_and_write_amplification
	int print_Ev; // = print_Restart

	int nstep_Ev = 0;
	int nstep_Restart = 0;
	// int print_Ev = 100;
	int k0 = 100;
	int dk0 = 1;
	double lNLS = 0.002;
	fftw_complex R0 = 0.1;
	fftw_complex V0;
	double Length;
	double grav = 9.81;
	double Dk;
	double norm;
	double du;
	double Epot;
	double Ekin;
	double Etot;
	double Momx;
	double Momy;
	int nstep = 0;
	int contin = -1;
//============================ const ===================
	double Diff[N];

// ========================= Distribution ==============
	double EtaMax, EtaMin;
	double EpotMax, EpotMin;
	double EkinMax, EkinMin;
	double MomxMax, MomxMin;
	double MomyMax, MomyMin;
	double dEtaMax, dEtaMin;

	double EtaMaxTot, EtaMinTot;
	double EpotMaxTot, EpotMinTot;
	double EkinMaxTot, EkinMinTot;
	double MomxMaxTot, MomxMinTot;
	double MomyMaxTot, MomyMinTot;
	double dEtaMaxTot, dEtaMinTot;
	
	double Epotu[N];
	double Ekinu[N];
	double Momxu[N];
	double Momyu[N];
	fftw_complex HdPsik[N];
	fftw_complex HdPsiu[N];
	fftw_plan HdPsik_2_HdPsiu_;
// ========================= Distribution ==============

	fftw_complex Etak[N];
	fftw_complex PsikIni[N];
	fftw_plan Psiu_2_Psik;
	fftw_plan dPhik_2_dPhiu_;


	fftw_complex RkCur[N];
	fftw_complex RkOld[N];
	fftw_complex RkNew[N];
	fftw_complex RuCur[N];
	fftw_complex VkCur[N];
	fftw_complex VkOld[N];
	fftw_complex VkNew[N];
	fftw_complex VuCur[N];

	fftw_complex dPhiu[N];
	fftw_complex dPhik[N];
	fftw_complex Phiu[N];
	fftw_complex Phik[N];

	fftw_complex dPsiu[N];
	fftw_complex dPsik[N];
	fftw_complex Psiu[N];
	fftw_complex Psik[N];

	fftw_plan dPhiu_2_dPhik;
	fftw_plan Phik_2_Phiu_;
	fftw_plan dPsiu_2_dPsik;
	fftw_plan Psik_2_Psiu_;
	fftw_plan dZk_2_dZu_;
	fftw_plan Yu_2_Yk;
	fftw_plan dXTk_2_dXTu_;
	fftw_plan XTk_2_XTu_;
	fftw_plan Yk_2_Yu_;
	fftw_plan RkCur_2_RuCur_;
	fftw_plan VkCur_2_VuCur_;
	fftw_plan RuCur_2_RkCur;
	fftw_plan VuCur_2_VkCur;
	fftw_plan dRk_2_dRu_;
	fftw_plan dVk_2_dVu_;
	fftw_plan dVVcu_2_dVVck;
	fftw_plan dVVck_2_dVVcu_;
	fftw_plan Uu_2_Uk;
	fftw_plan Uk_2_Uu_;
	fftw_plan dUk_2_dUu_;
	fftw_plan rhsRu_2_rhsRk;
	fftw_plan rhsVu_2_rhsVk;


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
	// fftw_plan dXTu_2_dXTk;
	// fftw_plan dYu_2_dYk;
//=======================================================	
	fftw_plan dYk_2_dYu_;
	fftw_complex Steepness_u[N];
	
//=======================================================

	
	
	fftw_complex dRu[N];
	fftw_complex dRk[N];
	fftw_complex dVu[N];
	fftw_complex dVk[N];
	

	fftw_complex Uu[N];
	fftw_complex Uk[N];
	fftw_complex dUu[N];
	fftw_complex dUk[N];
	fftw_complex dVVck[N];
	fftw_complex dVVcu[N];

	
	
	

	fftw_complex rhsRu[N];
	fftw_complex rhsRk[N];
	fftw_complex rhsVu[N];
	fftw_complex rhsVk[N];



	fftw_complex rhsR1[N];
	fftw_complex rhsR2[N];
	fftw_complex rhsR3[N];
	fftw_complex rhsR4[N];
	fftw_complex rhsV1[N];
	fftw_complex rhsV2[N];
	fftw_complex rhsV3[N];
	fftw_complex rhsV4[N];

	FILE *energy;
	FILE *pfileR, *pfileZ, *pfileV, *pfileEnMom;
	FILE *surf, *psi_data,  *pfilePhi, *temp;
	FILE *data, *restart, *pfileU;

	FILE *data_Ru;
	FILE *data_Vu;

