#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>
#include<unistd.h>  
#include <sys/stat.h>
#include <sys/types.h>

#include "MT19937-64.c" //Random numerbs
#define PI 3.1415926
#define size 1024 //System size
#define sample 0.1


#include "Modular_network_generator.c"
#include "normal_rand.c"
//double K = 0.1;
//double r = 1e-0;
double C; //Noise level
int nm=32, av_k=28; // number of modules and average degree
double const dt = 1e-2; //tmie step
double const t_f = 20, t_i = 0;
double n_steps = (t_f-t_i)/dt; ///number of steps
double V[size][size]; //adjacency matrix
double g_m, m_m;
void network_generator(const int N, int nm, int av_k, double r, double V[N][N]);
double std_norm();
void assign_freq(const int N, double f[N], int mode);
void initialize(const int N, double ph[N]);
void x_dot(const int N, double V[N][N], double ph[N], double ph_r[N], double f[N], double K);
void solver(const int N, double ph[N], double ph_r[N], double f[N], double V[N][N], double r, double K, double E);
double mag(const int N, double ph[N]);
double mod_mag(const int N, int nm, double ph[N]);


int main(int argc, char *argv[])
{
	const int N = size; // system size
	double phase[N], phase_rate[N], freq[N]; // phase variables and intrinsic frequencies
	double r, E_i, K = 100; 
	const int n_nets = 3, n_init = 4; // no of networks and initial conditions
	const int n_iter = n_nets*n_init; 
	double avg_g[n_iter], avg_m[n_iter];
	int seed = (unsigned) time(NULL);
	init_genrand64(seed);
    //double r_l[1] = {1.0};
    double base, power;
	//sscanf(argv[1],"%lf",&r);
	//sscanf(argv[2],"%lf",&E_i);
	r = atof(argv[1]);
	E_i = atof(argv[2]);
	printf("r = %.4f, E = %.4f\n", r, E_i);    
	double E_f = E_i+0.1, dE = 1.0;
	//for (int r_i = 0; r_i<1; r_i++)
	//{    
	    //r = r_l[r_i];
	    //printf("r : %.4f \n", r);
	    char fname[100];
	    sprintf(fname, "N_%d_nm_%d_r_%.4f_E_%.4f_iter_%d.dat", N, nm, r, E_i, n_iter);
	    FILE *fp;
	    fp = fopen(fname, "w");
	    for (double E = E_i; E<=E_f; E = E+dE)
	    {
        	//printf("E : %.4f \n", E);
        	int idx=0;
	        for (int i = 0; i<n_nets; i++)
	        {
		        network_generator(N, nm, av_k, r, V);	
		        for (int j=0; j<n_init; j++)
		        {
                    g_m = 1.5, m_m = 1.5;
                    assign_freq(N, freq, 0);
			        initialize(N, phase);
			        printf("idx : %d \n", idx);
			        solver(N, phase, phase_rate, freq, V, r, K, E);
			        avg_g[idx] = g_m;
			        avg_m[idx] = m_m;
			        //printf("cp2; g_m : %.4f", g_m);
                	//printf("cp2; m_m : %.4f", m_m);
			        idx = idx+1; 
		        }
	        }
  	        //Statistics
            double G = 0, M = 0, G2 = 0, M2 = 0, var_G = 0, var_M = 0;
  	        for (int k = 0; k<n_iter; k++)
  	        {
                G = avg_g[k]/n_iter + G;
                M = avg_m[k]/n_iter + M;
                G2 = avg_g[k]*avg_g[k]/n_iter + G2;
                M2 = avg_m[k]*avg_m[k]/n_iter + M2;
  	        }
            var_G = G2 - G*G;
            var_M = M2 - M*M;  
            fprintf(fp, "%.4e %.12e %.12e %.12e %.12e \n", E, G, var_G, M, var_M);          
        }
        fclose(fp);   
    //}
	
	/*
	char fname[100];
	FILE *f1;
	sprintf(fname, "Time_r_%.4f_k_%d_N_%d_nm_%d_K_%.4f_t_%.2f_dt_%.4f.dat", r, av_k, (int)size, nm, K, (t_f-t_i), dt);
	f1 = fopen(fname, "w");
	double n_s = (t_f-t_i)/dt, t=0;
	for (int i=0; i<n_s; i++)
	{
		fprintf(f1, "%.4f\n", t);
		t = t+dt;
	}
	fclose(f1);
	*/
	printf("C : \n");
	printf("M_ar : \n");
	return 0;
}

void assign_freq(const int N, double f[N], int mode) // assign zero intrinsic frequncy
{
	if (mode == 0)
		for (int i=0; i<N; i++)
			f[i]=0;
	else
		printf("Enter appropriate mode\n");
}

void initialize(const int N, double ph[N]) //initial conditons
{
	for (int i=0; i<N; i++)
		ph[i] = genrand64_real2()*2.0*PI;
}

void x_dot(const int N, double V[N][N], double ph[N], double ph_r[N], double f[N], double K) //compute rate of change
{
	double tmp;
	for (int i=0; i<N; i++)
	{
		tmp = 0;
		for (int j=0; j<N; j++)
		{
				tmp = tmp + V[i][j]*(sin(ph[j]-ph[i]));
		}
		ph_r[i] = f[i] + K*tmp/av_k;
	}
}


/*
void MP_Euler_step(double *x, double *x_new, double dt, double K)
{
	double x_tmp[size], x_r[size];
	x_dot(x, x_r, freq, K);
	for (int i=0; i<size; i++)
	{
		x_tmp[i] = x[i] + x_r[i]*dt*0.5; //Half step
	}
	//printf("%f", x_r[1]);
	x_dot(x_tmp, x_r, freq, K);
	for (int i=0; i<size; i++)
	{
		x_new[i] = x[i] + x_r[i]*dt;
	}
	//printf("%f", x_r[1]);
}
*/


void solver(const int N, double ph[N], double ph_r[N], double f[N], double V[N][N], double r, double K, double E) // Stochastic solver
{
	int n = N/nm;
	double phase_tmp[N], t;
	double m_steps = sample*n_steps;
	double time = t_f-t_i;
    double noise_dt = sqrt(dt);
	t=0;
	char fn1[100], fn2[100];
	//sprintf(fn1, "GM_t_series_r_%.4f_E_%.4f.dat", r, E);
	//sprintf(fn2, "MM_t_series_r_%.4f_E_%.4f.dat", r, E);
	//FILE *fp1, *fp2;
	//fp1 = fopen(fn1, "w");
	//fp2 = fopen(fn2, "w");
	for (int i=0; i<(int)(n_steps-m_steps); i++)
	{
		x_dot(N, V, ph, ph_r, f, K);
		for (int j=0; j<N; j++)
		{
		    phase_tmp[j] = ph[j] + ph_r[j]*dt + E*std_norm()*noise_dt;
		}
		
		for (int j=0; j<N; j++)
		{
			ph[j] = phase_tmp[j];
		}
		//fprintf(fp1, "%.4f %.4f \n", t, mag(N));
     	//fprintf(fp2, "%.4f %.4f \n", t, mod_mag(N, 16));
		t=t+dt;
	}
    g_m = 0; m_m = 0;	
	for (int i=0; i<(int)(m_steps); i++)
	{
		g_m = g_m + mag(N, ph);
		m_m = m_m + mod_mag(N, nm, ph);		
		x_dot(N, V, ph, ph_r, f, K);
		for (int j=0; j<N; j++)
		{
		    phase_tmp[j] = ph[j] + ph_r[j]*dt + E*std_norm()*noise_dt;
		}
		
		for (int j=0; j<N; j++)
		{
			ph[j] = phase_tmp[j];
		}
		//fprintf(fp1, "%.4f %.4f \n", t, mag(N));
     	//fprintf(fp2, "%.4f %.4f \n", t, mod_mag(N, 16));
		t=t+dt;
	}
	g_m = g_m/m_steps;
	m_m = m_m/m_steps;
	//printf("cp1; avg_g : %.4f", g_m);
	//printf("cp1; avg_m : %.4f", m_m);	
	//fclose(fp1);
	//fclose(fp2);
}

double mag(const int N, double ph[N])
{	
	double complex tmp=0;
	for (int i=0; i<N; i++)
	{
		tmp = tmp + cexp(I*ph[i]);
	}
	double avg = cabs(tmp)/N;
	return (avg);
}

double mod_mag(const int N, int nm, double ph[N])
{
	double mod_size = (double)N/(double)nm;
	int initial=0, final=mod_size;
	double avg=0;
	double complex tmp;
	for (int i=0; i<nm; i++)
	{	
		tmp=0;
		for (int j=initial; j<final; j++)
		{
			tmp = tmp + cexp(I*ph[j]);
		}
		avg = avg + cabs(tmp);
		initial=initial+mod_size;
		final = final+mod_size;
	}
	avg=avg/N;
	return (avg);
}
