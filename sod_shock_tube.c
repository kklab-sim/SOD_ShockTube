/* 1D SOD SHOCK TUBE CASE

   AUTHOR: Florian NGUYEN - The University of Tokyo, Jan. 2017.

   DESCRIPTION: 1D code solving Euler equations in the case of
   Sod's shock tube case, as defined by:
   P1=1   ; RHO1=1     ; U1=0
   P2=0.1 ; RHO2=0.125 ; U2=0

   NUMERICAL METHOD: 
   AUSM-DV (Advection Upstream Splitting Method) Flux Splitting
   MUSCL interpolation at cell interfaces for higher order accuracy (with flux limiters)
   Backward Euler time differencing scheme (1st order accuracy)

*/


// HEADER
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>

#include "sod_functions.h"

// DEFINITION OF FUNCTIONS
extern double minmod(double r);
extern double superbee(double r);
extern double vanLeer(double r);
extern double vanAlbada(double r);
extern void MUSCL(double *xleft,double *xright,double xm1,
	double x,double xp1,double xp2,double muscl_bias, char* limiter);
extern void AUSM(double *flux,double rhoL, double rhoR,
	double uL,double uR,double pL,double pR,double hL,
	double hR,double gamma, double bias,double ENTRO_BIAS);


// BEGINNING OF MAIN ROUTINE
int main()
{
	// User parameters
	double L = 1.0;                        // Length of computational domain (m)
	int N = 200;                           // Number of cells (should be even number)
	double T_SIM = 0.2;                    // Simulation time (sec)
	double DT = 1e-4;                      // Time step (sec)
	int USE_MUSCL = 1;                     // Space accuracy: 1st order if 0, 2nd order if 1
	double X0 = 0.0;                       // Absciss of the left boundary of the calculation domain (m)
	double GAMMA = 1.4;                    // Heat capacity ratio
	double R = 287.0;                      // Specific gas constant (J/kg/K)
	char* fileName = "SOD_RESULTS.DAT";    // Name of the file in which to save results
	double MUSCL_BIAS = 1.0/3;             // MUSCL bias for 3rd order space accuracy
	char* MUSCL_LIMITER = "minmod";        // Limiter function for gradients: none,minmod,superbee,vanLeer,vanAlbada
	double AUSM_BIAS = 10.0;               // Bias coef. for pressure gradient in AUSM
	double AUSM_ENTROPY_FIX = 0.125;       // Parameter for the entropy fix in AUSM


	// Deduced parameters
	double DX = L/N;                       // Grid size (m)
	double X_SHOCK = X0 + L/2;             // Position of discontinuity (m)
	int N_STEP_SIM = 1 + (int)(T_SIM/DT);  // Number of time iterations required


	// Print parameters for checking
	printf("Parameters of the current simulation:");
	printf("\nT_SIM = %f sec.",T_SIM);
	printf("\nDT = %f ms.",DT*1000);
	printf("\nDX = %f m.",DX);
	printf("\nN = %d cells over the domain.",N);
	printf("\nSteps to perform: %d",N_STEP_SIM);


	// Cell-centred mesh creation (including ghost cells)
	int N_GHOST = USE_MUSCL + 1;           // Number of ghost cells at each boundary
	int N_CELL = N + 2*N_GHOST;            // Total number of cells over domain (including ghost cells)
	double x_cell[N_CELL];                 // Cell-centre coordinates
	double x_node[N_CELL + 1];             // Cell-boundary coordinates

	for (int i = 0 ; i < N_CELL ; ++i)
	{
		x_cell[i] = X0 + (i - N_GHOST)*DX + DX/2 ;
		x_node[i] = X0 + (i - N_GHOST)*DX;
	}
	x_node[N_CELL] = x_node[N_CELL - 1]+DX;

	printf("\n\nMesh computed...");


	// Generate initial conditions on main domain
	double rho[N_CELL], p[N_CELL], u[N_CELL], T[N_CELL], e[N_CELL], h[N_CELL];
	for (int i = N_GHOST; i < N_CELL - N_GHOST; ++i)
	{
		if (x_cell[i] <= X_SHOCK)
		{
			// Primitive variables
			rho[i] = 1.0;
			u[i] = 0.0;
			p[i] = 1.0;

			// Redundant variables
			T[i] = p[i]/(rho[i]*R);
			e[i] = p[i]/(GAMMA - 1) + 0.5*rho[i]*pow(u[i],2.0);
			h[i] = (e[i] + p[i])/rho[i];

		}
		else if (x_cell[i] > X_SHOCK)
		{
			// Primitive variables
			rho[i] = 0.125;
			u[i] = 0.0;
			p[i] = 0.1;

			// Redundant variables
			T[i] = p[i]/(rho[i]*R);
			e[i] = p[i]/(GAMMA - 1) + 0.5*rho[i]*pow(u[i],2.0);
			h[i] = (e[i] + p[i])/rho[i];
		}
	}


	// Variables initialization before time loop
	double flux_left[3],flux_right[3];
	double rhol,rhor,ul,ur,pl,pr,hl,hr,el,er;
	double rhoN[N_CELL], pN[N_CELL], uN[N_CELL], TN[N_CELL], eN[N_CELL], hN[N_CELL];


	// TIME LOOP STARTING...
	for (int nt = 0 ; nt < N_STEP_SIM ; ++nt)
	{
		// BOUNDARY CONDITIONS AT LEFT AND RIGHT BOUNDARIES (MIRROR WALLS)
		for (int i = N_GHOST ; i < N_CELL - N_GHOST ; ++i)
		{
			for (int i = 0; i < N_GHOST ; ++i)
			{
				// Left boundary (High pressure and density): 0 equals 3 and 1 equals 2
				rho[N_GHOST - 1 - i] = rho[N_GHOST + i];
				u[N_GHOST - 1 - i] = -u[N_GHOST + i];
				p[N_GHOST - 1 - i] = p[N_GHOST + i];

				T[N_GHOST - 1 - i] = p[N_GHOST - 1 - i]/(rho[N_GHOST - 1 - i]*R);
				e[N_GHOST - 1 - i] = p[N_GHOST - 1 - i]/(GAMMA - 1) + 0.5*rho[N_GHOST - 1 - i]*pow(u[N_GHOST - 1 - i],2.0);

				// Right boundary (Low pressure and density): N-2 equals N-3 and N-1 equals N-4
				rho[N_CELL - N_GHOST + i] = rho[N_CELL - N_GHOST-1-i];
				u[N_CELL - N_GHOST + i] = - u[N_CELL - N_GHOST-1-i];
				p[N_CELL - N_GHOST + i] = p[N_CELL - N_GHOST-1-i];

				T[N_CELL - N_GHOST + i] = p[N_CELL - N_GHOST + i]/(rho[N_CELL - N_GHOST + i]*R);
				e[N_CELL - N_GHOST + i] = p[N_CELL - N_GHOST + i]/(GAMMA-1) + 0.5*rho[N_CELL - N_GHOST + i]*pow(u[N_CELL - N_GHOST + i],2.0);
			}
		}


		// FLUX SPLITTING (BASED ON AUSM-DV AND HANEL'S DISSIPATIVE SCHEME)
		for (int i = N_GHOST; i < N_CELL - N_GHOST ; ++i)
		{
			// Left cell interface parameters (i - 1/2)
			if (USE_MUSCL)
			{
				// Interpolation of primitive variables (rho,u,p)
				MUSCL(&rhol,&rhor,rho[i-2],rho[i-1],rho[i],rho[i+1],MUSCL_BIAS,MUSCL_LIMITER);
				MUSCL(&ul,&ur,u[i-2],u[i-1],u[i],u[i+1],MUSCL_BIAS,MUSCL_LIMITER);
				MUSCL(&pl,&pr,p[i-2],p[i-1],p[i],p[i+1],MUSCL_BIAS,MUSCL_LIMITER);

				// Deduction of other parameters
				el = pl/(GAMMA-1) + 0.5*rhol*pow(ul,2.0); // Only needed to calculate hl
				er = pr/(GAMMA-1) + 0.5*rhor*pow(ur,2.0); // Only needed to calculate hr
				hl = (pl + el)/rhol;
				hr = (pr + er)/rhor;
			}
			else if (!USE_MUSCL)
			{
				rhol = rho[i-1];
				rhor = rho[i];
				ul = u[i-1];
				ur = u[i];
				pl = p[i-1];
				pr = p[i];
				hl = h[i-1];
				hr = h[i];
			}

			// Splitting flux on LEFT face
			AUSM(&flux_left[0],rhol,rhor,ul,ur,pl,pr,hl,hr,GAMMA,AUSM_BIAS,AUSM_ENTROPY_FIX);

			// Right cell interface parameters (i + 1/2)
			if (USE_MUSCL)
			{
				// Interpolation of primitive variables (rho,u,p)
				MUSCL(&rhol,&rhor,rho[i-1],rho[i],rho[i+1],rho[i+2],MUSCL_BIAS,MUSCL_LIMITER);
				MUSCL(&ul,&ur,u[i-1],u[i],u[i+1],u[i+2],MUSCL_BIAS,MUSCL_LIMITER);
				MUSCL(&pl,&pr,p[i-1],p[i],p[i+1],p[i+2],MUSCL_BIAS,MUSCL_LIMITER);

				// Deduction of other parameters
				el = pl/(GAMMA-1) + 0.5*rhol*pow(ul,2.0); // Only needed to calculate hl
				er = pr/(GAMMA-1) + 0.5*rhor*pow(ur,2.0); // Only needed to calculate hr
				hl = (pl + el)/rhol;
				hr = (pr + er)/rhor;
			}
			else if (!USE_MUSCL)
			{
				rhol = rho[i];
				rhor = rho[i+1];
				ul = u[i];
				ur = u[i+1];
				pl = p[i];
				pr = p[i+1];
				hl = h[i];
				hr = h[i+1];
			}

			// Splitting flux on LEFT face
			AUSM(&flux_right[0],rhol,rhor,ul,ur,pl,pr,hl,hr,GAMMA,AUSM_BIAS,AUSM_ENTROPY_FIX);


			// COMPUTATION OF SOLUTION
			rhoN[i] = rho[i] - DT/DX*(flux_right[0] - flux_left[0]);
			uN[i] = 1.0/rhoN[i]*(rho[i]*u[i] - DT/DX*(flux_right[1] - flux_left[1]));
			eN[i] = e[i] - DT/DX*(flux_right[2] - flux_left[2]);


			// UPDATE OF SOLUTION BEFORE NEXT ITERATION
			pN[i] = (GAMMA-1)*(eN[i]-0.5*rhoN[i]*pow(uN[i],2.0));
			TN[i] = pN[i]/(rhoN[i]*R);
			hN[i] = (eN[i] + pN[i])/rhoN[i];
		}


		// Update all parameters before the next iteration
		for (int i = N_GHOST ; i < N_CELL - N_GHOST ; ++i)
		{
			rho[i] = rhoN[i];
			u[i] = uN[i];
			p[i] = pN[i];
			e[i] = eN[i];
			T[i] = TN[i];
			h[i] = hN[i];
		}
	} // END OF TIME LOOP	


	// Print a message after the last iteration
	printf("\n\nTime loop finished: exporting results...");


	// EXPORTING ALL RESULTS IN A SINGLE .DAT FILE
	FILE *fp=fopen(fileName,"w");
    if(fp == NULL)
    {
        sprintf("ERROR: File '",fileName,"' could not be opened for writing!");
    }
    else
    {
        for (int i = N_GHOST ; i < N_CELL - N_GHOST ; ++i)
        {
            fprintf(fp,"%f %f %f %.8f %.8f %.8f %.8f %.8f %.8f\n",x_node[i],x_node[i+1],x_cell[i],rho[i],u[i],p[i],T[i],e[i],h[i]);
        }
    }
    fclose(fp);

    printf("\n\nRoutine ended successfully! -Florian\n");
}

// END OF FILE
