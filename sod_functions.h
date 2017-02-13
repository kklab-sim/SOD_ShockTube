#include<stdio.h>
#include<math.h>
#include<stdarg.h>
#include<stdlib.h>
#include<string.h>
#include<dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#define getName(var)  #var

// MINMOD limiter function
double minmod(double r)
{
	return fmax(0.0,fmin(1.0,r));
}

// SUPERBEE limiter function
double superbee(double r)
{
	return fmax(fmax(0.0,fmin(2*r,1.0)),fmax(0.0,fmin(2.0,r)));
}

// VAN LEER limiter function
double vanLeer(double r)
{
	return (r + fabs(r))/(1+fabs(r));
}

// VAN ALBADA limiter function
double vanAlbada(double r)
{
	return (2*r)/(1+r*r);
}

// MUSCL operator: xleft and xright are the outputs, and xm1,x,xp1,xp2 correspond to
// the (i-1,i,i+1,i+2) indexes when calculating at face (i+1/2), and to the indexes
// (i-2,i-1,i,i+1) on the (i-1/2) face 
void MUSCL(double *xleft,double *xright,double xm1,double x,double xp1,double xp2,
double muscl_bias, char* limiter)
{
    double deltaplusL,deltaminL,deltaplusR,deltaminR;
    double rL,rR;

    // Depending on MUSCL_BIAS:
    // double BIAS = (3 - muscl_bias)/(1 - muscl_bias);

    // Choice of flux limiter function
    if (strcmp(limiter,"none")==0)
    {
	    deltaminL = x - xm1;
	    deltaplusL = xp1 - x;

	    deltaminR = xp1 - x;
	    deltaplusR = xp2 - xp1;

	    *xleft = x + 0.25*((1 - muscl_bias)*deltaminL + (1 + muscl_bias)*deltaplusL);
	    *xright = xp1 - 0.25*((1 - muscl_bias)*deltaplusR + (1 + muscl_bias)*deltaminR);

	    return;
    }
    else if (strcmp(limiter,"minmod")==0)
    {	
	    deltaminL = x - xm1;
	    deltaplusL = xp1 - x;
	    rL = deltaminL/deltaplusL;

	    deltaminR = xp1 - x;
	    deltaplusR = xp2 - xp1;
	    rR = deltaminR/deltaplusR;

	    *xleft = x + 0.25*minmod(rL)*((1 - muscl_bias)*deltaminL + (1 + muscl_bias)*deltaplusL);
	    *xright = xp1 - 0.25*minmod(1.0/rR)*((1 - muscl_bias)*deltaplusR + (1 + muscl_bias)*deltaminR);

	    return;
    } 
    else if (strcmp(limiter,"superbee")==0)
    {	
	    deltaminL = x - xm1;
	    deltaplusL = xp1 - x;
	    rL = deltaminL/deltaplusL;

	    deltaminR = xp1 - x;
	    deltaplusR = xp2 - xp1;
	    rR = deltaminR/deltaplusR;

	    *xleft = x + 0.25*superbee(rL)*((1 - muscl_bias)*deltaminL + (1 + muscl_bias)*deltaplusL);
	    *xright = xp1 - 0.25*superbee(rR)*((1 - muscl_bias)*deltaplusR + (1 + muscl_bias)*deltaminR);

	    return;
    } 
    else if (strcmp(limiter,"vanLeer")==0)
    {	
	    deltaminL = x - xm1;
	    deltaplusL = xp1 - x;
	    rL = deltaminL/deltaplusL;

	    deltaminR = xp1 - x;
	    deltaplusR = xp2 - xp1;
	    rR = deltaminR/deltaplusR;

	    *xleft = x + 0.25*vanLeer(rL)*((1 - muscl_bias)*deltaminL + (1 + muscl_bias)*deltaplusL);
	    *xright = xp1 - 0.25*vanLeer(rR)*((1 - muscl_bias)*deltaplusR + (1 + muscl_bias)*deltaminR);

	    return;
    }
    else if (strcmp(limiter,"vanAlbada")==0)
    {	
	    deltaminL = x - xm1;
	    deltaplusL = xp1 - x;
	    rL = deltaminL/deltaplusL;

	    deltaminR = xp1 - x;
	    deltaplusR = xp2 - xp1;
	    rR = deltaminR/deltaplusR;

	    *xleft = x + 0.25*vanAlbada(1.0/rL)*((1 - muscl_bias)*deltaminL + (1 + muscl_bias)*deltaplusL);
	    *xright = xp1 - 0.25*vanAlbada(1.0/rR)*((1 - muscl_bias)*deltaplusR + (1 + muscl_bias)*deltaminR);

	    return;
    } 
}

// AUSM Flux-Splitting (DV scheme including entropy and shock fixes)
void AUSM(double *flux,double rhoL, double rhoR,double uL,double uR,double pL,
double pR,double hL,double hR,double gamma, double bias,double ENTRO_BIAS)
{
	// Variable definition
    double s=0.0;
    double rhou2MV=0.0,rhou2MD=0.0,rhouhalf=0.0,uhalf=0.0,phalf=0.0;
    double uLplus=0.0,uRminus=0.0,pLplus=0.0,pRminus=0.0;

    // Left and right parameters of AUSM scheme
    double cL=sqrt(gamma*pL/rhoL),cR=sqrt(gamma*pR/rhoR);
    double alphaL=pL/rhoL,alphaR=pR/rhoR;
    double ALPHA_L=2*alphaL/(alphaL+alphaR), ALPHA_R=2*alphaR/(alphaR+alphaL);
    double cM=fmax(cL,cR);

    // Definition of UL(plus), PL(plus), UR(minus) and PR(minus)
    if (fabs(uL)<=cM)
    {
        uLplus=0.5*(uL+fabs(uL))+ALPHA_L*(0.25/cM*pow(uL+cM,2.0)-0.5*(uL+fabs(uL)));
        pLplus=0.25*pL*pow(uL/cM+1,2.0)*(2-uL/cM);
    }
    else
    {
        uLplus=0.5*(uL+fabs(uL));
        pLplus=0.5*pL*(1+fabs(uL)/uL);
    }
    if (fabs(uR)<=cM)
    {
        uRminus=0.5*(uR-fabs(uR))+ALPHA_R*(-0.25/cM*pow(uR-cM,2.0)-0.5*(uR-fabs(uR)));
        pRminus=0.25*pR*pow(uR/cM-1,2.0)*(2+uR/cM);
    }
    else
    {
        uRminus=0.5*(uR-fabs(uR));
        pRminus=0.5*pR*(1-fabs(uR)/uR);
    }

    // Bias function for pressure gradient
    s = 0.5*fmin(1.0,bias*fabs(pR-pL)/fmin(pR,pL));

    // AUSM-V and AUSM-D momentum flux terms
    uhalf=uLplus+uRminus;
    phalf=pLplus+pRminus;
    rhouhalf=0.5*(uhalf*(rhoR+rhoL)-fabs(uhalf)*(rhoR-rhoL));
    rhou2MD=0.5*(rhouhalf*(uL+uR)-fabs(rhouhalf)*(uR-uL));
    rhou2MV=uLplus*rhoL*uL+uRminus*rhoR*uR;

    // AUSM-DV expression of flux
    flux[0] = uLplus*rhoL+uRminus*rhoR;
    flux[1] = (0.5+s)*rhou2MV+(0.5-s)*rhou2MD+phalf;
    flux[2] = 0.5*(rhouhalf*(hL+hR)-fabs(rhouhalf)*(hR-hL));

    // Entropy fix (numerical dissipation for single expansion waves)
    if((uL-cL<0.0 && uR-cR>0.0)&&(!(uL+cL<0.0 && uR+cR>0.0)))
    {
        flux[0]=flux[0]-ENTRO_BIAS*(uR-cR-uL+cL)*(rhoR-rhoL);
        flux[1]=flux[1]-ENTRO_BIAS*(uR-cR-uL+cL)*(rhoR*uR-rhoL*uL);
        flux[2]=flux[2]-ENTRO_BIAS*(uR-cR-uL+cL)*(rhoR*hR-rhoL*hL);
    }
    else if((!(uL-cL<0.0 && uR-cR>0.0))&&(uL+cL<0.0 && uR+cR>0.0))
    {
        flux[0]=flux[0]-ENTRO_BIAS*(uR+cR-uL-cL)*(rhoR-rhoL);
        flux[1]=flux[1]-ENTRO_BIAS*(uR+cR-uL-cL)*(rhoR*uR-rhoL*uL);
        flux[2]=flux[2]-ENTRO_BIAS*(uR+cR-uL-cL)*(rhoR*hR-rhoL*hL);
    }
    return;
}