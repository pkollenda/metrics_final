/*
**	assignmain_v0.ox  
**
**  Purpose:
**		Comparing empirical rejection probabilities of Common Factor Tests in
**			nonlinear autoregressive model by testing based on F-Test, full
**			bootstrap and approximate bootstrap according to Davidson, MacKinnon
**			(1999).
** 
** 		v0	Set up structure, incorporate simulation of data 
**
**  Date:
**   	Sun Nov 22 20:47:13 2015
**
**  Author:
**    	Philipp Kollenda
**		MPhil Tinbergen Institute Amsterdam
**
*/
#include <oxstd.h>
//#include <oxfloat.h>   	// Includes predefined constant numbers
//#include <oxprob.h>	 	// Includes Ox probability package
//#include <oxdraw.h>	  	// Includes Ox graphing package
//#include <packages/gnudraw/gnudraw.h>
//#include <packages/oxutils/oxutils.h>


#include "./support_files/datasim_v0.ox"
#include "./support_files/ftest.ox"

main()
{
	decl iN, vb0, vrho01, vrho02, vs0, vp0, vpar01, vpar02;
	decl iS, iB, ima, imb, dalpha;
	decl mX, vY, vpar_unrest;
	decl br;

	// Magic Numbers
		// Setting the seed
		ranseed(136);			// Setting the seed for the random number generator

		// Data generating process
		iN = 20;				// 	# of observations
		vb0 = <1; 2; 3; 4>;		// 	DGP parameters for beta
		vrho01 = 0.9;			//	DGP parameter for rho (case 1)
		vrho02 = <-0.3; 0.1>;	//	DGP parameter for rho (case 2) 
		vs0 = 1;				//	DGP variance of disturbances (epsilon)
		vp0 = 0.5;				// 	DGP parameter for AR(1) process of X's 

		// Estimation parameters
		iS = 100;				//	# of repetitions
		iB = 99;				//	# of bootstrap samples
		ima = 2;				//	step size in approximate bootstrap (case a)
		imb = 3;				//	step size in approximate bootstrap (case b)
		dalpha = 0.05;			//	desired level of the test
		vpar01 = vb0|vrho01|-vrho01*vb0[1:];
		vpar02 = vb0|vrho02|-vrho02[0]*vb0|-vrho02[1]*vb0;
	
	// Initialisation
	br = simulation(iN, vb0, vrho01, vs0, vp0 , &mX, &vY);

	// Estimation
	br = ftest(mX, vY, &vpar_unrest) && br;
	
	// Show Output
		// Make Graphs

		// Debugging & Checks
		print(vpar_unrest~vpar01);
}
