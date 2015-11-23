//	QUESTIONS!
//	Should I call this function each time a new model is estimated? I.e. in each repition? Or better once with iN*iS rows?
//	We need to drop the first p (i.e. 1 or 2) observations. Should I then use iN+1 bzw. iN+2 observations in the first place?
//	Where and when do I set the common root? Before each estimation technique, the root should be reset. (i.e 4 times)
//		->	ranseed(-1) resets to the initial seed. 
//	How will I do this for the two different cases?

// 	ASSUMPTIONS!
//	The DGP assumes that errors (epsilon) are normally distributed. Assumes that all the parameters are as given in the input.
//	Assumes that when calculating vY with vY_{t-1}, setting first X observation to 0 is fine. Assumes AR(1) is implemented
//	correctly by cumsum (and that's what is meant by it). 	

/*
**	simulation(const iN, const vb0, const vrho0, const vs0, const vp0 , const amX, const avY)
**
**  Purpose:
**    	Simulate the data, to be used in the estimation
**		Watch out!!! I only simulate iN observations. So this function needs to be called in the loops each time again.
**				-> Alternatively I could make iN * iS observations and then call different indices in the loops. Whats better?
**
**  Inputs:
**		iN 		# of observations
**		vb0 	DGP parameters for beta
**		vrho0 	DGP parameter for rho (depending on the case input is either vrho01 or vrho02)
**		vs0		DGP variance of disturbances (epsilon)
**		vp0  	DGP parameter for AR(1) process of X's 
**
**  Output:
**	  	mX 		(iN x 4),	matrix of exogenous regressors
**		vY 		(iN x 1),	vector of simulated y-values (see (21) "D&McK").
**
**  Return value:
**	  	br 		returns 0 if any simulated value is NaN
*/
simulation(const iN, const vb0, const vrho0, const vs0, const vp0 , const amX, const avY)
{
	decl vdisturb, iK, iAdd, mX_noconst, vY_nolag, br;

	iK = rows(vb0);				//	# of columns of X						
	iAdd = sizerc(vrho0)+2;		//	# of obs. need to be added to have iN after deleting first rho. 
	vdisturb = sqrt(vs0) .* rann(iN+iAdd, iK-1);	

	mX_noconst = cumsum(vdisturb, vp0);
	amX[0] = 1~mX_noconst;

	vY_nolag = amX[0]*vb0 - vrho0*lag(amX[0], 1)*vb0 + sqrt(vs0) .* rann(iN+iAdd, 1);		// What about this Lag0?  
	avY[0] = vY_nolag + vrho0 * lag(vY_nolag, 1);

	// amX[0] = amX[0][iAdd:][];
	amX[0] = (amX[0]~lag(avY[0], 1)~lag(amX[0][][1:], 1))[iAdd:][];			// Also what about this Lag0?
	avY[0] = avY[0][iAdd:][];

	br = !(isnan(amX[0]~avY[0]));
	return br;

}
