/*---------------------------------------------------------------------------*/
// This file contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution through the nonlinear crystal.
// * 
// * In particular, this file should be used for any three-wave mixing problem, i.e.,
// * sum or difference frequency generation (SFG or DFG) for a single-pass 
// * (just one round trip) or inside a cavity (multiple round trips).
/*---------------------------------------------------------------------------*/


#ifndef _CWES3CUH
#define _CWES3CUH

#pragma once


/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x,y are different fields) */
// INPUTS 
// dAp, dAs, dAi : evolved electric fields
//  Ap,  As,  Ai: electric fields
//  lp,  ls,  li: wavelengths
//  kp,  ks,  ki: kappas
//        z: crystal position
//     SIZE: size vector
// OUTOPUT
// save in dAp, dAs, dAi the evolved electric fields
void dAdz( complex_t *dAp, complex_t *dAs,  complex_t *dAi, complex_t *Ap, complex_t *As, complex_t *Ai, real_t lp, real_t ls, real_t li, real_t kp, real_t ks, real_t ki, real_t dk, real_t z )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	
	for (uint idx = 0; idx < SIZE; idx++){
		dAp[idx]  = Im * kp * As[idx] * Ai[idx] * CpxExp(-dk*z) ;
		dAs[idx]  = Im * ks * Ap[idx] * CpxConj(Ai[idx]) * CpxExp(+dk*z);
		dAi[idx]  = Im * ki * Ap[idx] * CpxConj(As[idx]) * CpxExp(+dk*z);
	}
	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
// INPUTS 
// auxp, auxs, auxi : auxiliary vectors
//   Ap,   As,   Ai : electric fields
//   lp,   ls,   li : wavelengths
//   kp,   ks,   ki : vectors for Runge-Kutta
//                 s: scalar for Runge-Kutta
//              SIZE: size vector
// OUTOPUT
// save in auxp, auxs, auxi the evolved electric fields
void LinealCombination( complex_t *auxp, complex_t *auxs, complex_t *auxi, complex_t *Ap, complex_t *As, complex_t *Ai, complex_t *kp, complex_t *ks, complex_t *ki, double s )
{
	
	for (uint idx = 0; idx < SIZE; idx++){
		auxp[idx] = Ap[idx] + kp[idx] * s;
		auxs[idx] = As[idx] + ks[idx] * s;
		auxi[idx] = Ai[idx] + ki[idx] * s;
	}
	
	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
// INPUTS 
//  Ap,  As, Ai : electric fields
//  kXp,kXs, kXi: vectors for Runge-Kutta
//       dz: step size
//     SIZE: size vector
// OUTOPUT
// Update electric fields using Runge-Kutta after one step size, dz
void rk4(complex_t *Ap, complex_t *As,  complex_t *Ai, complex_t *k1p, complex_t *k1s, complex_t *k1i, complex_t *k2p, complex_t *k2s, complex_t *k2i, complex_t *k3p, complex_t *k3s, complex_t *k3i, complex_t *k4p, complex_t *k4s, complex_t *k4i, real_t dz )
{
	
	for (uint idx = 0; idx < SIZE; idx++){
		Ap[idx] = Ap[idx] + (k1p[idx] + 2*k2p[idx] + 2*k3p[idx] + k4p[idx]) * dz / 6;
		As[idx] = As[idx] + (k1s[idx] + 2*k2s[idx] + 2*k3s[idx] + k4s[idx]) * dz / 6;
		Ai[idx] = Ai[idx] + (k1i[idx] + 2*k2i[idx] + 2*k3i[idx] + k4i[idx]) * dz / 6;
	}
	
	return ;
}


/** Computes the linear part: Ax = Ax.exp(i.f(Ω)*z), where f(Ω) is a frequency dependant functions
 * including the group velocity and the group velocity dispersion parameters. */
// INPUTS 
//     auxp,  auxs,  auxi : auxiliary vectors
//      Apw,   Asw,   Aiw : electric fields in frequency domain
//                      w : angular frequency 
//       lp,   ls,     li : wavelengths
//       vp,   vs,     vi : group-velocity
//      b2p,  b2s,    bsi : group-velocity dispersion
// alphap, alphas, alphai : linear absorpion
//           SIZE: size vector
//              z: crystal position
// OUTOPUT
// save in auxp, auxs the electric fields after appling dispersion
void LinearOperator(complex_t *auxp, complex_t *auxs, complex_t *auxi, complex_t *Apw, complex_t* Asw, complex_t* Aiw, real_t *w, real_t lp, real_t ls, real_t li, real_t vp, real_t vs, real_t vi, real_t b2p, real_t b2s, real_t b2i,  real_t b3p, real_t b3s, real_t b3i, real_t z)
{
	
	real_t attenp = expf(-0.5*alpha_crp*z);
	real_t attens = expf(-0.5*alpha_crs*z);
	real_t atteni = expf(-0.5*alpha_cri*z);
	
	for (uint idx = 0; idx < SIZE; idx++){
		
		auxp[idx] = Apw[idx] * CpxExp(z*w[idx]*((1/vs-1/vp)+0.5*w[idx]*b2p + w[idx]*w[idx]*b3p/6));
		auxs[idx] = Asw[idx] * CpxExp(z*w[idx]*((1/vs-1/vs)+0.5*w[idx]*b2s + w[idx]*w[idx]*b3s/6));
		auxi[idx] = Aiw[idx] * CpxExp(z*w[idx]*((1/vs-1/vi)+0.5*w[idx]*b2i + w[idx]*w[idx]*b3i/6));
		
		Apw[idx] = auxp[idx] * attenp;
		Asw[idx] = auxs[idx] * attens;
		Aiw[idx] = auxi[idx] * atteni;
		
	}
	
	return ;
}


/** Compute the evolution of the electric fields for a single-pass using the SSFM.
 * The nonlinear crystal is divided into NZ slices with a size of dz. 
 * The SSFM is performed with the following steps:
 * 	1 - The nonlinear part is solved using RK4 for the first semistep dz/2
 * 	2 - The linear part is solved in the frequency domain for the full step dz.
 * 	3 - Repeat 1 for dz/2. 
 * 	4-  Repeat steps 1-3 until finishing the crystal
 */
void EvolutionInCrystal( real_t *w_ext, complex_t *Ap, complex_t *As, complex_t *Ai, complex_t *Apw, complex_t *Asw, complex_t *Aiw, complex_t *k1p, complex_t *k1s, complex_t *k1i, complex_t *k2p, complex_t *k2s, complex_t *k2i, complex_t *k3p, complex_t *k3s, complex_t *k3i, complex_t *k4p, complex_t *k4s, complex_t *k4i, complex_t *auxp, complex_t *auxs, complex_t *auxi, real_t lp, real_t ls, real_t li, real_t vp, real_t vs, real_t vi, real_t b2p, real_t b2s, real_t b2i, real_t b3p, real_t b3s, real_t b3i, real_t dk, real_t kp, real_t ks, real_t ki, real_t dz )
{
	
	// Set plan for cuFFT //
	fftwf_plan plan_1 = NULL; // c2c for input field
	
	real_t z = 0;
	for (uint s = 0; s < NZ; s++){
		/* First RK4 for dz/2 */
		//k1 = dAdz(kappas,dk,z,A)
		dAdz( k1p, k1s, k1i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z );
		//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k1p, k1s, k1i, 0.5 );
		dAdz( k2p, k2s, k2i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/4 );
		// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k2p, k2s, k2i, 0.5 );
		dAdz( k3p, k3s, k3i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/4 );
		// k4 = dAdz(kappas,dk,z+dz,A+k3)
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k3p, k3s, k3i, 1.0 );
		dAdz( k4p, k4s, k4i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/2 );
		// A = A+(k1+2*k2+2*k3+k4)*dz/6
		rk4( Ap, As, Ai, k1p, k1s, k1i, k2p, k2s, k2i, k3p, k3s, k3i, k4p, k4s, k4i, dz/2 );
		
		// Linear operator for dz
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Ap), reinterpret_cast<fftwf_complex*>(Apw), FFTW_BACKWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		FFTScale (Apw, SIZE);
		
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(As), reinterpret_cast<fftwf_complex*>(Asw), FFTW_BACKWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		FFTScale (Asw, SIZE);
		
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Ai), reinterpret_cast<fftwf_complex*>(Aiw), FFTW_BACKWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		FFTScale (Aiw, SIZE);
		
		LinearOperator( auxp, auxs, auxi, Apw, Asw, Aiw, w_ext, lp, ls, li, vp, vs, vi, b2p, b2s, b2i, b3p, b3s, b3i, dz );
		
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Apw), reinterpret_cast<fftwf_complex*>(Ap), FFTW_FORWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Asw), reinterpret_cast<fftwf_complex*>(As), FFTW_FORWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		
		plan_1 = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Aiw), reinterpret_cast<fftwf_complex*>(Ai), FFTW_FORWARD, FFTW_MEASURE);
		fftwf_execute(plan_1);
		
		
		/* Second RK4 for dz/2 */
		//k1 = dAdz(kappas,dk,z,A)
		dAdz( k1p, k1s, k1i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z );
		//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k1p, k1s, k1i, 0.5 );
		dAdz( k2p, k2s, k2i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/4 );
		// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k2p, k2s, k2i, 0.5 );
		dAdz( k3p, k3s, k3i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/4 );
		// k4 = dAdz(kappas,dk,z+dz,A+k3)
		LinealCombination( auxp, auxs, auxi, Ap, As, Ai, k3p, k3s, k3i, 1.0 );
		dAdz( k4p, k4s, k4i, Ap, As, Ai, lp, ls, li, kp, ks, ki, dk, z+dz/2 );
		// A = A+(k1+2*k2+2*k3+k4)*dz/6
		rk4( Ap, As, Ai, k1p, k1s, k1i, k2p, k2s, k2i, k3p, k3s, k3i, k4p, k4s, k4i, dz/2 );
		
		z+=dz;
	}
	
	fftwf_destroy_plan(plan_1);
	
	return ;
}



#endif // -> #ifdef _CWES3CUH
