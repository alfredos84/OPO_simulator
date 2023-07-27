/* Author Alfredo Daniel Sanchez: alfredo.daniel.sanchez@gmail.com */

#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <chrono>
#include <sys/time.h>
#include <cufft.h>

#include <fftw3.h>

// Complex data type
using complex_t = cufftComplex;
using real_t = float;

// Define constants
const real_t PI   = 3.14159265358979323846;    // pi
const real_t C    = 299792458*1E6/1E12;        // speed of ligth in vacuum [um/ps]
const real_t EPS0 = 8.8541878128E-12*1E12/1E6; // vacuum pertivity [W.ps/V²μm] 

const uint SIZE   = 1 << 12; // vector size
const uint NZ     = 150;     // steps over z direction
const uint NRT    = 500;   // number of round trips 

// Package headers
#include "headers/common.h"
#include "headers/operators.h"
#ifdef PPLN // Mgo:PPLN nonlinear crystal
#include "headers/ppln.h"
#endif
#ifdef SPPLT // Mgo:sPPLT nonlinear crystal
#include "headers/spplt.h"
#endif
#include "headers/functions.h"
#ifdef THREE_EQS // Define 2 or 3 coupled-wave equations
#include "headers/cwes3.h"
#else
#include "headers/cwes2.h"
#endif
#include "headers/files.h"



int main(int argc, char *argv[]){
	
	////////////////////////////////////////////////////////////////////////////////////////
	// Set GPU and timing
	
	std::cout << "\n\n\n#######---Welcome to OPO calculator---#######\n\n\n" << std::endl;
	
	time_t current_time; // timing the code
	time(&current_time);
	std::cout << "Date: " << ctime(&current_time) << std::endl;
	double iStart = seconds();
	
	
	#ifdef CW_OPO
	std::cout << "Regime: continuous wave" << std::endl;
	#endif
	#ifdef NS_OPO
	std::cout << "Regime: nanosecond" << std::endl;
	#endif
	#ifdef THREE_EQS
	std::cout << "Three equations" << std::endl;
	#else
	std::cout << "Two equations" << std::endl;
	#endif
	
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	//* Define simulation parameters, physical quantities and set electric fields */
	
	// Grids, crystal and cavity parameters //
	real_t lp        = atof(argv[1])*1e-3;  // pump wavelength   [μm]
	real_t ls        = 1.060; //2*lp;           // signal wavelength [μm]
	real_t li        = lp*ls/(ls-lp);  // idler wavelength  [μm]
	
	// 	real_t deff      = 14.77e-6;       // effective second-order susceptibility [um/V]
	real_t Temp      = atof(argv[2]); // crystal temperature [ºC]
	real_t Lambda    = atof(argv[3]); // grating period for QPM [μm]  
	real_t Lcr       = 5e3;            // crystal length [um]
	
	real_t np        = n(lp, Temp);         // pump ref. index
	real_t vp        = group_vel(lp, Temp); // pump group velocity [μm/ps]
	real_t b2p       = gvd(lp, Temp);       // pump GVD [ps²/μm] 
	real_t b3p       = 0.*TOD(lp, Temp);    // pump TOD [ps³/μm]	
	real_t kp        = 2*PI*deff/(np*lp);        // pump kappa [1/V]
	
	real_t ns        = n(ls, Temp);         // signal ref. index
	real_t vs        = group_vel(ls, Temp); // signal group velocity [μm/ps]
	real_t b2s       = gvd(ls, Temp);       // signal GVD [ps²/μm] 
	real_t b3s       = 0.*TOD(ls, Temp);    // signal TOD [ps³/μm]		
	real_t ks        = 2*PI*deff/(ns*ls);        // signal kappa [1/V]
	
	real_t ni        = n(li, Temp);         // idler ref. index
	real_t vi        = group_vel(li, Temp); // idler group velocity [μm/ps]
	real_t b2i       = gvd(li, Temp);       // idler GVD [ps²/μm]
	real_t b3i       = 0.*TOD(li, Temp);    // idler TOD [ps³/μm]	
	real_t ki        = 2*PI*deff/(ni*li);        // idler kappa [1/V]
	
	real_t dk        = 2*PI*( np/lp-ns/ls-ni/li-1/Lambda ); // mismatch factor
	real_t dkp       = 1/vp-1/vs;                           // group velocity mismatch	
	real_t Lcav      = 5 * Lcr;                             // cavity length [um]
	real_t Rs        = atof(argv[5])*0.01;                  // Reflectivity at signal wavelength 
	real_t alphas    = 0.5*((1-Rs)+alpha_crs*Lcr);          // Total losses for threshold condition signal
	#ifdef THREE_EQS
	real_t Ri        = 0.98;                                 // Reflectivity at idler wavelength 
	real_t alphai    = 0.5*((1-Ri)+alpha_cri*Lcr);          // Total losses for threshold condition idler
	#endif
	real_t t_rt      = (Lcav+Lcr*(ns-1))/C;                 // round-trip time [ps]
	real_t FSR       = 1/t_rt;	                          // free-spectral range
	real_t finesse   = 2*PI/(1-Rs);                         // cavity finesse
	real_t lw        = FSR/finesse*1e6;                     // cavity Linewidth [MHz]
	real_t delta     = atof(argv[6]);                       // cavity detuning [rad] 
	real_t epsilon   = atof(argv[7])*0.01;                  // dispersion compensation index
	real_t GDD       = -epsilon*b2s*Lcr;                    // GDD [ps²]
	real_t TODscomp  = -0.01*atof(argv[8])*b3s*Lcr;        // TOD compensation [ps³]
	real_t TODicomp  = -0.01*atof(argv[8])*b3i*Lcr;        // TOD compensation [ps³]
	
	
	// z discretization, time and frequency discretization
	real_t dz        = Lcr/NZ;   // number of z-steps in the crystal
	real_t dT        = t_rt/SIZE; // time step in [ps]
	real_t dF        = 1/t_rt;    // frequency step in [THz]
	
	
	bool video       = false;
	uint Nrts;        // number of last round trips to save (only for cw)
	if(video){Nrts = 100;}
	else{Nrts = 16;}
	
	
	#ifdef CW_OPO
	uint SIZEL = SIZE*Nrts; // size of large vectors for full simulation
	real_t T_width = t_rt*Nrts; // total time for the saved simulation
	#endif
	#ifdef NS_OPO
	uint SIZEL = SIZE*NRT; // size of large vectors for full simulation
	real_t T_width = t_rt*NRT; // total time for the saved simulation
	#endif	
	
	// Time vector T for one round trip
	real_t *T = (real_t*) malloc(sizeof(real_t) * SIZE);
	linspace( T, SIZE, -0.5*t_rt, 0.5*t_rt);
	
	// Time vector Tp for full simulation
	real_t *Tp = (real_t*) malloc(sizeof(real_t) * SIZEL);
	linspace( Tp, SIZEL, -0.5*T_width, 0.5*T_width);
	
	// Time vector Fp for full simulation
	real_t dFp  = 1/T_width;
	real_t *Fp = (real_t*) malloc(sizeof(real_t) * SIZEL);
	linspace( Fp, SIZEL, -0.5*SIZEL*dFp, +0.5*SIZEL*dFp);
	
	// Frequency and angular frequency vectors f and Ω
	real_t *F = (real_t*) malloc(sizeof(real_t) * SIZE);
	linspace( F, SIZE, -0.5*SIZE*dF, +0.5*SIZE*dF);
	real_t *w = (real_t*) malloc(sizeof(real_t) * SIZE);
	fftshift(w,F, SIZE);
	for (uint i=0; i<SIZE; i++){
		w[i] = 2*PI*w[i]; // angular frequency [2*pi*THz]
	}
		
	// Define memory size for complex host vectors
	uint nBytes   = sizeof(complex_t)*SIZE;
	
	// Difine which fields are resonant (SRO, DRO or TRO)
	bool is_Ap_resonant = false;
	bool is_As_resonant = true;
	#ifdef THREE_EQS
	bool is_Ai_resonant = true;
	#endif
	
	// Define input pump parameters
	real_t waist = 55;             // beam waist radius [um]
	real_t spot  = PI*waist*waist; // spot area [μm²]
	real_t Ith, Pth;               // Power and intensity threshold 
	#ifdef THREE_EQS
	// Power and intensity threshold non-degenerate OPO 
	if (!is_Ai_resonant){
		std::cout << "SRO: As is resonant" << std::endl;
		Ith   = EPS0*C*np*ns*ni*ls*li*pow((1/deff/Lcr/PI),2)*alphas/2;
	}
	if (!is_As_resonant){
		std::cout << "SRO, Ai is resonant" << std::endl;
		Ith   = EPS0*C*np*ns*ni*ls*li*pow((1/deff/Lcr/PI),2)*alphai/2;
	}
	if (is_As_resonant and is_Ai_resonant){
		std::cout << "DRO, As and Ai are resonant" << std::endl;
		Ith   = EPS0*C*np*ns*ni*ls*li*pow((1/deff/Lcr/PI),2)*alphas*alphai/8;
	}	
	Pth   = Ith*spot;
	#else
	// Power and intensity threshold degenerate DRO 
	Ith   = EPS0*C*np*powf((ns*ls*alphas/deff/Lcr/PI), 2)/8;
	Pth   = Ith*spot;
	#endif
	
	real_t Nth   = atof(argv[4]);             // Times over the threshold
	real_t Inten = atof(argv[4])*Ith;         // Pump intensity in [W/um²]
	real_t Power = Inten*spot;                // Pump power in [W]
	real_t Ap0   = sqrt(2*Inten/(np*EPS0*C)); // Input pump field strength [V/μm]
	
	// Define input pump vector
	#ifdef CW_OPO
	complex_t *Ap    = (complex_t*)malloc(nBytes); // input pump vector
	complex_t *Ap_in = (complex_t*)malloc(nBytes); // input pump vector
	input_field_T(Ap_in, Ap0, SIZE );              // set input pump vector (cw)
	#endif
	
	#ifdef NS_OPO
	real_t FWHM      = 10000;                              // intensity FWHM for input [ps]
	real_t sigmap    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2))); // σ of electric field gaussian pulse [ps]
	complex_t *Ap    = (complex_t*)malloc(nBytes); // input pump vector
	complex_t *Ap_in = (complex_t*)malloc(sizeof(complex_t)*SIZEL);// input pump vector
	input_field_T(Ap_in, Ap0, Tp, sigmap, SIZEL); // set input pump vector (gaussian pulse)
	#endif
	
	
	// Define input signal vector (NOISE)
	complex_t *As = (complex_t*)malloc(nBytes);
	NoiseGeneratorCPU ( As, SIZE );
	
	#ifdef THREE_EQS	
	// Define input idler vector (NOISE)
	complex_t *Ai = (complex_t*)malloc(nBytes);
	NoiseGeneratorCPU ( Ai, SIZE );
	#endif
	
	
	// Intracavy phase modulator
	bool using_phase_modulator = atoi(argv[9]);
	real_t mod_depth, fpm, df;
	if(using_phase_modulator){
		mod_depth       = atof(argv[10])*PI;
		df              = atof(argv[11])*sqrtf(Nth-1)*alphas/(PI*mod_depth)*FSR;
		fpm             = FSR - df;
	}
	
	
	// Define string variables for saving files
	std::string Filename, SAux, Extension = ".dat";
	bool save_input_fields = false;  // Save input fields files
	if (save_input_fields){
		#ifdef CW_OPO
		Filename = "pump_input";	SaveVectorComplex (Ap_in, SIZE, Filename);
		#endif
		#ifdef NS_OPO
		Filename = "pump_input";	SaveVectorComplex (Ap_in, SIZEL, Filename);
		#endif
		Filename = "signal_input";	SaveVectorComplex (As, SIZE, Filename);
		#ifdef THREE_EQS	
		Filename = "idler_input";	SaveVectorComplex (Ai, SIZE, Filename);	
		#endif
	}
	
	
	bool print_param_on_screen = true;	// Print parameters on screen
	if ( print_param_on_screen ){
		std::cout << "\n\nSimulation parameters:\n\n " << std::endl;
		std::cout << "Number of round trips   = " << NRT  << std::endl;
		std::cout << "Pump wavelength         = " << lp*1e3 << " nm" << std::endl;
		std::cout << "Signal wavelength       = " << ls*1e3 << " nm" << std::endl;
		#ifdef THREE_EQS
		std::cout << "Idler wavelength        = " << li*1e3 << " nm" << std::endl;
		#endif
		std::cout << "Temp                    = " << Temp << " ºC" << std::endl;
		std::cout << "np                      = " << np << std::endl;
		std::cout << "ns                      = " << ns << std::endl;
		std::cout << "ni                      = " << ni << std::endl;
		std::cout << "\u03BD⁻¹ pump                = " << 1.0/vp << " ps/\u03BCm" << std::endl;
		std::cout << "\u03BD⁻¹ signal              = " << 1.0/vs << " ps/\u03BCm" << std::endl;
		#ifdef THREE_EQS
		std::cout << "\u03BD⁻¹ idler               = " << 1.0/vi << " ps/\u03BCm" << std::endl;		
		#endif
		std::cout << "\u0394k                      = " << dk << " \u03BCm⁻¹" << std::endl;
		std::cout << "\u0394k'                     = " << dkp << " ps/\u03BCm" << std::endl;	
		std::cout << "GVD pump                = " << b2p << " ps²/\u03BCm" << std::endl;
		std::cout << "GVD signal              = " << b2s << " ps²/\u03BCm" << std::endl;
		#ifdef THREE_EQS
		std::cout << "GVD idler               = " << b2i << " ps²/\u03BCm" << std::endl;		
		#endif
		std::cout << "TOD pump                = " << b3p << " ps³/\u03BCm" << std::endl;
		std::cout << "TOD signal              = " << b3s << " ps³/\u03BCm" << std::endl;		
		std::cout << "Net GVD                 = " << (1-epsilon)*b2s << " ps²/\u03BCm" << std::endl;
		std::cout << "GVD compensation        = " << atoi(argv[7]) << " %"  << std::endl;
		std::cout << "Net TOD                 = " << (1-0.01*atoi(argv[8]))*b3s*Lcr*1e3 << " fs³"  << std::endl;
		std::cout << "TOD compensation        = " << atof(argv[8]) << " %"  << std::endl;		
		std::cout << "deff                    = " << deff*1e6 << " pm/V"  << std::endl;
		std::cout << "\u039B                       = " << Lambda << " \u03BCm"  << std::endl;
		std::cout << "\u03B1cp                     = " << alpha_crp << " \u03BCm⁻¹"  << std::endl;
		std::cout << "\u03B1cs                     = " << alpha_crs << " \u03BCm⁻¹" << std::endl;
		std::cout << "\u03B1s                      = " << alphas << std::endl;
		#ifdef THREE_EQS
		std::cout << "\u03B1ci                      = " << alpha_cri << " \u03BCm⁻¹" << std::endl;
		std::cout << "\u03B1i                      = " << alphas << std::endl;
		#endif
		
		std::cout << "Crystal length          = " << Lcr*1e-3 << " mm"  << std::endl;
		std::cout << "Cavity  length          = " << Lcav*1e-3 << " mm"  << std::endl;
		std::cout << "\u0394z                      = " << dz << " \u03BCm"  << std::endl;
		std::cout << "Reflectivity (signal)   = " << Rs*100 << " %"  << std::endl;
		#ifdef THREE_EQS
		std::cout << "Reflectivity (idler)    = " << Ri*100 << " %"  << std::endl;	
		#endif
		std::cout << "Cavity Finesse          = " << finesse << std::endl;	
		std::cout << "Cavity lw (FWHM)        = " << lw << " MHz"  << std::endl;	
		std::cout << "Round-trip time         = " << std::setprecision(15) << t_rt << " ps"  << std::endl;	
		std::cout << "FSR                     = " << std::setprecision(15) << FSR*1e3 << " GHz"  << std::endl;
		std::cout << "Cavity detuning (\u03B4)     = " << delta << "\u03C0"  << std::endl;	
		std::cout << "Using N                 = " << SIZE << " points" << std::endl;
		std::cout << "dT                      = " << dT << " ps" << std::endl;
		std::cout << "SIZEL                   = " << SIZEL << std::endl;
		std::cout << "dFp                     = " << dFp << " THz" << std::endl;
		std::cout << "Max frequency           = " << Fp[SIZEL-1] << " THz" << std::endl;
		std::cout << "Ap0                     = " << Ap0 << " V/um" << std::endl; 
		std::cout << "waist                   = " << waist << " \u03BCm" << std::endl;
		std::cout << "spot                    = " << spot << " \u03BCm²" << std::endl;
		std::cout << "Power threshold         = " << Pth << " W" << std::endl;
		std::cout << "Power                   = " << Power << " W" << std::endl;
		std::cout << "Times above the thres.  = " << Nth << std::endl;
		if(using_phase_modulator){
			std::cout << "Using a phase modulator:" << std::endl;
			std::cout << "Mod. depth (\u03B2)          = " << atof(argv[10]) << "\u03C0 rad = " << mod_depth << " rad" << std::endl;
			std::cout << "Freq. detuning (\u03B4f)     = " << df*1e6 << " MHz" << std::endl;
			std::cout << "Mod. frequency(fm)      = " << fpm*1e3 << " GHz" << std::endl;
			std::cout << "\n\nPoint in the space of parameters:\n" << std::endl;
			std::cout << "(N,\u03B2,\u03B4f,\u03B5) = ( " << Nth << ", " << atof(argv[10]) << ", "  << std::setprecision(4) << df*1e6 << ", " << epsilon << " )\n\n" << std::endl;			
		}
		else{std::cout << "No phase modulator" << std::endl;
			std::cout << "\n\nPoint in the space of parameters:\n" << std::endl;
			std::cout << "( N, \u03B2, \u03B4f, \u03B5 ) = ( " << Nth << ", 0, 0, " << std::setprecision(2) << epsilon << " )\n\n" << std::endl;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	//* Define CPU vectors //

	complex_t *As_total = (complex_t*)malloc(sizeof(complex_t)*SIZEL);
	complex_t *Ap_total = (complex_t*)malloc(sizeof(complex_t)*SIZEL);
	complex_t *Asw = (complex_t*)malloc(nBytes);
	complex_t *Apw = (complex_t*)malloc(nBytes);

	
	// RK4 (kx) and auxiliary (aux) CPU vectors 
	complex_t *k1p = (complex_t*)malloc(nBytes);
	complex_t *k2p = (complex_t*)malloc(nBytes);
	complex_t *k3p = (complex_t*)malloc(nBytes);
	complex_t *k4p = (complex_t*)malloc(nBytes);
	complex_t *k1s = (complex_t*)malloc(nBytes);
	complex_t *k2s = (complex_t*)malloc(nBytes);
	complex_t *k3s = (complex_t*)malloc(nBytes);
	complex_t *k4s = (complex_t*)malloc(nBytes);
	complex_t *auxp = (complex_t*)malloc(nBytes);
	complex_t *auxs = (complex_t*)malloc(nBytes);
	

	#ifdef THREE_EQS
	complex_t *Aiw = (complex_t*)malloc(nBytes);
	complex_t *Ai_total = (complex_t*)malloc(sizeof(complex_t)*SIZEL);
	complex_t *k1i = (complex_t*)malloc(nBytes);
	complex_t *k2i = (complex_t*)malloc(nBytes);
	complex_t *k3i = (complex_t*)malloc(nBytes);
	complex_t *k4i = (complex_t*)malloc(nBytes);
	complex_t *auxi = (complex_t*)malloc(nBytes);
	bool idler_pm = true; // phase modulator applies on idler
	#endif
	
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////	
	//* Main loop. Fields in the cavity */
	
	// Set plan for FFT //
	fftwf_plan plan = NULL; // c2c for input field

	std::cout << "Starting main loop on CPU & GPU...\n" << std::endl;
	uint mm = 0; // counts for cw saved round trips
	for (uint nn = 0; nn < NRT; nn++){
		if( nn%10 == 0 or nn == NRT-1 )
			std::cout << "#round trip: " << nn << std::endl;
		
		#ifdef CW_OPO
		// update the input pump in each round trip
		memcpy( Ap, Ap_in, nBytes );
		#endif
		#ifdef NS_OPO
		// read the input pump in nn-th round trip
		ReadPump ( Ap, Ap_in, NRT, nn, SIZE );
		#endif
		
		#ifdef THREE_EQS
		if (!is_Ai_resonant){	// For non-resonant field, it resets Ai in every round trip
			NoiseGeneratorCPU ( Ai, SIZE );
		}
		#endif
		
		if (!is_As_resonant){	// For non-resonant field, it resets As in every round trip
			NoiseGeneratorCPU ( As, SIZE );
		}
		
		#ifdef THREE_EQS // Single pass for coupled wave equations (2 or 3)
		EvolutionInCrystal( w, Ap, As, Ai, Apw, Asw, Aiw, k1p, k1s, k1i, k2p, k2s, k2i, k3p, k3s, k3i, k4p, k4s, k4i, auxp, auxs, auxi, lp, ls, li, vp, vs, vi, b2p, b2s, b2i, b3p, b3s, b3i, dk, kp, ks, ki, dz );
		#else
		EvolutionInCrystal( w, Ap, As, Apw, Asw, k1p, k1s, k2p, k2s, k3p, k3s, k4p, k4s, auxp, auxs, lp, ls, vp, vs, b2p, b2s, b3p, b3s, dk, kp, ks, dz );
		#endif
		
		
		if(GDD!=0){ // adds dispersion compensation
			AddGDD (Asw, auxs, w, GDD);
			plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Asw), reinterpret_cast<fftwf_complex*>(As), FFTW_FORWARD, FFTW_MEASURE);
			fftwf_execute(plan);
			#ifdef THREE_EQS
			AddGDD (Aiw, auxi, w, GDD);
			plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Aiw), reinterpret_cast<fftwf_complex*>(Ai), FFTW_FORWARD, FFTW_MEASURE);
			fftwf_execute(plan);

			#endif
		}		
		
		if( using_phase_modulator ){ // use an intracavy phase modulator of one o more fields
			PhaseModulatorIntraCavity(As, auxs, mod_depth, fpm, T);
			#ifdef THREE_EQS
			if(idler_pm){
				PhaseModulatorIntraCavity(Ai, auxi, mod_depth, fpm, T);
			}
			#endif
		}
		
		if (is_As_resonant){ // if As is resonant, adds phase and losses
			AddPhase(As, auxs, Rs, delta, nn);
		}
		
		#ifdef THREE_EQS
		if (is_Ai_resonant){  // if Ai is resonant, adds phase and losses
			AddPhase(Ai, auxi, Ri, delta, nn);
		}
		#endif
		
		#ifdef CW_OPO	// saves systematically every round trip
		if (video){  
			if (nn % 100 == 0){ // this branch is useful if the user want to save the round trips every 100 ones
				std::cout << "Saving the " << nn << "-th round trip" << std::endl;
				SaveRoundTrip(As_total, As, mm, Nrts ); // saves signal
				SaveRoundTrip(Ap_total, Ap, mm, Nrts ); // saves pump
				#ifdef THREE_EQS
				SaveRoundTrip(Ai_total, Ai, mm, Nrts ); // saves idler
				#endif
				mm += 1;
			}			
		}
		else{  // this branch is useful if the user want to save the last NRT-Nrts round trips
			if (nn >= NRT -Nrts){                
				SaveRoundTrip( Ap_total, Ap, mm, Nrts ); // saves pump
				SaveRoundTrip( As_total, As, mm, Nrts ); // saves signal
				#ifdef THREE_EQS
				SaveRoundTrip( As_total, As, mm, Nrts ); // saves idler
				#endif
				mm += 1;
			}
		}
		#endif
		#ifdef NS_OPO	// save the simulation in the NS regime
		SaveRoundTrip(Ap_total, Ap, nn, NRT ); // saves pump
		SaveRoundTrip(As_total, As, mm, NRT ); // saves signal
		#ifdef THREE_EQS
		SaveRoundTrip(Ai_total, Ai, mm, NRT ); // saves idler
		#endif
		#endif
		
	} // End of main loop
	
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////	
	// Saving results in .dat files using the function SaveFileVectorComplex() //
	 
	bool save_vectors = true; // Decide whether or not save these vectors
	if (save_vectors){
		std::cout << "\nSaving time and frequency vectors...\n" << std::endl;
		Filename = "Tp"; SaveVectorReal (Tp, SIZEL, Filename+Extension);
		Filename = "freq"; SaveVectorReal (Fp, SIZEL, Filename+Extension);
		Filename = "T"; SaveVectorReal (T, SIZE, Filename+Extension);
	}
	else{ std::cout << "\nTime and frequency were previuosly save...\n" << std::endl;
	}
	
	// Save the full simulation
	Filename = "signal_output";	SaveVectorComplex ( As_total, SIZEL, Filename );
	Filename = "pump_output";	SaveVectorComplex ( Ap_total, SIZEL, Filename );
	#ifdef THREE_EQS
	Filename = "idler_output";	SaveVectorComplex ( Ai_total, SIZEL, Filename );
	#endif
	
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	// Deallocating memory from CPU and destroying plans
	
	free(Tp); free(T); free(Fp); free(w); free(F);
	free(Ap); free(As); free(Ap_in);
	free(k1p); free(k2p); free(k3p); free(k4p);
	free(k1s); free(k2s); free(k3s); free(k4s);
	free(auxp); free(auxs);
	
	#ifdef THREE_EQS
	free(Ai); free(Ai_total);
	free(k1i); free(k2i); free(k3i); free(k4i);
	free(auxi);
	#endif
	
	// Destroy FFT context
	fftwf_destroy_plan(plan); 
	
	////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////	
	// Finish timing: returns the runtime simulation
	
	double iElaps = seconds() - iStart;
	if(iElaps>60){std::cout << "\n\n...time elapsed " <<  iElaps/60.0 << " min\n\n " << std::endl;}
	else{std::cout << "\n\n...time elapsed " <<  iElaps << " seconds\n\n " << std::endl;}
	
	time(&current_time);
	std::cout << ctime(&current_time) << std::endl;
	
	return 0;
}
