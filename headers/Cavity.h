/*---------------------------------------------------------------------------*/
// * This file contains functions to -----
/*---------------------------------------------------------------------------*/

#ifndef _CAVITYH
#define _CAVITYH


/** Add phase (phase) and mirror losses (R) after a single-pass */
void roundTripChange( cVech_t &A, cVech_t aux, real_t R, real_t phase)
{
	
	for ( uint idx = 0; idx < SIZE; idx++){
		// aux[idx] = (sqrtf(R) * std::exp(Im*phase)) * A[idx];
		// A[idx]   = aux[idx];
		A[idx] = sqrtf(R) * A[idx];
	}
	
	return ;
}


template<typename Crystal>
class Cavity
{	// Difine the class Cavity

public:	
	real_t R, Lcav, trt, fsr; // cavity 
	real_t gamma;             // GDD comp
	bool gdd;
	real_t beta, fpm;		  // EOM	
	bool eom;
	Crystal *Cr;

	Cavity( Crystal *_Cr, real_t _Lcav, real_t _R ) : Cr(_Cr), Lcav(_Lcav), R(_R)
	{	// Constructor
		trt = ( Lcav + (Cr->Lcr)*((Cr->ns) - 1) )/C;
		fsr = 1.0f/trt;

		printf("\nInstance of the class Cavity.\n");
	}

	~Cavity()
	{	// Destructor
		printf("Cavity Destructor.\n");
	}

	// Methods definition
	void applyReflectivity(cVech_t &Ax, real_t phase);
	void setEOM(real_t _beta, real_t _fpm);
	void setGDD (real_t _gamma);

};


// Methods declaration

template<typename Crystal>
void Cavity<Crystal>::applyReflectivity(cVech_t &Ax, real_t phase)
{
	cVech_t aux(Ax.size());
	roundTripChange(Ax, aux, this->R, phase);
	
	return ;
}

template<typename Crystal>
void Cavity<Crystal>::setEOM(real_t _beta, real_t _fpm)
{	
	this->eom = true; 
	this->fpm = _fpm; this->beta = _beta; 
	std::cout << "Set intracavity EOM with \u03B2 = " << this->beta << "- fpm = " << this->fpm << std::endl;
	return ;
}

template<typename Crystal>
void Cavity<Crystal>::setGDD (real_t _gamma)
{
	// set gamma between 0 and 1
	this->gdd = true;
	this->gamma = _gamma;
	std::cout << "Set " << static_cast<int>(100*this->gamma)<< "% intracavity GDD compensation.\n" << std::endl;
	return ;
}


#endif // -> #ifdef _CAVITYH