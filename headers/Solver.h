/*---------------------------------------------------------------------------*/
// * This file contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution through the nonlinear Cr.
// * 
// * In particular, this file should be used when only three equation describes the 
// * problem, i.e., sum or difference frequency generation (SFG or DFG).
// *
// * For any specific process, please check the form of the Couple wave equations
// * in the first function called dAdz(). 
/*---------------------------------------------------------------------------*/



#ifndef _SOLVERH
#define _SOLVERH


/** Size-scaling vectors after Fourier transform */
inline void realScaling ( cVech_t &A )
{
	real_t s = static_cast<real_t>(SIZE);
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		A[idx] = A[idx] / s ;
	}
	
	return ;
}



/** Scales all vectors after Fourier transform */
inline void CUFFTScale ( cVech_t &Awp, cVech_t &Aws, cVech_t &Awi )
{
	real_t s = static_cast<real_t>(SIZE);
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		Awp[idx] = Awp[idx] / s ; Aws[idx] = Aws[idx] / s ;	Awi[idx] = Awi[idx] / s ;
	}
	
	return ;
}


/** This function compensates the GVD after a single-pass */
inline void addGDD(cVech_t &A, cVech_t aux, rVech_t w, real_t GDD)
{
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		aux[idx] = A[idx] * std::exp( Im*0.5f*w[idx]*w[idx]*GDD );
		A[idx] = aux[idx];
	}
	
	return ;
}

/** Applies an electro optical modulator to an electric field after one round trip. */
inline void intracavityEOM( cVech_t &A, cVech_t aux, real_t beta, real_t fpm, rVech_t t)
{
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		aux[idx] = A[idx] * std::exp(Im*beta*sinf(2.0f*PI*fpm*t[idx]));
		A[idx] = aux[idx];}
		
	return ;
}


/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x represents the different fields) */
inline void dAdz( cVech_t &dAp, cVech_t &dAs,  cVech_t &dAi, cVech_t Ap, cVech_t As, cVech_t Ai,
			real_t kp, real_t ks, real_t ki, real_t dk, real_t z )
{
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		dAp[idx]  = Im * kp * As[idx] * Ai[idx] * std::exp(-Im*0.0f*dk*z) ;
		dAs[idx]  = Im * ks * Ap[idx] * std::conj(Ai[idx]) * std::exp(+Im*0.0f*dk*z);
		dAi[idx]  = Im * ki * Ap[idx] * std::conj(As[idx]) * std::exp(+Im*0.0f*dk*z);
	}
	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
inline void LinealCombination( cVech_t &auxp, cVech_t &auxs, cVech_t &auxi, 
						cVech_t Ap, cVech_t As, cVech_t Ai, 
						cVech_t kp, cVech_t ks, cVech_t ki, real_t s ){

	for ( uint idx = 0; idx < SIZE; idx++ ){
		auxp[idx] = Ap[idx] + kp[idx] * s;
		auxs[idx] = As[idx] + ks[idx] * s;
		auxi[idx] = Ai[idx] + ki[idx] * s;
	}
	
	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
inline void rk4(cVech_t &Ap, cVech_t &As,  cVech_t &Ai, cVech_t k1p, cVech_t k1s, cVech_t k1i, 
					cVech_t k2p, cVech_t k2s, cVech_t k2i, cVech_t k3p, cVech_t k3s, cVech_t k3i,
					cVech_t k4p, cVech_t k4s, cVech_t k4i, real_t dz ){
		
	for ( uint idx = 0; idx < SIZE; idx++ ){
		Ap[idx] = Ap[idx] + (k1p[idx] + 2.0f*k2p[idx] + 2.0f*k3p[idx] + k4p[idx]) * dz / 6.0f;
		As[idx] = As[idx] + (k1s[idx] + 2.0f*k2s[idx] + 2.0f*k3s[idx] + k4s[idx]) * dz / 6.0f;
		Ai[idx] = Ai[idx] + (k1i[idx] + 2.0f*k2i[idx] + 2.0f*k3i[idx] + k4i[idx]) * dz / 6.0f;
	}
	
	return ;
}


/** Computes the linear part: Ax = Ax.exp(i.f(Ω)*z), where f(Ω) is a frequency dependant functions
 * including the group velocity and the group velocity dispersion parameters. */
inline void LinearOperator(cVech_t &auxp, cVech_t &auxs, cVech_t &auxi, cVech_t &Awp, cVech_t &Aws, cVech_t &Awi,
								rVech_t w, real_t vp, real_t vs, real_t vi, real_t b2p, real_t b2s, real_t b2i, 
								real_t b3p, real_t b3s, real_t b3i, real_t alpha_crp,
								real_t alpha_crs, real_t alpha_cri, real_t z)
{
	
	real_t attenp = std::exp(-0.5f*alpha_crp*z);
	real_t attens = std::exp(-0.5f*alpha_crs*z);
	real_t atteni = std::exp(-0.5f*alpha_cri*z);
	
	for ( uint idx = 0; idx < SIZE; idx++ ){
		auxp[idx] = Awp[idx] * std::exp(Im*z*w[idx]*((1/vp-1/vi)+0.5f*w[idx]*b2p + w[idx]*w[idx]*b3p/6.0f));
		auxs[idx] = Aws[idx] * std::exp(Im*z*w[idx]*((1/vs-1/vs)+0.5f*w[idx]*b2s + w[idx]*w[idx]*b3s/6.0f));
		auxi[idx] = Awi[idx] * std::exp(Im*z*w[idx]*((1/vi-1/vs)+0.5f*w[idx]*b2i + w[idx]*w[idx]*b3i/6.0f));
	
		Awp[idx] = auxp[idx] * attenp;
		Aws[idx] = auxs[idx] * attens;
		Awi[idx] = auxi[idx] * atteni;
	}
	
	return ;
}


template<typename Crystal>
class Solver
{	// Difine the class Solver for modelling the fields propagation and heating
public:	

	cVech_t k1p, k2p, k3p, k4p;
	cVech_t k1s, k2s, k3s, k4s;
	cVech_t k1i, k2i, k3i, k4i;
	cVech_t auxp, auxs, auxi;

	Crystal *Cr;
	EFields<Crystal> *A;
	Cavity<Crystal> *Cav;

	Solver(Crystal *_Cr, EFields<Crystal> *_A) : Cr(_Cr), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

		printf("\nInstance of the class Solver.\n");
	}


	Solver(Crystal *_Cr, Cavity<Crystal> *_Cav, EFields<Crystal> *_A) : Cr(_Cr), Cav(_Cav), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

		printf("\nInstance of the class Solver.\n");
	}
	
	~Solver(){ printf("Solver Destructor.\n"); }

	// Methods definition
	void EOM ( cVech_t &Ax );
	void GDD( cVech_t &Ax );
	void dispersion ( );
	void solverRK4( real_t z );
	void SSFM ( );
	void runOPO(const std::vector<int>& save_roundtrips);
	void runOPO( );
	void runSinglePass( );

};


// Methods declaration

template<typename Crystal>
inline void Solver<Crystal>::EOM ( cVech_t &Ax )
{
	cVech_t aux(Ax.size());			
	intracavityEOM( Ax, aux, Cav->beta, Cav->fpm, A->t );
	return ;
}


template<typename Crystal>
inline void Solver<Crystal>::GDD( cVech_t &Ax )
{
	
	cVech_t Awx(Ax.size());	cVech_t aux(Ax.size());

	fftwf_plan plan = NULL; // c2c for input field

    plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Ax.data()), reinterpret_cast<fftwf_complex*>(Awx.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	realScaling( Awx );
	
	real_t GDD = -(Cav->gamma)*(Cr->b2s)*(Cr->Lcr);
	addGDD(Awx, aux, A->w, GDD);

    plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(Awx.data()), reinterpret_cast<fftwf_complex*>(Ax.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	
	fftwf_destroy_plan(plan);

	return ;
}


template<typename Crystal>
inline void Solver<Crystal>::dispersion ( )
{	// A->Applies the dispersion term to the electric fields
	// Parameters for kernels
	
	// Set plan for cuFFT //
	fftwf_plan plan = NULL;	

	real_t vp  = Cr->vp,  vs  = Cr->vs,  vi  = Cr->vs;
	real_t b2p = Cr->b2p, b2s = Cr->b2s, b2i = Cr->b2i; 
	real_t b3p = Cr->b3p, b3s = Cr->b3s, b3i = Cr->b3i;
	real_t alpha_crp = Cr->alpha_crp, alpha_crs = Cr->alpha_crs, alpha_cri = Cr->alpha_cri;
	real_t dz = Cr->dz;

	// Linear operator for dz
    plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(A->Ap.data()), reinterpret_cast<fftwf_complex*>(A->Awp.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	
    plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(A->As.data()), reinterpret_cast<fftwf_complex*>(A->Aws.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);

	plan = fftwf_plan_dft_1d(SIZE, reinterpret_cast<fftwf_complex*>(A->Ai.data()), reinterpret_cast<fftwf_complex*>(A->Awi.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	
	CUFFTScale(A->Awp, A->Aws, A->Awi);

	LinearOperator( this->auxp, this->auxs, this->auxi, A->Awp, A->Aws, A->Awi, A->w,
					vp, vs, vi, b2p, b2s, b2i, b3p, b3s, b3i, alpha_crp, alpha_crs, alpha_cri ,dz);

	plan = fftwf_plan_dft_1d(SIZE,  reinterpret_cast<fftwf_complex*>(A->Awp.data()), reinterpret_cast<fftwf_complex*>(A->Ap.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	
	plan = fftwf_plan_dft_1d(SIZE,  reinterpret_cast<fftwf_complex*>(A->Aws.data()), reinterpret_cast<fftwf_complex*>(A->As.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);

	plan = fftwf_plan_dft_1d(SIZE,  reinterpret_cast<fftwf_complex*>(A->Awi.data()), reinterpret_cast<fftwf_complex*>(A->Ai.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
	
	fftwf_destroy_plan(plan);
	
	return ;
}


template<typename Crystal>
inline void Solver<Crystal>::solverRK4( real_t z )
{	
	// A->Applies the Fourh-order Runge-Kutta Method with fixed step size dz
	// This function apply the fourth-order Runge-Kutta method	
		
	real_t dz = (Cr->dz)*0.5f;
	real_t dk = A->dk; 
	real_t kp = A->kp, ks = A->ks, ki = A->ki;

	//k1 = dAdz(kappas,dk,z,A)
	dAdz( this->k1p, this->k1s, this->k1i, A->Ap, A->As, A->Ai, kp, ks, ki, dk, z );
 
	//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
	LinealCombination( this->auxp, this->auxs, this->auxi, A->Ap, A->As, A->Ai, this->k1p, this->k1s, this->k1i, 0.5f );
	dAdz( this->k2p, this->k2s, this->k2i, this->auxp, this->auxs, this->auxi, kp, ks, ki, dk, z+dz/2.0f );

	// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
	LinealCombination( this->auxp, this->auxs, this->auxi, A->Ap, A->As, A->Ai, this->k2p, this->k2s, this->k2i, 0.5f );
	dAdz( this->k3p, this->k3s, this->k3i, this->auxp, this->auxs, this->auxi, kp, ks, ki, dk, z+dz/2.0f );

	// k4 = dAdz(kappas,dk,z+dz,A+k3)
	LinealCombination( this->auxp, this->auxs, this->auxi, A->Ap, A->As, A->Ai, this->k3p, this->k3s, this->k3i, 1.0f );
	dAdz( this->k4p, this->k4s, this->k4i, this->auxp, this->auxs, this->auxi, kp, ks, ki, dk, z+dz );

	// A = A + (k1+k2+k3+k4)/6
	rk4( A->Ap, A->As, A->Ai, this->k1p, this->k1s, this->k1i, this->k2p, this->k2s, this->k2i, 
		this->k3p, this->k3s, this->k3i, this->k4p, this->k4s, k4i, dz );


	return ;
	
}


template<typename Crystal>
inline void Solver<Crystal>::SSFM ( )
{
	for (real_t z = 0; z <= Cr->Lcr; z+=(Cr->dz))
	{	
		solverRK4(z); // RK4 in dz/2
		dispersion(); // Dispersion in dz
		solverRK4(z); // RK4 in dz/2
	}

	return ;
}


template<typename Crystal>
void Solver<Crystal>::runOPO(const std::vector<int>& save_roundtrips)
{
	
	std::unordered_set<int> save_set(save_roundtrips.begin(), save_roundtrips.end());

	real_t iStart = seconds();	// Timing code
	
	for( uint r = 0; r < NRT; r++ )
	{

		runOPO_status (r, NRT/10); // print simulation status on screen
		// std::cout << "# Round trip: " << r << " - Completed " << r*100/NRT << "%" << "\t\r" << std::flush;
		SSFM(); // compute couple-wave equations

		// intracavity elements
		if(Cav->eom){EOM(A->As); EOM(A->Ai);}
		if(Cav->gdd){GDD(A->As); GDD(A->Ai);}

		// Apply boundary conditions to resonant fields	
		if(std::get<1>(A->resFields)){Cav->applyReflectivity(A->As, 0.0f);}
		if(std::get<2>(A->resFields)){Cav->applyReflectivity(A->Ai, 0.0f);}
		A->cwField( A->Power );

		// save specific round trips
		if (save_set.find(r) != save_set.end()) {
			SaveVectorComplex(A->Ap, "output_pump_"+std::to_string(r));
			SaveVectorComplex(A->As, "output_signal_"+std::to_string(r));
			SaveVectorComplex(A->Ai, "output_idler_"+std::to_string(r));
        }
	}

	real_t iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;
}


// Overloaded runOPO() without saving data
template<typename Crystal>
void Solver<Crystal>::runOPO() {runOPO(std::vector<int>());}


template<typename Crystal>
void Solver<Crystal>::runSinglePass()
{
	
	real_t iStart = seconds();	// Timing code

	SSFM();

	real_t iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;
}



#endif // -> #ifdef _SOLVERH