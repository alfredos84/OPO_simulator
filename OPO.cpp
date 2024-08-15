// * Author name: Alfredo Daniel Sanchez
// * email:       alfredo.daniel.sanchez@gmail.com

// * Compile with: 
// * g++ OPO.cpp -o OPO -O3 -lfftw3 -lfftw3f -lm

#include "headers/Libraries.h"		// Required libraries
#include "headers/PackageLibraries.h"	// Required package libraries

int main(int argc, char *argv[]){
	
	///////////////////////////////////////////////////////////////////////////////////
	// 1. Build the OPO

	// a. set parameters
	real_t lp=0.532, ls=2*lp, li=ls*lp/(ls-lp), waist=55.0f, reflectivity = 0.98f;
	real_t Lcr = 20.e3f, T = 27.0f, Lambda = 6.97f;


	// b. set crystal
	MgOsPPLT * cr1 = new MgOsPPLT(Lcr, T, Lambda, lp, ls, li);
	cr1->getCrystalProp(); // print on screen crystal properties


	// c. set cavity
	Cavity<MgOsPPLT> * cav1 = new Cavity<MgOsPPLT>( cr1, 4*cr1->Lcr, reflectivity );
	// cav1->setEOM( 0.8*PI, cav1->fsr ); // set intracavity EOM
	cav1->setGDD( 1.0 ); // set GDD compensation: value in [0,1] for 0 to 100%

	// set pump power from threshold definition
	real_t alphas = 0.5*( (1-cav1->R)+cr1->alpha_crs*cr1->Lcr ), alphai = alphas; 
	real_t Ith   = EPS0*C*cr1->np*cr1->ns*cr1->ni*cr1->ls*cr1->li*powf((1/cr1->dQ/cr1->Lcr/PI),2)*alphas*alphai/8;
	real_t Power = 4*Ith*(waist*waist*PI); 	printVarOnScreen("Power = ", Power);
	
	// d. set electric fields
	EFields<MgOsPPLT> * efs = new EFields<MgOsPPLT>(lp, ls, li, Power, waist, cr1);
	efs->resFields = std::make_tuple(false, true, true); //are (pump,signal,idler) resonant?
	efs->cwField( Power ); efs->noiseGenerator(efs->As); efs->Ai = efs->As;
	efs->setTimeFreqVectors( cav1->trt );
	// real_t FWHM = 0.2; efs->setTimeFreqVectors( 50*FWHM );	// set time and freq vectors
	// efs->gaussianField( Power, FWHM );	efs->noiseGenerator(efs->As); efs->noiseGenerator(efs->Ai); 


	///////////////////////////////////////////////////////////////////////////////////
	// 2. Run the OPO
	
	Solver<MgOsPPLT> * solv1 = new Solver<MgOsPPLT>(cr1, cav1, efs);
	std::vector<int> save_roundtrips(64); // set of roundtrips to save
	for (uint i = 0; i < save_roundtrips.size(); i++)
		save_roundtrips[i] = NRT - save_roundtrips.size() + i;
	solv1->runOPO(save_roundtrips);


	// save relevant vectors
	SaveVectorReal(solv1->A->t, "time");
	SaveVectorReal(solv1->A->F, "freq");
	SaveVectorComplex(solv1->A->Ap, "output_pump");
	SaveVectorComplex(solv1->A->As, "output_signal");
	SaveVectorComplex(solv1->A->Ai, "output_idler");
	
	// delete objects
	delete solv1; delete efs, delete cav1; delete cr1;
  	
	///////////////////////////////////////////////////////////////////////////////////
	std::cout << "\t---> Simulation finished <---" << std::endl;
	
	return 0;
	
}