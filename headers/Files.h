/*---------------------------------------------------------------------------*/
// * This file contains four functions that save files in .dat extension
// * 1 - SaveFileVectorReal()       : save CPU real vectors 
// * 2 - SaveFileVectorComplex()    : save CPU complex vectors

// Inputs:
// - Vector   : vector to save (stored on CPU or GPU)
// - Filename : name of the saved file
/*---------------------------------------------------------------------------*/


#ifndef _FILES
#define _FILES


void SaveVectorReal (rVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension = ".dat";
	myfile.open(Filename+extension);
	for (int i = 0; i < V.size(); i++)
		myfile << std::setprecision(20) << V[i] << "\n";
	myfile.close();
	
	return;
}


void SaveVectorComplex (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].real() << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].imag() << "\n";
	myfile.close();
	
	return;
	
}

#endif // -> #ifdef _FILES