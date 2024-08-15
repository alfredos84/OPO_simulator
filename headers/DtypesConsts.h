#ifndef _DTYPESCONSTSH
#define _DTYPESCONSTSH

/** 
 * Complex data type: a set of datatypes are
 * defined to make the code more readable.
 *
 * Definitions for numbers
 * real_t   : datatype for real numbers
 * complex_t: datatype for complex numbers
 * 
 * Definitions for vectors:
 * 
 * rVech_t  : real vector
 * cVech_t  : complex vector
 */


using real_t    = float;
using complex_t = std::complex<real_t>;
using rVech_t   = std::vector<real_t>;
using cVech_t   = std::vector<complex_t>;


using namespace std::complex_literals;
const complex_t Im = 1if; // imaginary number


// Define global constants
const uint SIZE   = 1 << 10;	// vector size
const uint NZ     = 100;		// size discretization
const uint NRT    = 10000;		// number of round trips    

// Define relevant physical constants
const real_t PI   = std::acos(-1);		// pi
const real_t C    = 299792458*1E6/1E12;			// speed of ligth in vacuum [um/ps]
const real_t EPS0 = 8.8541878128E-12*1E12/1E6;	// vacuum pertivity [W.ps/V²μm] 


#endif // -> #ifdef _DTYPESCONSTSH
