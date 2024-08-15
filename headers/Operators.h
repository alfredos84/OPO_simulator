/*---------------------------------------------------------------------------*/
// * This file contains a set of overloaded operators to deal with complex numbers.
/*---------------------------------------------------------------------------*/


#ifndef _OPERATORS
#define _OPERATORS

#pragma once


/////////////////////////////////////     OPERATORS     ////////////////////////////////////////
inline complex_t  operator+(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a   + b[0][0];
	c[0][1] =     + b[0][1];
	
	return c;
}


inline complex_t  operator+(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c[0][0] = a   + b[0][0];
	c[0][1] =     + b[0][1];
	
	return c;
}


inline complex_t  operator+(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a[0][0] + b[0][0];
	c[0][1] = a[0][1] + b[0][1];
	
	return c;
}

inline complex_t  operator-(const complex_t &a) {
	
	complex_t c;    
	c[0][0] = -a[0][0];
	c[0][1] = -a[0][1];
	
	return c;
}

inline complex_t  operator-(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a   - b[0][0];
	c[0][1] =     - b[0][1];
	
	return c;
}


inline complex_t  operator-(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c[0][0] =  b[0][0] - a ;
	c[0][1] =  b[0][1] ;
	
	return c;
}


inline complex_t  operator-(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a[0][0] - b[0][0];
	c[0][1] = a[0][1] - b[0][1];
	
	return c;
}


inline complex_t  operator*(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a * b[0][0] ;
	c[0][1] = a * b[0][1] ;
	
	return c;
}


inline complex_t  operator*(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c[0][0] = a * b[0][0] ;
	c[0][1] = a * b[0][1] ;
	
	return c;
}


inline complex_t  operator*(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c[0][0] = a[0][0] * b[0][0] - a[0][1] * b[0][1] ;
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[0][0] ;
	
	return c;
}


inline complex_t  operator/(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c[0][0] = b[0][0] / a ;
	c[0][1] = b[0][1] / a ;
	
	return c;
}


inline complex_t operator/(const real_t& a, const complex_t& b) {

	real_t denominator = b[0][0] * b[0][0] + b[0][1] * b[0][1];
	complex_t c;
	c[0][0] = (+a * b[0][0]) / denominator;
	c[0][1] = (-a * b[0][1]) / denominator;

	return c;
}


inline complex_t operator/(const complex_t& a, const complex_t& b) {

	real_t denominator = b[0][0] * b[0][0] + b[0][1] * b[0][1];
	complex_t c;
	c[0][0] = (a[0][0] * b[0][0] + a[0][1] * b[0][1]) / denominator;
	c[0][1] = (a[0][1] * b[0][0] - a[0][0] * b[0][1]) / denominator;

	return c;
}


/////////////////////////////////////////////////////////////////////////////////////////////////

//* Complex exponential e^(i*a) */
complex_t CpxExp (real_t a)
{
	complex_t b;
	b[0][0] = cosf(a) ;	b[0][1] = sinf(a) ;
	
	return b;
}


//* Complex exponential e^(a+i*b) */
complex_t CpxExp (complex_t a)
{
	complex_t b;
	b[0][0] = expf(a[0][0])*cosf(a[0][1]) ;	b[0][1] = expf(a[0][0])*sinf(a[0][1]) ;
	
	return b;
}


//* Complex conjugate */
complex_t CpxConj (complex_t a)
{
	complex_t b;
	b[0][0] = +a[0][0] ; b[0][1] = -a[0][1] ;
	
	return b;
}


//* Complex absolute value  */
real_t CpxAbs (complex_t a)
{
	real_t b;
	b = sqrtf(a[0][0]*a[0][0] + a[0][1]*a[0][1]);
	
	return b;
}


//* Complex square absolute value */
real_t CpxAbs2 (complex_t a)
{
	real_t b;
	b = a[0][0]*a[0][0] + a[0][1]*a[0][1];
	
	return b;
}




//////////////////////////////
struct MulRealVecRealScalar // Functor V*s: V real vector and s real scalar
{
	real_t s;
	MulRealVecRealScalar(real_t N) {s = N;};
	
	real_t operator()(real_t x1)
	{
		return x1*s;
	}
};


struct MulCpxVecRealScalar // This functor performs the a*A, a real
{
	real_t a;
	MulCpxVecRealScalar(real_t A)  {a=A;};

	
	complex_t operator()(complex_t V1)
	{
		complex_t result; 
		result[0][0] = a*V1[0][0]; result[0][1] = a*V1[0][1]; 
		return result;
	}
};


struct MulCpxVecCpxScalar // This functor scales by N, with N as a complex constant
{
	complex_t Norm;
	MulCpxVecCpxScalar(complex_t N) {Norm = N;};
	
	complex_t operator()(complex_t V1)
	{
		return make_float2(V1[0][0]*Norm[0][0] - V1[0][1]*Norm[0][1],
						   V1[0][0]*Norm[0][1] + V1[0][1]*Norm[0][0]);
		
	}
};


// //////////////////////////////
// // Returns A *= c -> A = c*A, with A real vector and c real constant
rVech_t operator*=(rVech_t &rhs, const real_t realscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulRealVecRealScalar(realscalar));
    return rhs;
}


rVech_t operator*=(rVech_t &rhs, const real_t realscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulRealVecRealScalar(realscalar));
    return rhs;
}


// Returns A *= c -> A = c*A, with A complex vector and c real constant
cVech_t operator*=(cVech_t &rhs, const real_t realscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecRealScalar(realscalar));
    return rhs;
}


cVech_t operator*=(cVech_t &rhs, const real_t realscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecRealScalar(realscalar));
    return rhs;
}


// Returns B = c*A, A complex vector and c complex constant
cVech_t operator*(cVech_t &rhs, const complex_t complexscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecCpxScalar(complexscalar));
    return rhs;
}


cVech_t  operator*(cVech_t &rhs, const complex_t complexscalar) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecCpxScalar(complexscalar));
    return rhs;
}

cVech_t operator*(const complex_t complexscalar, cVech_t &rhs) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecCpxScalar(complexscalar));
    return rhs;
}


cVech_t operator*(const complex_t complexscalar, cVech_t &rhs) {
    std::transform(rhs.begin(), rhs.end(),
                      rhs.begin(), MulCpxVecCpxScalar(complexscalar));
    return rhs;
}
//////////////////////////////


#endif // -> #ifdef _OPERATORS