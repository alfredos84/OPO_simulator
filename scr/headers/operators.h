/*---------------------------------------------------------------------------*/
// * This file contains a set of operatos useful for the execution of the main file.
/*---------------------------------------------------------------------------*/


#ifndef _OPERATORSCUH
#define _OPERATORSCUH

#pragma once

/** Sinc funcion: sin(x)/x */
real_t sinc(  real_t x  )
{
	// SINC function
	if (x == 0){return 1.0;} else{ return sinf(x)/x;}
}

/////////////////////////////////////     OPERATORS     ////////////////////////////////////////
inline complex_t  operator+(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a   + b.x;
	c.y =     + b.y;
	
	return c;
}


inline complex_t  operator+(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c.x = a   + b.x;
	c.y =     + b.y;
	
	return c;
}


inline complex_t  operator+(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	
	return c;
}


inline complex_t  operator-(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a   - b.x;
	c.y =     - b.y;
	
	return c;
}


inline complex_t  operator-(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c.x =  b.x - a ;
	c.y =  b.y ;
	
	return c;
}


inline complex_t  operator-(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	
	return c;
}


inline complex_t  operator*(const real_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a * b.x ;
	c.y = a * b.y ;
	
	return c;
}


inline complex_t  operator*(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c.x = a * b.x ;
	c.y = a * b.y ;
	
	return c;
}


inline complex_t  operator*(const complex_t &a, const complex_t &b) {
	
	complex_t c;    
	c.x = a.x * b.x - a.y * b.y ;
	c.y = a.x * b.y + a.y * b.x ;
	
	return c;
}


inline complex_t  operator/(const complex_t &b, const real_t &a) {
	
	complex_t c;    
	c.x = b.x / a ;
	c.y = b.y / a ;
	
	return c;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

complex_t CpxExp (real_t a)
{
	complex_t b;
	b.x = cos(a) ;	b.y = sin(a) ;
	
	return b;
}


complex_t CpxConj (complex_t a)
{
	complex_t b;
	b.x = +a.x ; b.y = -a.y ;
	
	return b;
}


real_t CpxAbs (complex_t a)
{
	real_t b;
	b = sqrtf(a.x*a.x + a.y*a.y);
	
	return b;
}


real_t CpxAbs2 (complex_t a)
{
	real_t b;
	b = a.x*a.x + a.y*a.y;
	
	return b;
}


#endif // -> #ifdef _OPERATORSCUH
