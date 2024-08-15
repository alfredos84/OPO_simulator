#ifndef _COMMON_H
#define _COMMON_H


double seconds() {
    static auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - start;
    return elapsed.count();
}


void TimingCode( double iElaps)
{
	if ( iElaps < 60. ){std::cout << "\n\nTime elapsed " <<  iElaps << " seconds\n\n " << std::endl;}
	else if ( iElaps >= 60. and iElaps < 3600. ){std::cout << "\n\nTime elapsed " <<  iElaps/60 << " minutes\n\n " << std::endl;}
	else{std::cout << "\n\nTime elapsed " <<  iElaps/3600 << " hours\n\n " << std::endl;}
	
	return ;
}


// Linear spacing for time vectors
template<typename T>
void linspace( T& Vec, real_t xmin, real_t xmax)
{
    uint size = Vec.size();
	for (int i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return ;
}


template<typename T>
T linspace( real_t xmin, real_t xmax, uint size)
{
    T Vec (size);
	for (int i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return Vec ;
}


struct PrintComplex
{
    
    void operator()(const complex_t& x) const
    {
        printf("(%f, %f)\n", x.real(), x.imag());
    }
};


struct PrintReal
{
    
    void operator()(const real_t& x) const
    {
        printf("%f\n", x);
    }
};


void runOPO_status( uint r, uint print_each)
{
    if( (r%print_each == 0) or (r == NRT-1) )
        std::cout << "# Round trip: " << r << " - Completed " << r*100/NRT << "%" << "\t\r" << std::flush;
    return ;
}

template<typename T>
void printVarOnScreen( std::string text,  T var)
{
    std::cout << text << var << std::endl;
    return ;
}


#endif // _COMMON_H