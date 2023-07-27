# cpuOPO package

is a C-based toolkit for simulating optical parametric oscillators using the coupled-wave equations (CWEs) that well-describe the three wave mixing (TWM) processes in a second-order nonlinear media.

The provided software implements a solver for the CWEs including dispersion terms, linear absorption and intracavity element if they are required. It also includes flags to solve nanosecond or continuous wave time regimes. However, the user is free to incorporate picosecond or femtosecond regimes by making the proper corrections.

This code is useful for simulations based on three-wave mixing proccesses such as optical parametric oscillators (OPOs).
It solves the coupled-wave equations (CWEs) for signal, idler and pump using a parallel computing scheme based on Split-Step Fourier Method (SSFM).

## Setup and execution

To run simulations using the package clone this project typing in a terminal
```
git clone https://github.com/alfredos84/cpuOPO.git
```
Once the project was cloned, the user will find a parent folder `OPO` containing the folder `src`that in turn contains the header folder  `headers` with all the required `<files.h>`. The bash file `OPO.sh` used to compile and execute the package by passing several simulations parameters.

### Bash file `src/OPO.sh`

The bash file is mainly used to massively perform simulations by passing the main file different parameters such as pump power, cavity detuning, etc. Before starting the user has to allow the system execute the bash file. To do that type in the terminal
```
chmod 777 OPO.sh # enable permissions for execution
```

Finally, to execute the file execute the following command line
```
./OPO.sh         # execute the files
```

In the `OPO.sh` file you will find the command line for the compilation:
```
nvcc OPO.cu -D<REGIME> -DCW_OPO -DPPLN -I /usr/local/include -L /usr/local/include/fftwf3 -lfftw3 -lfftw3f -lm -o OPO
```
where the preprocessor variable `<REGIME>` could be either `CW_OPO` or `NS_OPO`. This compiles the package using two coupled-wave equations. If the user wants to use three coupled-wave equations, add the additional preprocessor variable `-DTHREE_EQS` in the compilation line. 
Notice here `nvcc` compiler is used and the main file extension is `.cu` instead of traditional `gcc` and `.c`, respectively, since this code was implemented to compare the computational speed up with its CUDA counterpart. Please visit my project `cuOPO` for users interested in parallel computing.

Finally, the execution is done using the command line in the `OPO.sh` file is
```
./OPO $L $TEMP $GRPER $N $R $DELTAS $GDD $TOD $U $MODDEP $FREQMOD | tee -a $FILE
```
where `$ARGx` and others are variables externaly passed to the main file `cuOPO.cu`. It was written in this way to make easier performing simulations massively.

### Outputs

This package returns a set of `.dat` files with the signal, idler and pump electric fields, separated into real and imaginary parts. It also returns time and frequency vectors

### Contact me
For any questions or queries, do not hesitate to contact the developer by writing to alfredo.daniel.sanchez@gmail.com
