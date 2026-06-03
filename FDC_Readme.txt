#####################################################################
#
# Full Diffusion Code (FDC) README File
#
#
#
#
# Created April 5 2016 by Hui Zhu 
# Revised Jan 1 2017 by Hui Zhu
# Revised Feb 21 2017 By Hui Zhu

#####################################################################
#
#Code Description
#
#####################################################################

Full Diffusion Code is developed to calculated quasi-linear bounce-averaged diffusion coefficients for right-hand whistler mode waves in pure-proton plasma in dipole field.
It is MPI parallel architecture written by Fortran language. So We strongly recommend to run it on supercomputer rathan than personal computer.



#####################################################################
#
# Code Organization
#
#####################################################################

/Src		- Fortran source files

	B_field.f90
	Bessel_FirstKind.f90
	BesselFunction_firstkind.f90
	B_field.f90
	Bw_wave.f90
	D_ba.f90
	D_local_n.f90
	ElectronDensity.f90
	errorfuction.f90
	FDC_MPI.f90
	Freq_Spectral.f90
	f_function.f90
	get_Loss_Cone.f90
	int_D_loc.f90
	Integration_Lambda.f90
	int_N_omega.f90
	Lambda_mirror_new.f90
	Local_alpha.f90
	Module_Constants.f90
	Module_Variables.f90
	MPI_Parallel.f90
	ResonanceOrder.f90
	Roots_finder.f90
	Sum_N.f90
	Theta_wave.f90
	WeightFunction.f90
	zrhqr1.f90
	zrhqr2.f90
	zrhqr3.f90
/Plots		-Matlab plotting scripts


/Inputs		-Input file examples
	TypicalHiss.txt
	TypicalNightsideChorus.txt
	TypicalDaysideChorus.txt
	FaTypicalHiss.txt
        FaTypicalNightsideChorus.txt
        FaTypicalDaysideChorus.txt
	KS14Hiss.txt
	JupiterChorus.txt
	2016DaysideChorus.txt
	2016NightsideChorus.txt
        2016Hiss.txt

/Execute	-submission file

	submit_yellowstone.sh

#####################################################################
#
# How to compile and run it
#
#####################################################################

UCLA Hoffman2 Cluster

	1. copy the whole folder into the path you want
	
	2. goto Src path

	3. compile it with two commands below
	
		mpiifort Module_Constants.f90 Module_Variables.f90 -llapack -lblas -O3 -c 
	
		mpiifort *.f90 -llapack -lblas -O3 -o ../Execute/FDC
	
	4. goto Execute path and cp a proper input example file from Inputs path to current path 

	5. properly revise input file according to your need

	6. type 'intelmpi.q' to submit job

NCAR Yellowstone Supercomputer

        1. copy the whole folder into the path you want

        2. goto Src path

        3. compile it with three commands below
	
		module load mkl

                mpif90 Module_Constants.f90 Module_Variables.f90 -c

                mpif90 *.f90 -o ../Execute/FDC

        4. goto Execute path and cp a proper input example file from Inputs path to current path  

        5. properly revise input file according to your need

	6. revise dummyinputfile into the name of inout file you want in submit_yellowstone.sh file
        
	7. type 'bsub<submit.sh<submit_yellowstone.sh' to submit job

GFZ
	GFZ mp machines at our section:

	1. At first, download the FDC code from our gitlab.

	2. If you do not need to modify the input file, then we step to compile it, which is shown in step 4.

	3. If you need to modify the input file or want to get familiar with the input file, please see the README file for FDC.

	4. Compile FDC in server, because later it will run in the server, so ssh mag14(or 15, 16)

	5. When we compile it, we need the following command, and do it in the src folder which contains the .f90 files
       		mpif90 Module_Constants.f90 Module_Variables.f90 -llapack -lblas -O3 -c
       		mpif90 *.f90 -o main -I/usr/lib/x86_64-linux-gnu/openmpi/include -ffree-line-length-none -L/usr/lib/ -l:libblas.so.3 -l:liblapack.so.3

	6. Then we run it
       		mpirun -np 14 ./main ../Inputs/OA2015_Night.txt

       		The .txt file is the input file and the number after -np specifies how many nodes do you need for your run.

	7. After the run is finished. You will find two new generated files in the src directory.
	The new .plt file is the result of calculation; the new .txt file is the logfile.

	GFZ cluster:
	1. At first, apply an accound in GFZ server
	2. ssh glic
	3. compile:
		get into the folder Src
		module load iimpi
		module load OpenBLAS FlexiBLAS
		make all
	4. run:
		get into Excecute folder
		edit the input file part in the.sh file 
	5. sbatch name of .sh file (including .sh)
	6. Use sacct -S YEAR-MO-DA -u username to check if the code is running; you can also write -S teatime to check jobs submitted after 4 pm on the current day	

Using the Makefile

	You can compile the code using the makefile located in the Src folder. 
	There are several compilation options available, as described below. 
	To compile the code, simply use the make command followed by the specified options:

	Compile Options

		- Compile with mpiifort: 
			make all
		- Compile with mpif90: 
			make mpif90
		- Compile with MKL library support using mpiifort: 
			make mkl
		- Clean the build (delete *.mod, *.o, and the compiled program): 
			make clean

On Ubuntu

	It is necessary to install an MPI library on your laptop, even if you do not plan to use it for parallel execution.
	To install the required MPI library, run the following commands:

		sudo apt-get update
		sudo apt-get install libopenmpi-dev openmpi-bin
		
	Make sure that Lapack and BLAS are installed
		sudo apt-get install libblas-dev liblapack-dev

	Compile using the mpif90 option: 
		make mpif90

	FDC cannot be executed on a single process. Use minimum 2 processed to run the code:
		mpiexec -n 2 ./FDC_mpi	
	
#####################################################################
#
# How to initialize the Input text file
#
#####################################################################

Line 12: specify the name of Output file.
Line 15: specify the planet.
		Planet Options:		Earth : B0 = 2.997527836235779e-5 Tesla
					Jupiter : B0 = 4.2818e-4 Tesla
					Dipole : B0 is determined by the following B0, and the configuration is dipolar field.
Line 18: specify the magnetic intensity. The intensity works only if Planet option is specified as Dipole.
Line 21: specify the points of pitch angle, energy and Lshell.
Line 24: specify the minimum and maximum of Lshell.
Line 27: specify the minimum and maximum of Energy.
Line 30: specify the minimum and maximum of pitch angle.
Line 33: specify the values of Kp10 and MLT.
Line 36: specify the resonance order type (NresType), minimum and maximum of resonance order. The value works only if NresType is specified as constant.
		NresType Options:	Constant : minimum and maximum are determined by input values.
				        MR12/Mr12/mR12/mr12 : minimum and maximum are determined by Equation 22 in Mourenas and Ripoll 2012.
Line 39: specify the electron density model (DenModel), number density value and index of density latitudinal variation (DenLatINdex). The number denisyt value works only if DenModel is specified as constant.
		DenModel Options:	Constant : electron density is determined by input value.
					CA12/Ca12/cA12/ca12 : electron density is determined by plasmasphere model in Carpenter and Anderson 1992.
					SHEELY01/Sheely01/sheely01 : electron density is determined by plasmatrough model in Sheely et al. 2001.
					DENTON06/Denton06/denton06 : electron density is determined by plasmasphere model in Denton et al. 2004 and 2006.
		DenLatINdex forumla:	N = Neq * cos(latitude)^DenLatINdex , Neq is equatorial density and N is density at particular latitude.
			DenLatINdex = 0 means that the electron density is constant along the latitude.
Line 42: specify the frequency spectral distribution (FreqDistribution).
		FreqDistribution option:	Gaussian/GAUSSIAN : Gaussian distribution
						2016Hiss/2016HISS : Spectral is determined in Li et al. 2015.
						2017Chorus/2017CHORUS : spectral is determined 	

Line 45: specify the frequency unit (FreqUnit), frequency center, frequency width, frequency lower cutoff and frequency upper cutoff.
		FreqUnit Options:	HZ/Hz/hZ/hz : Input frequency center, frequency width, frequency lower cutoff and frequency upper cutoff are specified with unit Hz.
					RATIO/Ratio/ratio : Input frequency center, frequency width, frequency lower cutoff and frequency upper cutoff are specified with 
							    the ratio to equatorial electron gyrofrequency.
					2017Chorus/2017CHORUS : Input frequency lower cutoff and frequency upper cutoff are specifeid
Line 48: specify the points of latitude, normal angle and normalized factor N.
Line 51: specify the minimum and maximum of wave latitude range
Line 54: specify the wave amplitude type (BwType) and amplitude value. The amplitude value works only if BwType is specified as constant.
		BwType Options:		Constant : wave amplitude is determined by input value.
					TypicalDaysideChorus : typical dayside chorus amplitude in Kim et al. 2011
					TypicalNightsideChorus : typical nightside chorus amplitude in Kim et al. 2011
					TypicalHiss : typical plasmaspheric hiss amplitude in Kim et al. 2011
					KS14NightsideHiss : plasmaspheric hiss amplitude model in Orlova et al. 2014
					NbbLBC : Lower band Chorus amplitude model in Binbin et al. 2011
					NbbUBC : Upper band Chorus amplitude model in Binbin et al. 2011
					2016DaysideChorus : dayside chorus in Hui et al. 2017
					2016NightsideChorus : nightside chorus in Hui et al. 2017
					2017Chorus/2017CHORUS : chorus wave for all MLTs
Line 57: specify the normal angle model type (ThetaType), normal angle center for first Gaussian, width for first Gaussian, center for second Gaussian, width for second Gaussian, lower cutoff, upper cutoff and ratio between first Gaussian and second Gaussian. The values work only if ThetaType is specified as constant
		ThetaType Options:	Constant : normal angle center, width, lower cutoff and upper cutoff are determined by input values.
					ObliqueHiss : normal angle center, width, lower cutoff and upper cutoff vary with latitude see Orlova et al. 2014.
					OA2015: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Agapitov et al.2015
					NbbLBC: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Binbin et al.2011
					NbbUBC: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Binbin et al.2011
#####################################################################
#
# Reference
#
#####################################################################

Carpenter and Anderson (1992), An ISEE/whistler model of equatorial electron density in the magnetosphere, J. Geophys. Res. 
Sheely et al. (2001), An empirical plasmasphere and trough density model CRRES observations, J. Geophys. Res. 
Denton et al. (2004), Electron density in the magnetosphere, J. Geophys. Res.
Denton et al. (2006), Distribution of density along magnetospheric field lines, J. Geophys. Res.
Kim et al. (2011), Understanding the dynamic evolution of the relativistic electron slot region including radial and pitch angle diffusion, J. Geophys. Res.
Mourenas and Ripoll (2012), Analytical estimates of quasi-linear diffusion coefficients and electron lifetimes in the inner radiation belt, J. Geophys. Res.
Orlova et al. (2014), Activity-dependent global model of electron loss inside the plasmasphere, Geophys. Res. Lett.
BinBin et al. (2011),Resonant scattering of plasma sheet electrons leading to diffuse auroral precipitation:2. Evaluation for whistler mode chorus waves, J. Geophys. Res.
Li et al. (2015), Statiscal properties of plasmaspheric hiss derived from Van Allen Probes data and their effects on radiation belt electron dynamics, J. Geophys. Res.
Zhu et al. (2017), Long-term relativistic electron radiation-belt simulations with VERB code: Application of new diffusion coefficient parameterizations, J. Geophys. Res., submitted
#####################################################################	
#
#Contact
#
#####################################################################

Please contact Hui Zhu		zhuhui328@gmail.com 
	       Yuri Shprits	shprits@igpp.ucla.edu


