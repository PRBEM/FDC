## User defined variables

`Line 12`: specify the name of Output file.

`Line 15`: specify the planet.

*Planet Options*:		
- `Earth` : B0 = 2.997527836235779e-5 Tesla
- `Jupiter` : B0 = 4.2818e-4 Tesla
- `Dipole` : B0 is determined by the following B0, and the configuration is dipolar field.

`Line 18`: specify the magnetic intensity. 

>The intensity works only if Planet option is specified as Dipole.

`Line 21`: specify the points of pitch angle, energy and Lshell.

`Line 24`: specify the minimum and maximum of Lshell.

`Line 27`: specify the minimum and maximum of Energy.

`Line 30`: specify the minimum and maximum of pitch angle.

`Line 33`: specify the values of Kp10 and MLT.

`Line 36`: specify the resonance order type (NresType), minimum and maximum of resonance order. 
>The value works only if NresType is specified as constant.

*NresType Options*:	
- `Constant` : minimum and maximum are determined by input values.
- `MR12/Mr12/mR12/mr12`: minimum and maximum are determined by Equation 22 in Mourenas and Ripoll 2012.

`Line 39`: specify the electron density model (DenModel), number density value and index of density latitudinal variation (DenLatINdex). 
>The number denisyt value works only if DenModel is specified as constant.

*DenModel Options*:	
- `Constant` : electron density is determined by input value.
- `CA12/Ca12/cA12/ca12` : electron density is determined by plasmasphere model in Carpenter and Anderson 1992.
- `SHEELY01/Sheely01/sheely01` : electron density is determined by plasmatrough model in Sheely et al. 2001.
- `DENTON06/Denton06/denton06` : electron density is determined by plasmasphere model in Denton et al. 2004 and 2006.
- `DenLatINdex forumla`:	$N = Neq \times(cos(latitude)^{DenLatINdex})$ , Neq is equatorial density and N is density at particular latitude. DenLatINdex = 0 means that the electron density is constant along the latitude.

`Line 42`: specify the frequency spectral distribution (FreqDistribution).
		
*FreqDistribution option*:	
- `Gaussian/GAUSSIAN` : Gaussian distribution
- `2016Hiss/2016HISS` : Spectral is determined in Li et al. 2015.
- `2017Chorus/2017CHORUS` : spectral is determined 	

`Line 45`: specify the frequency unit (FreqUnit), frequency center, frequency width, frequency lower cutoff and frequency upper cutoff.
		
*FreqUnit Options*:	
- `HZ/Hz/hZ/hz` : Input frequency center, frequency width, frequency lower cutoff and frequency upper cutoff are specified with unit Hz.
- `RATIO/Ratio/ratio` : Input frequency center, frequency width, frequency lower cutoff and frequency upper cutoff are specified with the ratio to equatorial electron gyrofrequency.
- `2017Chorus/2017CHORUS` : Input frequency lower cutoff and frequency upper cutoff are specifeid

`Line 48`: specify the points of latitude, normal angle and normalized factor N.

`Line 51`: specify the minimum and maximum of wave latitude range

`Line 54`: specify the wave amplitude type (BwType) and amplitude value. 
>The amplitude value works only if BwType is specified as constant.
		
*BwType Options*:		
- `Constant` : wave amplitude is determined by input value.
- `TypicalDaysideChorus` : typical dayside chorus amplitude in Kim et al. 2011
- `TypicalNightsideChorus` : typical nightside chorus amplitude in Kim et al. 2011
- `TypicalHiss` : typical plasmaspheric hiss amplitude in Kim et al. 2011
- `KS14NightsideHiss` : plasmaspheric hiss amplitude model in Orlova et al. 2014
- `NbbLBC` : Lower band Chorus amplitude model in Binbin et al. 2011
- `NbbUBC` : Upper band Chorus amplitude model in Binbin et al. 2011
- `2016DaysideChorus` : dayside chorus in Hui et al. 2017
- `2016NightsideChorus` : nightside chorus in Hui et al. 2017
- `2017Chorus/2017CHORUS` : chorus wave for all MLTs

`Line 57`: specify the normal angle model type (ThetaType), normal angle center for first Gaussian, width for first Gaussian, center for second Gaussian, width for second Gaussian, lower cutoff, upper cutoff and ratio between first Gaussian and second Gaussian. 
>The values work only if ThetaType is specified as constant
		
*ThetaType Options*:	
- `Constant` : normal angle center, width, lower cutoff and upper cutoff are determined by input values.
- `ObliqueHiss` : normal angle center, width, lower cutoff and upper cutoff vary with latitude see Orlova et al. 2014.
- `OA2015`: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Agapitov et al.2015
- `NbbLBC`: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Binbin et al.2011
- `NbbUBC`: normal angle center, width, lower cutoff and upper cutoff vary with latitude see Binbin et al.2011

