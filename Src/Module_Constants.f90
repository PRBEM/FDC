module constants
implicit none

real*8,parameter::PI=3.1415926
real*8,parameter::RAD=PI/180.0
real*8,parameter::me=9.1093897e-31,mi=1.6726231e-27,c=2.99792458e8,C0=8.84194128288e-12,qe=1.60217733e-19,sqe=-1
integer,parameter::unit_input=10,unit_log=20,unit_output=114

character(len=30)::PT  !!planet
REAL*8::B0

character(len=100)::Output_file,log_file

integer::NL,NE,NPA
real*8::PA_min,PA_max,E_min,E_max,L_min,L_max

real*8::Kp10,Kp,MLT

character(len=30)::NresType
integer::nres_min,nres_max

character(len=30)::DenModel
real*8::DenValue,DenLatIndex

character(len=30)::FreqUnit
character(len=30)::FreqDistribution
real*8::omega_m_per,d_omega_per,o_lc,o_uc

integer::nint,nint_zero_res,nint_X,nint_N,nint_N_res_cone

real*8::lam_min,lam_max

character(len=30)::BwType
real*8::Bw_Input

character(len=30)::ThetaType
real*8::Theta_Input(7)

real::time_start,time_finish

end module constants
