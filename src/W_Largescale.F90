module W_Largescale_mod



use mpp_mod,               only: input_nml_file
use time_manager_mod,      only: time_type
use constants_mod,         only: grav, rdgas, rvgas, pi, cp_air, kappa
use fv_arrays_mod,         only: fv_grid_type
use fms_mod,               only: error_mesg, FATAL, file_exist, open_namelist_file,  &
                                 check_nml_error, mpp_pe, mpp_root_pe, close_file, &
                                 write_version_number, stdlog, mpp_error

use fv_mp_mod,         only: mp_reduce_sum, mp_reduce_min, mp_reduce_max, is_master

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end, &
                                   vert_turb_driver_restart


use diag_manager_mod,      only: register_diag_field, send_data

implicit none


  public ::  W_Largescale_Driver, W_Largescale_init, g0_sum

  logical:: g0_sum_initialized = .false.
!  real, allocatable, dimension(:,:) :: l_area
!  real:: global_area = -1.
  character(len=16) :: mod_name = 'W_Largescale_mod'
!  real, allocatable, dimension(:,:) :: l_area

!  real:: global_area = -1.

!
! uanber:  2016:
! 
! This module is used for parameterizing the large scale vertical
! velocity (W_LargeScale) (or similarly the mass flux) in the GFDL cloud 
! resolving limited-domain model.
! Three methods can be used to do this:
!
! 1- Weak Temperature Gradient (WTG) method: relaxing domain mean virtual
!  temperature towards a prescribed profile over a  constant relaxation time
!  scale (tau_wtg). This method is implemented using a simple finite difference
!  scheme.
!
! 2- Damped Gravity Waves (DGW) method: simulating a single gravity wave doing
!  the actual adjustment of the density anomaly (or virtual temperature
!  anomaly). Anomaly is difference from a prescribed profile as above. This
!  method is implemented using a finite element scheme that solves the linear
!  tridiagonal system. Note: this is advection diffusion equation with variable
!  source term (RHS) and zero advection velocity. It requires two boundary
!  conditions at surface and tropopause W =0.
!
! 3- Spectral Weak Temperature Gradient (SWTG) method: spectral (Fourier)
!  decomposition of the virtual temperature anomaly (from the prescribed profile).
!  sine functions are chosen as bases for the decomposition as they represent
!  wave adjustment in the vertical. This method assumes variable relaxation time
!  with hight. 
!=====================================================!
!=====================================================!
! REMEMBER: 1 is top layer, km is the lowest layer    !
!=====================================================!
!=====================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==================================!
! Weak Temperature Gradient (WTG): !
!
!  W*dtheta_ref/dz  = (theta - theta_ref) /tau
! 
!  In the PBL:
!  W = interpolation from above PBL to zero at the surface
!===========================================================!


!============================!
! Damped Gravity Waves: DGW  !
!============================!

! d^2W/d^2Z = const*(T_v - T_v_ref)/T = RHS
!
! (1/dz) * W_k-1  - (2/dz) * W_k + (1/dz) * W_k+1 = RHS* dz
!      W_k-1  -2* W_k +  W_k+1 = RHS* dz*dz =D
! or:
! 
! B*W = D
!
!              | -2  1  0  0 0 . . .  0 |
!              |  1 -2  1  0 0 . . .  0 |
! B=           |  0  1 -2  1 0 . . .  0 |
!              |  . . . . . .  . . . .  |
!              |  0  0  0  0 .. 1 -2  1 |
!              |  0  0  0  0 . .0  1 -2 |
!
! a = -2 ; b=c= 1
!
! Construction of L and U:
! 
! L = 
!
!      |1  0  0  0 . . . . 0|
!      |e2 1  0  0 . . . . 0|
!      |0  e3 1  0 . . . . 0|
!      |.  .  .  .  .  . . 0|
!      |0  0  0  . e_n-1 1 0| 
!      |0  0  0  .  .  e_n 1|
! U =
!       |f1  b1  0  0  . . .    0 |     
!       |0   f2  b2  0  . . .   0 |
!       |0   0   f3  b3  . .    0 |
!       |.   .   .   .   . . .    |
!       |0   0   0   .f_n-1  b_n-1|
!       |0   0   0  . .  .    fn  |
!
!
! Now: BW=D  ==> (LU)W=D ==> L(UW)=D
!      Ly = D ...(1) solve for y
!      UW = y ...(2) solve for W
!=============================================================================!
!=============================================================================!


!============================!
! Spectral WTG: SWTG         !
!============================!
! A- loop on vertical levels k:
! 1- calculate: d_theta_dz 
! 2- calculate theta_W = (theta_avg-theta_ref)/d_theta_dz
! 3- Calculate Brunt-Vaisala frequency (NN)
! 
! B - vertical integrated average of Brunt-Vaisala frequency N ~= 1e-2 s-1
! 
! C- loop on vertical modes j (1D arrays):
! 1- calculate m_j
! 2- calculate spectral theta and vertically integrate
! 3- calculate spectral tau
!
! D- loop over vertical levels k
!  - loop over vertical modes j:
!  - calculate W_Largescale by summing over dummy j
!
!
!=============================================================================!
!=============================================================================!
!========================
!----- namelist ---------
!========================
! if a variable is in the namelist it should not be in the subroutine argumnet
  logical:: do_W_Largescale      = .true. !.false.
  logical:: q_horiz_adv          = .false. !.false.


  integer ::  id_zfull, id_pfull, id_theta, id_W_Largescale
!uanber
  integer :: id_z_avg, id_theta_avg, id_T_avg, id_q_avg, id_utmp
!uanber
  integer :: tau_WTG = 2. 
  integer :: W_Largescale_flag = 1.
  real ::  wavenumber = 1e-5 
  real ::  damping = 1.157e-5
  real, dimension(51):: theta_ref, moisture_ref  ! 51 is the model levels
  real :: L = 200000.  ! horizontal scale for gravity waves adjustment in meters

  real    :: missing_value = -1.e10


 namelist /W_Largescale_nml/do_W_Largescale, W_Largescale_flag, q_horiz_adv,theta_ref, moisture_ref, & 
                           tau_WTG,  wavenumber, damping, L!, &

                           !zero_winds, tau_zero


 ! write(*,*)  do_W_Largescale
 
 contains

! subroutine  W_Largescale_Driver(npz, is, ie, js, je, ng, nq, pdt, time, hydrostatic, &
!                                 pt, delz, q, pe, delp, peln, t_dt, q_dt, pbltop, W_Largescale3D, z_avg, theta_avg, &
!! WTG Variables                        
!                                 theta, T_avg, q_avg,& ! theta_ref, moisture_ref,& ! tau_WTG,PBL_H,&! K_PBL,&
!! DGW Variables                                  
!                                 T_v, T_v_ref, Tv_ref, Tv_avg,&
!! SWTG Variables
!                                 !L,& ! N, h, c, m, j

!!advection variables                                      
!                                 d_theta_ref_dz,&
!                                 W_dT_dz, W_dq_dz, q_H_dW_dz,&
!! gridstruct
!                                 gridstruct)


! subroutine  W_Largescale_Driver(npz, is, ie, js, je, ng, nq, pdt, time, hydrostatic, &
!                                 pt, pfull, zfull, delz, q, pe, delp, peln, t_dt, q_dt, pbltop, W_Largescale3D,&
!                                 area)



 subroutine  W_Largescale_Driver(npz, is, ie, js, je, ng, nq, pdt, time, hydrostatic, &
                                 pt, pfull, phalf, zfull, zhalf, q, t_dt, q_dt, pbltop, W_Largescale3D,&
                                 area)


!==========================!
! declare varables         !
!==========================!

 integer, INTENT(IN) ::  npz
! integer, INTENT(IN) :: 
 integer, INTENT(IN) :: is, ie, js, je, ng, nq !, nwat
 real   , INTENT(IN) :: pdt  !what is this? physics timestep
! real   , INTENT(IN) :: agrid(is-ng:ie+ng,js-ng:je+ng, 2)
! real   , INTENT(IN) :: ak(npz+1), bk(npz+1)
! real, INTENT(IN):: p_ref
! real, INTENT(IN):: oro(is:ie,js:je)       ! land fraction
 logical, INTENT(IN):: hydrostatic

 type(time_type), intent(in) :: Time
! real, INTENT(IN), optional:: time_total

! real   , INTENT(INOUT) :: u(is-ng:ie+  ng,js-ng:je+1+ng,npz)
! real   , INTENT(INOUT) :: v(is-ng:ie+1+ng,js-ng:je+  ng,npz)
! real, INTENT(INOUT)::  pk (is:ie,js:je,npz+1)
! real, INTENT(INOUT)::  pkz(is:ie,js:je,npz)
! real, INTENT(INOUT)::   pt(is-ng:ie+ng,js-ng:je+ng,npz)!, theta(is-ng:ie+ng,js-ng:je+ng,npz) !theta is potential temperature 
! real, INTENT(INOUT):: delp(is-ng:ie+ng,js-ng:je+ng,npz)
! real, INTENT(INOUT)::    q(is-ng:ie+ng,js-ng:je+ng,npz, nq)  !
! real, INTENT(INOUT)::   pe(is-1:ie+1 ,1:npz+1,js-1:je+1)
 !real, INTENT(INOUT):: peln(is  :ie   ,1:npz+1,js  :je  )
! real, INTENT(INOUT):: delz(is-ng :ie+ng  ,js-ng  :je+ng  ,npz)

 real, INTENT(INOUT)::   pt(is:ie,js:je,npz)!, theta(is-ng:ie+ng,js-ng:je+ng,npz) !theta is potential temperature 
 !real, INTENT(INOUT):: delp(is:ie,js:je,npz)
 real, INTENT(INOUT)::    q(is:ie,js:je,npz,nq)  !
 !real, INTENT(INOUT):: delz(is:ie,js:je,npz)
 real,            intent(in),    dimension(is:ie,js:je,npz) :: pfull, zfull
 real, intent(in), dimension(is:ie,js:je,npz+1) :: phalf, zhalf
! real, intent(inout):: w(is-ng:ie+ng,js-ng:je+ng,npz)
! real, INTENT(INOUT):: sst(is:ie,js:je)
 real, intent(in), dimension(:,:)   :: area

 !type(fv_grid_type), intent(INOUT), target :: gridstruct
! type(fv_grid_type) :: gridstruct
! area => gridstruct%area

! Tendencies:
! real, INTENT(INOUT):: u_dt(is-ng:ie+ng,js-ng:je+ng,npz)
! real, INTENT(INOUT):: v_dt(is-ng:ie+ng,js-ng:je+ng,npz)
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)
! real, INTENT(IN):: ua(is-ng:ie+ng,js-ng:je+ng,npz)
! real, INTENT(IN):: va(is-ng:ie+ng,js-ng:je+ng,npz)
! uanber: surface fluxes
! local
! real, dimension(is:ie,js:je):: flux_t, flux_q, flux_u, flux_v, delm !uanber: surface fluxes of heat, moisture and momentum, respectively.
 !real, dimension(is:ie,js:je,npz):: theta  ! uanber: potential temperature
!================================================!
!uanber: W_Largescale. define variables as local !
!================================================!
 real, dimension(npz) :: W_Largescale, W_dT_dz, W_dq_dz, q_H_dW_dz
 real, dimension(npz) :: d_theta_ref_dz
! real, dimension(npz) , INTENT(INOUT) :: z_avg, theta_avg, T_avg, q_avg !domain average z, theta, T, and q
 ! real, dimension(is-ng:ie+ng,js-ng:je+ng,npz):: theta
 real, dimension(is:ie,js:je,npz):: theta, temp
 real, dimension(is:ie,js:je,npz,nq):: q3
!  real, dimension(npz):: theta_ref, q_ref  ! theta_reference and  q reference (respect.)
 !real, dimension(npz):: d_theta_ref_dz  ! static stability
 real, dimension(npz):: dW_dz !vertical advection of T, q, dW/dz,horizontal advection of q (respect.)
 real, INTENT(IN), dimension(is:ie, js:je) :: pbltop  ! top of the PBL   
!real:: PBL_H  ! pbl hight in (m), should be in the namelist
! real :: tau_WTG  ! relaxation time scale for WTG in hours. Should be in the namelist
 real :: H_tropopause ! tropopause height in m
 integer:: K_PBL!, W_Largescale_flag  !PBL_K is the first model level above the boundary layer.

!  real :: PBL_H
 !real, dimension(is:ie,js:je,npz):: theta!(is-ng:ie+ng,js-ng:je+ng,npz)
 real, parameter:: g= 9.81   !gravity acceleration
! wevecoupling
 real, dimension(is:ie,js:je,npz):: T_v!(is-ng:ie+ng,js-ng:je+ng,npz)
 real, dimension(is:ie,js:je,npz):: T_v_ref!(is-ng:ie+ng,js-ng:je+ng,npz)
 real, dimension(npz):: Tv_ref, Tv_avg

 !type(fv_grid_type) :: gridstruct
 ! local variables for finite element method to solve for DGW
! real :: wavenumber, damping !  parameters for namelist
 real, dimension(npz)::a, b, c, e, f, y,RHS, D
! real, intent(IN):: wavenumber, damping !  parameters for namelist
! real, dimension(is:ie, js:je, npz):: dz

real const
!= wavenumber*wavenumber*g/damping
! logical used
! logical :: q_horiz_adv
!=================================================!
!     for spectral WTG                            !
!=================================================!
real, dimension(npz) :: theta_W, NN, int_N, int_theta ! integral N and integral theta
real :: N
integer , parameter:: JJ=10  ! number of vertical modes
!J=10
real, dimension(JJ) :: theta_mode, m, tau_mode

!=============================================================!
! Uanber New_Cool and cloud fraction:
! Local
 !real, dimension(is:ie,js:je):: flux_t, flux_q, flux_u, flux_v, delm
! logical:: phys_hydrostatic = .true.
! real, parameter:: f1 = 2./3.
! real, parameter:: f2 = 1./3.
! real, dimension(is:ie,js:je,npz):: u3, v3, t3, p3, dz, zfull !remove dz
 real, dimension(is:ie,js:je,npz) :: p3, dz
! real, dimension(is:ie,js:je,npz) :: zfull, pfull
! real, dimension(is:ie,js:je,npz+1):: zhalf
! real, dimension(is:ie,js:je,npz,nq):: q3
! real, dimension(is:ie,js:je):: rain, snow, ice, graup, land, mu
 real, dimension(is:ie,js:je):: ps !, qs, rho, clouds
! real, dimension(is:ie,js:je):: olr, lwu, lwd, sw_surf, wet_t, net_rad
 real, INTENT(OUT), dimension(is:ie,js:je,npz):: W_Largescale3D
! Flux diag:
! real, dimension(is:ie,js:je):: rflux, qflux
 real, dimension(is:ie,npz):: den   ! density
 real :: tvm, rrg, zvir
! real, dimension(is:ie,npz):: t_dt_rad
! real, dimension(npz):: utmp, vtmp
! real:: sday, rrg, tvm, olrm, swab, sstm, clds, hflux1, hflux2, hflux3, precip
! real:: tmp, cooling, heating
! real:: fac_sm, rate_w, rate_u, rate_v, rate_t, rate_q
! real:: prec
 integer  i,j,k, km, iq, k_mp, top, srf ! uanber: top is the height of the tropopause. srf =km-1
! integer  isd, ied, jsd, jed
 integer  seconds, days
 logical print_diag, used

!============================!
! global average variables   !
!============================!

real, dimension(npz):: z_avg, theta_avg, T_avg, q_avg, p_avg, dz_avg, dp_avg, tdt_avg, theta_diff
real :: pbltop_avg

!  if (.not. sim_phys_initialized) call fv_phys_init(nwat)

!   call get_time (time, seconds, days)

!   if ( mod(seconds, 3600*print_freq)==0 ) then
!        print_diag = .true.
!   else
!        print_diag = .false.
!   endif


   zvir = rvgas/rdgas - 1.
   rrg  = rdgas / grav
   km = npz
   !isd = is-ng;   ied = ie + ng
   !jsd = js-ng;   jed = je + ng


!---------------------------------------------!
!   theta_ref = (/1422.34899902,   988.44714355,   837.76971436,   756.19171143,&
!         704.89904785,   659.75274658,   624.69958496,   594.47735596,&
!         567.9487915 ,   543.65704346,   521.18688965,   500.29034424,&
!         480.86959839,   462.89511108,   446.25180054,   430.93978882,&
!         416.78634644,   403.6416626 ,   391.50469971,   380.15713501,&
!         369.61151123,   359.7885437 ,   350.62435913,   342.14447021,&
!         335.53790283,   332.26791382,   330.10726929,   328.77716064,&
!         327.11672974,   325.54666138,   323.65499878,   321.67593384,&
!         319.49676514,   317.1675415 ,   314.73623657,   312.2074585 ,&
!         309.72235107,   307.33660889,   305.25241089,   303.37039185,&
!         301.9006958 ,   300.50131226,   299.49298096,   298.48794556,&
!        297.9647522 ,   297.26559448,   297.08786011,   296.51623535,&
!        296.42642212,   295.98757935,   296.66699219/)

! theta_ref = (/1624.05685505,  1245.53470078,   996.37626489,   861.74902814,&
!         779.87797585,   721.35924631,   675.56622485,   636.40649359,&
!         601.37753688,   569.44904124,   539.84040639,   512.78841644,&
!         486.86093987,   461.74687325,   436.50710646,   413.71090084,&
!         393.66352812,   377.4016675 ,   364.17051294,   353.84999995,&
!         347.13724841,   344.50694725,   342.66770188,   341.34019297,&
!         339.74673707,   338.23427005,   336.651982  ,   335.03414179,&
!         333.26415978,   331.39221901,   329.28123401,   327.02147638,&
!         324.58662257,   322.11080911,   319.58346384,   317.06066052,&
!         314.62313579,   312.36285691,   310.2831217 ,   308.32340135,&
!         306.47961862,   304.80509869,   303.0293399 ,   300.97815708,&
!         299.18853247,   298.45175422,   298.21276299,   298.10618403,&
!         298.04904548,   298.00837989,   298.01508165/)



! moisture_ref = (/3.51206131e-06,   3.42603539e-06,   3.14272052e-06,&
!         3.04803370e-06,   2.87751459e-06,   2.79866413e-06,&
!         2.62646995e-06,   2.52098857e-06,   2.38344478e-06,&
!         2.33124820e-06,   2.27568466e-06,   2.34541767e-06,&
!         2.51643678e-06,   2.59493663e-06,   2.61891091e-06,&
!         3.10319979e-06,   5.15995907e-06,   7.68771406e-06,&
!         9.47700482e-06,   1.08837930e-05,   1.44016876e-05,&
!         2.40854370e-05,   4.13352318e-05,   6.72906390e-05,&
!         1.10226872e-04,   1.67888211e-04,   2.35562358e-04,&
!         3.27264454e-04,   4.49523475e-04,   5.98036451e-04,&
!         7.79748603e-04,   1.00475212e-03,   1.27316895e-03,&
!         1.58247189e-03,   1.93739706e-03,   2.38717603e-03,&
!         2.88995006e-03,   3.42448172e-03,   3.96322086e-03,&
!         4.53406246e-03,   5.19566052e-03,   6.00345246e-03,&
!         7.29089137e-03,   9.31563322e-03,   1.15839178e-02,&
!         1.27766104e-02,   1.32187102e-02,   1.33745428e-02,&
!         1.34309800e-02,   1.34878587e-02,   1.35725066e-02/)
!
!moisture_ref = (/3.55007046e-06,   3.40918473e-06,   3.07161667e-06,&
!        2.92635309e-06,   3.23223435e-06,   4.13958469e-06,&
!        4.82509040e-06,   5.43248962e-06,   1.12531961e-05,&
!        2.83258796e-05,   4.08864325e-05,   4.99448543e-05,&
!        5.32221165e-05,   5.47802629e-05,   5.52017264e-05,&
!        5.54829276e-05,   5.56736813e-05,   5.60807894e-05,&
!        5.64286747e-05,   5.66506278e-05,   5.88498624e-05,&
!        6.55457843e-05,   7.23589183e-05,   8.06584794e-05,&
!        9.79433680e-05,   1.40152712e-04,   2.13823616e-04,&
!        3.24329652e-04,   4.65922058e-04,   6.43791747e-04,&
!        8.61441076e-04,   1.11722760e-03,   1.41635642e-03,&
!        1.75165629e-03,   2.11879495e-03,   2.51960871e-03,&
!        2.96001858e-03,   3.44486325e-03,   3.92053602e-03,&
!        4.40626033e-03,   4.88584815e-03,   5.36293397e-03,&
!        5.90437558e-03,   6.51787454e-03,   7.29414122e-03,&
!        8.06170143e-03,   8.85904022e-03,   9.31687653e-03,&
!        9.87449382e-03,   1.02995755e-02,   1.45917712e-02/)



! tau_WTG = 2
! wavenumber = 1e-5 !1e-6
! damping = 1.157e-5
! L=200000





!do k=1,km
! write(mpp_pe()+2000,*) k,' ', theta_ref(k)
!enddo
!---------------------------------------------!
!area => gridstruct%area
!---------------------------------------------!
! Compute zfull, zhalf
!do j=js,je
!do i=is,ie
!         zhalf(i,j,npz+1) = 0.
!         ps(i,j) = pe(i,km+1,j)
!      enddo
!   enddo

!  if ( hydrostatic ) then
!   do k=km,1,-1
!      do j=js,je
!         do i=is,ie
!                     tvm = rrg*pt(i,j,k)*(1.+zvir*q(i,j,k,1))
!               dz(i,j,k) = -tvm*(peln(i,k+1,j)-peln(i,k,j))
!               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!            zfull(i,j,k) = zhalf(i,j,k+1) + tvm*(1.-pe(i,k,j)/p3(i,j,k))
!            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
!         enddo
!      enddo
!   enddo
! else
   do k=km,1,-1  !Linjoing: reversed the loop as in hydrostatic case above
     do j=js,je
        do i=is,ie
                dz(i,j,k) = zhalf(i,j,k) - zhalf(i,j,k+1)               
!               dz(i,j,k) = delz(i,j,k)
!               p3(i,j,k) = delp(i,j,k)/(peln(i,k+1,j)-peln(i,k,j))
!            zhalf(i,j,k) = zhalf(i,j,k+1) - dz(i,j,k)
!            zfull(i,j,k) = zhalf(i,j,k+1) - 0.5*dz(i,j,k)
         enddo
      enddo
   enddo
! endif
!---------------------------------------------!
!      calculate potential temperature        !
!---------------------------------------------!

do k=1, km
   do j=js,je
      do i=is,ie
        ! theta(i,j,k)=pt(i,j,k)*(phalf(i,j,km+1)/pfull(i,j,k))**(kappa)  !potential temperature (phalf(i,ikm+1,j) is surface pressure)
          theta(i,j,k)=pt(i,j,k)*( 1.0E5/pfull(i,j,k) )**(kappa)  !potential temperature (phalf(i,ikm+1,j) is surface pressure)
         temp(i,j,k)=pt(i,j,k)  ! temperature
        !  write(mpp_pe()+1000,*) k,' ', temp(is,js,k)
      enddo
    enddo
enddo

 do iq=1,nq
     do k=1,npz
        do j=js,je
           do i=is,ie
               q3(i,j,k,iq) = q(i,j,k,iq) !+ pdt*q_dt(i,j,k,iq) ! no need for q_dt since it will be updated below
            enddo
         enddo
     enddo
  enddo


!!===========================================================!!
!!===========================================================!!
!!   uanber: tendencies from W_Largescale Parameterizartion  !!
!!===========================================================!!
!!===========================================================!!
! debug
!write(mpp_pe()+1000,*) do_W_Largescale

!if (do_W_Largescale) then



!call W_Largescale_Driver(npx, npy, npz, is, ie, js, je, ng, nq, &
!                                          pt, delz, q, pe, delp, peln, W_Largescale, z_avg, theta_avg, &
! WTG Variables                        
!                                          theta, T_avg, q_avg, theta_ref, q_ref, tau_WTG,PBL_H,&! K_PBL,&
! DGW Variables                                  
!                                          T_v, T_v_ref, Tv_ref, Tv_avg, wavenumber, damping,&
! SWTG Variables
                                      !    N, L, h, c, m, j

! advection variables                                      
!                                          d_theta_ref_dz,&
!                                          W_dT_dz, W_dq_dz, q_H_dW_dz,&
! gridstruct
!                                          gridstruct)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! uanber: domain mean averages for zfull, potential temperature theta, temperature pt, and moisture q: !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1, km
!      do j = js, je
!        do i = is, ie
!          if (zfull(i,j,k) .ne. zfull(i,j,k)) then
!              print*,zfull(i,j,k)
!          endif
!        enddo
!      enddo
!      print*,maxval(zfull),minval(zfull)
       
  ! z_avg(k)     = g0_sum(zfull(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   z_avg(k)     = g0_sum(zfull(is,js,k), is, ie, js, je, 1,.true., area)
!   write(mpp_pe()+7000,*) k,' ', z_avg(k)
!   write(mpp_pe()+2000,*) k,' ', zfull(is,is,k)
!   write(mpp_pe()+8000,*) k,' ', pfull(is,is,k)  
!   write(mpp_pe()+5000,*) k,' ', delz(is,is,k)
!   theta_avg(k) = g0_sum(theta(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   theta_avg(k)     = g0_sum(theta(is,js,k), is, ie, js, je, 1,.true., area) 
 
 theta_diff(k) = theta_avg(k) - theta_ref(k)
!  write(mpp_pe()+1000,*) k,' ', theta_diff(k) 
!  write(mpp_pe()+2000,*) k,' ', theta_avg(k)
!   T_avg(k)     = g0_sum(temp(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   T_avg(k)     = g0_sum(temp(is,js,k), is, ie, js, je, 1,.true., area)
!   write(mpp_pe()+2000,*) k,' ', T_avg(k)
!   q_avg(k)     = g0_sum(q(is,js,k,1), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   q_avg(k)     = g0_sum(q(is,js,k,1), is, ie, js, je, 1,.true., area)
 !  write(mpp_pe()+2000,*) k,' ', q_avg(k)
!   dz_avg(k)    = g0_sum(dz(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   dz_avg(k)     = g0_sum(dz(is,js,k), is, ie, js, je, 1,.true., area)
    write(mpp_pe()+1000,*) k,' ', dz_avg(k)
   !dp_avg(k)    = g0_sum(delp(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   
 !  tdt_avg(k)    = g0_sum(t_dt(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)

!   write(mpp_pe()+1000,*) k,' ', tdt_avg(k)
enddo
  ! pbltop_avg    = g0_sum(pbltop(is,js), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   pbltop_avg    = g0_sum(pbltop(is,js), is, ie, js, je, 1,.true., area)
!   write(mpp_pe()+6000,*)  pbltop_avg
!   write(mpp_pe()+6000,*)  pbltop(is,js)

!!!!!!!!!    write(mpp_pe()+2000,*) k,' ', dp_avg(k)


!===================================================!
!  apply to troposphere only, i.e below tropopause  !
!===================================================!
H_tropopause = 15000 !put this in the namelist

do k= km,1, -1
   if (z_avg(k) > H_tropopause) then
       top = k
      exit
   endif
enddo


!write(mpp_pe()+2000,*) top

!=========================!
! initialize W_Largescale !
!=========================!

do k=1, km
   W_Largescale(k)=0
   !write(mpp_pe()+3000,*) k,' ',W_Largescale(k)
enddo

!W_Largescale_flag =3
!write(mpp_pe()+7000,*) W_Largescale_flag

! WTG flag
!========================================!
!  Weak Temperature Gradient (WTG):      !    
!========================================!


if (W_Largescale_flag == 1) then    ! flag for: 1- WTG

!write(mpp_pe()+1000,*) W_Largescale_flag

!write(mpp_pe()+8000,*) 'I am fed up with this bug'

!uanber: calculate the static stability of the reference profile d_theta_ref_dz:
! Calculations are done on full levels because dtheta/dz is on half level, 
! then multiply by 1/2 to bring it onto full (mass) level, where theta, q, etc. are defined, and where W_Largescale
! should be defined

do k = top, km ! always start from the top so that near surface you find value of W_Largescale at top of PBL to interpolate to surface
  !W_Largescale(k)=0.0   ! initializing W_Largescale 
  if (z_avg(k) > pbltop_avg) then ! apply WTG above the PBL, PBL_H is either set by the user (~1000-1500 m)if no PBL scheme or provided interactively by the PBL scheme
      K_PBL = k
     ! if (k==km) then
      if (k==top)  then  ! at the top of the atmosphere
          d_theta_ref_dz(k) = (theta_avg(k)*(1.0+0.608*q_avg(k)) - theta_avg(k+1)*(1.0+0.608*q_avg(k+1))) / (z_avg(k)-z_avg(k+1))
      elseif (k==km) then
          d_theta_ref_dz(k) = (theta_avg(k-1)*(1.0+0.608*q_avg(k-1)) - theta_avg(k)*(1.0+0.608*q_avg(k))) / (z_avg(k-1)-z_avg(k))
      endif

      if (k>top.and.k<km) then ! below the top :calculating dtheta/dz on full level (=0.5* sum of dthat/dz on previous and following half levels)

         d_theta_ref_dz(k) = 0.5*((theta_avg(k-1)*(1.0+0.608*q_avg(k-1)) - theta_avg(k)*(1.0+0.608*q_avg(k))) / (z_avg(k-1)-z_avg(k)) &
                                 + (theta_avg(k)*(1.0+0.608*q_avg(k)) - theta_avg(k+1)*(1.0+0.608*q_avg(k+1))) / (z_avg(k)-z_avg(k+1)))

        ! write(mpp_pe()+5000,*) k, ' ', d_theta_ref_dz(k)
      endif

!uanber: make sure that static stability is bounded:
      if (d_theta_ref_dz(k) <0) then
          d_theta_ref_dz(k) = min(d_theta_ref_dz(k),-1e-3)
      else
          d_theta_ref_dz(k) = max(d_theta_ref_dz(k),1e-3)
      endif
       !   write(mpp_pe()+5000,*) k, ' ', d_theta_ref_dz(k)
   
!uanber: now apply WTG (note W_Largescale is now on full levels)
      W_Largescale(k) =  ( theta_avg(k) *(1.0+0.608*q_avg(k)) - theta_ref(k) *(1.0+0.608*moisture_ref(k)) ) / (tau_WTG*3600.0) / d_theta_ref_dz(k)
! debug 
 !  write(mpp_pe()+7000,*) k, ' ', W_Largescale(k)*100

  else ! in the boundary layer W_largescale is interpolated from the surface (=0) to the top of hte PBL (=K_PBL)
      W_Largescale(k) =  (W_Largescale(K_PBL) - 0.0)* (z_avg(k) - z_avg(km))/(z_avg(K_PBL)-z_avg(km))
      d_theta_ref_dz(k) =0.0 ! mixing in the boundary layer so no vertical gradient of temperature.
  endif

  
 
 !   write(mpp_pe()+8000,*) k, ' ',theta_ref(k)
     write(mpp_pe()+9000,*) k, ' ',W_Largescale(k)

enddo

! correct the sign since dz is negative (top (1) to bottom (km))
!do k=1, km
!   W_Largescale(k)= - W_Largescale(k)
!enddo


endif

!------------!
! end of WTG !
!------------!

!========================================!
!  Dampmed Gravity Waves (DGW) method:   !    
!========================================!


if (W_Largescale_flag == 2) then

!--------------------------------------------------!
! uanber: calculate virtual potential temperature  !
!--------------------------------------------------!

do k=1, km
!do k=top, km
   do j=js,je
      do i=is,ie
         T_v(i,j,k)=theta(i,j,k)*(1.0+0.608*q(i,j,k,1)) ! of the model's temp.
        ! T_v_ref(i,j,k) = theta_ref(k) *( (ps(i,j)/p3(i,j,k))**(kappa) ) * (1.0+0.608*q_ref(k))  ! of the reference profile
         T_v_ref(i,j,k) = theta_ref(k) * (1.0+0.608*moisture_ref(k))
          
       enddo
    enddo
enddo


do k=1, km
  ! Tv_avg(k) = g0_sum(T_v(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   Tv_avg(k)     = g0_sum(T_v(is,js,k), is, ie, js, je, 1,.true., area)
   write(mpp_pe()+5000,*) k,' ', Tv_avg(k)
 !  Tv_ref(k) = g0_sum(T_v_ref(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
   Tv_ref(k)     = g0_sum(T_v_ref(is,js,k), is, ie, js, je, 1,.true., area)
   write(mpp_pe()+2000,*) k,' ', Tv_ref(k)
!   dp_avg(k) = g0_sum(delp(:,:,:), is, ie, js, je, 0, gridstruct%area(is:ie, js:je), 1)
enddo


!interior of the domin 
!real, integer :: srf, tt
srf = km-1 !layer above km (~surface) because we will glue boundary condition at km and top
!tt = top+1
!
! construct the RHS
!
const = - wavenumber*wavenumber*g/damping

 do k= top,km
    RHS(k) = const *(Tv_avg(k) - Tv_ref(k))/Tv_avg(k)
    D(k) = RHS(k)* dz_avg(k)*dz_avg(k)
 enddo

! B*W = D
! construction of the symmetric matrix B 
!
do k=top,srf
   a(k)=-2
   b(k)=1
   c(k)=1
enddo



! solve the system with LU decomposition:
!
! determine e_i and f_i for matrices: L and U respectively:
!f(1)=a(1)

f(top)= a(top)

do k=top+1, srf
   e(k) = c(k)/f(k-1)
   f(k) = a(k) - e(k)*b(k-1)
!   write(mpp_pe()+1000,*) k, ' ', f(k)
!   write(mpp_pe()+2000,*) k, ' ', f(k)
enddo

! now solve for y in Ly=D (solve 1):
y(top)=D(top)

do k=top+1, srf
   y(k) = D(k) - e(k)*y(k-1)
enddo

! now solve for w in Uw=y (solve 2):
W_Largescale(srf)=y(srf)/f(srf)

do k=srf-1,top, -1
   W_Largescale(k) = ( y(k) - b(k)*W_Largescale(k+1) )/f(k)
enddo

!BCs:
 W_Largescale(top)=0
 W_Largescale(km)=0



endif


! debug
!do k = top, km
!   write(mpp_pe()+9000,*) k,' ', W_Largescale(k)
!enddo


!------------!
! end of DGW !
!------------!
!======================================================!
!  Spectral Weak Temperature Gradient (SWTG) method:   !    
!======================================================!

if (W_Largescale_flag == 3) then

!1- calculate static stability: dtheta_ref_dz

do k = top, km ! always start from the top so that near surface you find value of W_Largescale at top of PBL to interpolate to surface
      if (k==top)  then  ! at the top of the atmosphere
          d_theta_ref_dz(k) = (theta_avg(k)*(1.0+0.608*q_avg(k)) - theta_avg(k+1)*(1.0+0.608*q_avg(k+1))) / (z_avg(k)-z_avg(k+1))
      elseif (k==km) then
          d_theta_ref_dz(k) = (theta_avg(k-1)*(1.0+0.608*q_avg(k-1)) - theta_avg(k)*(1.0+0.608*q_avg(k))) / (z_avg(k-1)-z_avg(k))
      endif
      if (k>top.and.k<km) then ! below the top :calculating dtheta/dz on full level (=0.5* sum of dthat/dz on previous and following half levels)

         d_theta_ref_dz(k) = 0.5*((theta_avg(k-1)*(1.0+0.608*q_avg(k-1)) - theta_avg(k)*(1.0+0.608*q_avg(k))) / (z_avg(k-1)-z_avg(k)) &
                                 + (theta_avg(k)*(1.0+0.608*q_avg(k)) - theta_avg(k+1)*(1.0+0.608*q_avg(k+1))) / (z_avg(k)-z_avg(k+1)))
      endif

     !uanber: make sure that static stability is bounded:
      if (d_theta_ref_dz(k) <0) then
          d_theta_ref_dz(k) = min(d_theta_ref_dz(k),-1e-3)
      else
          d_theta_ref_dz(k) = max(d_theta_ref_dz(k),1e-3)
      endif

       
      theta_W(k) = (theta_avg(k) - theta_ref(k))/ d_theta_ref_dz(k)
       
    !  write(mpp_pe()+2000,*) k, ' ',d_theta_ref_dz(k)
      write(mpp_pe()+2000,*) k, ' ',theta_ref(k)

!Calculate Brunt-Vaisala frequency (N):
! N^2 = g * (d_theta_ref/dz)/theta_ref
!do k=top, km

      NN(k) = sqrt( abs( d_theta_ref_dz(k) ) * g /theta_avg(k) )

!      write(mpp_pe()+1000,*) k, ' ',theta_W(k)
!      write(mpp_pe()+2000,*) k, ' ',NN(k)
enddo

    

! vertical integral average of NN on nonuniform vertical grid

N=0
do k=top,km
   if (k==top) then
      int_N(k) = NN(k) * (z_avg(k)-z_avg(k+1))
   elseif (k==km) then
      int_N(k) = NN(k) * (z_avg(k-1)-z_avg(k))

   elseif (k>top.and.k<km) then
      int_N(k) = 0.5*( NN(k+1)+NN(k))*(z_avg(k)-z_avg(k+1))
   endif

  ! write(mpp_pe()+8000,*) k, ' ', int_N(k)
   N=N+int_N(k)
enddo

N=N/H_tropopause

!write(mpp_pe()+9000,*) N

! spectral theta:

do j=1,JJ  ! J is input = 10

   theta_mode(j)=0
   m(j) = j*pi/H_tropopause

   ! vertical integral 
   do k=top, km
      if (k==top) then
         int_theta(k) = theta_W(k)*sin(m(j)*z_avg(k)) * (z_avg(k)-z_avg(k+1))
      elseif (k==km) then
         int_theta(k) = theta_W(k)*sin(m(j)*z_avg(k)) * (z_avg(k-1)-z_avg(k))
      elseif (k>top.and.k<km) then
         int_theta(k) = 0.5*( theta_W(k+1)*sin(m(j)*z_avg(k+1)) + theta_W(k)*sin(m(j)*z_avg(k)) )*(z_avg(k)-z_avg(k+1))
      endif

      theta_mode(j) = theta_mode(j) + int_theta(k)
   enddo

   theta_mode(j) = theta_mode(j)*2/H_tropopause

! tau_mode

   tau_mode(j) = j*pi*L/H_tropopause/N   ! L is namelist input ~ 200,000 m

   write(mpp_pe()+6000,*) j, ' ',m(j)
   write(mpp_pe()+7000,*) j, ' ',theta_mode(j)
   write(mpp_pe()+8000,*) j, ' ',tau_mode(j)

enddo



! Now calculate W_Largescale:
!W_Large=0
do k=top, km
   do j=1,JJ

      W_Largescale(k) = W_Largescale(k) + theta_mode(j)*sin(m(j)*z_avg(k))/tau_mode(j)
   enddo
   write(mpp_pe()+9000,*) k, ' ',W_Largescale(k)

enddo


endif
!-------------!
! end of SWTG !
!-------------!

!endif ! end W_Largescale mode
!===========================================================================!
! uanber: vertical advection of temperature and moisture by W_Largescale:   !
! W dtheta/dz and W dq/dz. Since dtheta/dz will be on half levels, W should !
! also be on half levels, that's why *0.5.                                  !
!===========================================================================!

!===================================!
! initialize advection terms        !
!===================================!

do k=1, km
   W_dT_dz(k)=0
   W_dq_dz(k)=0
enddo

! calculate advection terms:

do k=top, km

   if (k>top.and.k<km) then
      if (W_Largescale (k)<0) then ! if W_Largescale negative use upwind scheme, if W<0 (descending),then use points at k+1 and k. If W>0 (ascending) use points k, and k-1
       W_dT_dz(k)     = 0.5* (W_Largescale(k-1) + W_Largescale(k))*(theta_avg(k-1) - theta_avg(k))/(z_avg(k-1) - z_avg(k))
       W_dq_dz(k)     = 0.5* (W_Largescale(k-1) + W_Largescale(k))*(q_avg(k-1) - q_avg(k))/(z_avg(k-1) - z_avg(k))
      else
       W_dT_dz(k)     = 0.5* (W_Largescale(k) + W_Largescale(k+1))*(theta_avg(k) - theta_avg(k+1))/(z_avg(k) - z_avg(k+1))
       W_dq_dz(k)     = 0.5* (W_Largescale(k) + W_Largescale(k+1))*(q_avg(k) - q_avg(k+1))/(z_avg(k) - z_avg(k+1))
      endif
    endif
    if (k==km) then  ! surface
      W_dT_dz(k)     =   W_Largescale(k) * (theta_avg(k-1) - theta_avg(k))/(z_avg(k-1) - z_avg(k))
      W_dq_dz(k)     =   W_Largescale(k) * (q_avg(k-1) - q_avg(k))/(z_avg(k-1) - z_avg(k))
   elseif (k==top) then  ! top
      W_dT_dz(k)     =   W_Largescale(k) * (theta_avg(k) - theta_avg(k+1))/(z_avg(k) - z_avg(k+1))
      W_dq_dz(k)     =   W_Largescale(k) * (q_avg(k) - q_avg(k+1))/(z_avg(k) - z_avg(k+1))
   endif

! debug 
!     write(mpp_pe()+3000,*) k, ' ',W_dT_dz(k)
!     write(mpp_pe()+4000,*) k, ' ',W_dq_dz(k)
enddo

!==========================================================================!
! uanber: Horizontal advection of moisture: (q_ref - q)*dW_Largescale/dz   !
!==========================================================================!

if (q_horiz_adv) then

!===============================!
!  initialize
!===============================!

do k=1, km
   dW_dz(k)=0
   q_H_dW_dz(k) =0
enddo

do k=top, km
   if (k>top .and. k<km) then
      dW_dz(k) = 0.5*( (W_Largescale(k-1)-W_Largescale(k))/(z_avg(k-1)-z_avg(k)) &
                       +(W_Largescale(k)-W_Largescale(k+1))/(z_avg(k)-z_avg(k+1)) )
   endif
   if (k==top) then  ! top
      dW_dz(k) = (W_Largescale(k)-W_Largescale(k+1))/(z_avg(k)-z_avg(k+1))
   elseif (k==km) then ! surface
      dW_dz(k) = (W_Largescale(k-1)-W_Largescale(k))/(z_avg(k-1)-z_avg(k))
   endif

   if (W_Largescale(k)>0) then !bring moisture to the domain in the case of upward large-scale vertical motion
      q_H_dW_dz =  (moisture_ref(k) - q_avg(k))*dW_dz(k)
   else ! if W_Largescale <0
      q_H_dW_dz = 0.0
   endif
enddo

else
   q_H_dW_dz = 0.0

endif ! endif for horizontal moisture advection


!======================!
!    debug             !
!======================!

!!!do k=1, km
!   write(mpp_pe()+2000,*) k, ' ',theta_ref(k)
!   write(mpp_pe()+3000,*) k, ' ',theta_avg(k)
!   write(mpp_pe()+4000,*) k, ' ', d_theta_ref_dz(k)
!   write(mpp_pe()+5000,*) k, ' ',W_Largescale(k)
!   write(mpp_pe()+6000,*) k, ' ',W_dT_dz(k)
!   write(mpp_pe()+5000,*) k, ' ',W_dq_dz(k)
! enddo


do k=1, km
   do j=js,je
      do i=is,ie
         W_Largescale3D(i,j,k) = 0.0
      enddo
   enddo
enddo


!do k=1, km
!   write(mpp_pe()+1000,*) k, ' ',  t_dt(30,30,k)
!enddo

!do k=1, km
!   write(mpp_pe()+2000,*) k, ' ', W_dT_dz(k)
!enddo

!==============!
! Tendencies   !
!==============!
   do k=1, km
      do j=js,je
         do i=is,ie
            t_dt(i,j,k) = t_dt(i,j,k) - W_dT_dz(k)     ! keep in mind that this is the potential temperature tendency since we're using Z height!
            q_dt(i,j,k,1) = q_dt(i,j,k,1) - W_dq_dz(k) + q_H_dW_dz(k)
         enddo
      enddo
   enddo

!=========================================!
! update temperature iand moisture fields !
!=========================================!
! not necessary

 do k=1, km
    do j=js, je
       do i=is, ie
          temp(i,j,k)  = temp(i,j,k) + pdt*t_dt(i,j,k)
          q3 (i,j,k,1) = q3(i,j,k,1) + pdt*q_dt(i,j,k,1)
       enddo
    enddo
 enddo

!============================================! 
! note that zfull (so are z_avg, dz, delz, and dz_avg) is the geopotential height(g*Z = phi).
! multiply by g to convert to geometric height
do k=1, km
   do j=js,je
      do i=is,ie
         W_Largescale3D(i,j,k) = 9.81*100*W_Largescale(k)  ! convert to cm/s
         !write(mpp_pe()+9000,*) k, ' ', W_Largescale3D(32,32,k)
      enddo
   enddo
enddo

!do k=1, km
!   write(mpp_pe()+5000,*) k, ' ', t_dt(30,30,k)
!enddo

!do k=1, km
!   write(mpp_pe()+9000,*) k,' ', W_Largescale(k)
!   write(mpp_pe()+2000,*) k,' ', theta_avg(k)
!   write(mpp_pe()+3000,*) k,' ', theta_ref(k) 
!enddo


!endif  ! end of W_Largescal_module
!end subroutine W_Largescale_Driver

!------------------------------------------------------
!  uanber: output sensible and moisture fluxes, W_Largescale, ...ect.
!------------------------------------------------------


! output W_Largescale

if (id_W_Largescale > 0) then
    used=send_data(id_W_Largescale, W_Largescale3D(is:ie,js:je,:), Time)
endif


! debuging
!do k=1,km
!   write(mpp_pe()+2000,*)  W_Largescale(k)
!enddo


! height in m
if (id_zfull > 0) then
   used=send_data(id_zfull, zfull(is:ie,js:je,:), Time)
endif

! pressure in Pa
if (id_pfull > 0) then
   used=send_data(id_pfull, pfull(is:ie,js:je,:), Time)
endif


! potential temp in K
if (id_theta > 0) then
   used=send_data(id_theta, theta(is:ie,js:je,:), Time)
endif

end subroutine W_Largescale_Driver

!------------------------------------------------------
! domain mean wind relaxation to zero or other value
!-----------------------------------------------------

!if ( zero_winds ) then

! Zero out (doubly periodic) domain averaged winds with time-scale tau_zero:
! This loop can not be openMP-ed
! do k=1, km
!    utmp(k) = g0_sum(u3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
!    vtmp(k) = g0_sum(v3(is,js,k), is, ie, js, je, 0, gridstruct%area(is:ie,js:je), 1)
!#ifdef PRINT_W
!    if ( master ) then
!        if( sqrt(utmp(k)**2+vtmp(k)**2) > 1. ) then
!            write(*,*) k, 'Domain avg winds=', utmp, vtmp
!        endif
!    endif
!#endif
!  enddo

!
!uanber: debug global avg.

!do k=1, npz
!   write(mpp_pe()+20,*) k, ' ', utmp(k)
!enddo


!
!  rate_w = min(1., 1./(tau_zero*86400.))
!! $omp parallel do default(shared)
! do k=1, km
!    do j=js,je
!       do i=is,ie
!          u_dt(i,j,k) = u_dt(i,j,k) - utmp(k) * rate_w
!          v_dt(i,j,k) = v_dt(i,j,k) - vtmp(k) * rate_w
!       enddo
!    enddo
! enddo

!endif

!end subroutine W_Largescale_Driver


!=========================================!
!   domain average function def.          !
!=========================================!

! real function g0_sum(p, ifirst, ilast, jfirst, jlast, mode, reproduce, isd, ied, jsd, jed, area)

 real function g0_sum(p, ifirst, ilast, jfirst, jlast, mode, reproduce, area)
 use fv_mp_mod,         only: mp_reduce_sum 
! Fast version of global sum; reproduced if result rounded to 4-byte
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast
      integer, intent(IN) :: mode  ! if ==1 divided by global area
      logical, intent(IN) :: reproduce
      real,    intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
    !  integer, intent(IN) :: isd, ied, jsd, jed
!       real,    intent(IN) :: area(isd:ied,jsd:jed)
      real,    intent(IN) :: area(ifirst:ilast,jfirst:jlast)
      integer :: i,j
      real :: gsum, global_area

     global_area =0.
      do j=jfirst,jlast
         do i=ifirst,ilast
           ! write(*,*) area (i,j)
            global_area = global_area + area(i,j)
         enddo
      enddo
      call mp_reduce_sum(global_area)

      gsum = 0.
      do j=jfirst,jlast
         do i=ifirst,ilast
            gsum = gsum + p(i,j)*area(i,j)
         enddo
      enddo
!     call mpp_sum(gsum)    ! does this work?
      call mp_reduce_sum(gsum)

      if ( mode==1 ) then
           !g0_sum = gsum / (4.*pi*radius**2)
           g0_sum = gsum / global_area
      else
           g0_sum = gsum
      endif

      if ( reproduce ) g0_sum = real(g0_sum, 4) ! convert to 4-byte real

 end function g0_sum

!=========================================!
!   domain average function def.          !
!=========================================!
!real function g0_sum_old(p, ifirst, ilast, jfirst, jlast, ngc, area, mode)
! use mpp_mod,           only: mpp_sum
!      real, save :: global_area

! Fast version of globalsum
!      integer, intent(IN) :: ifirst, ilast
!      integer, intent(IN) :: jfirst, jlast, ngc
!      integer, intent(IN) :: mode  ! if ==1 divided by area
!      real, intent(IN) :: p(ifirst:ilast,jfirst:jlast)      ! field to be summed
!      real, intent(IN) :: area(ifirst-ngc:ilast+ngc,jfirst-ngc:jlast+ngc)
!      integer :: i,j
!      real gsum

!-------------------------
! Quick local sum algorithm
!-------------------------
!      if ( .not. g0_sum_initialized ) then
!         allocate (l_area(ifirst:ilast,jfirst:jlast))
!         global_area = 0.
!         do j=jfirst,jlast
!           do i=ifirst,ilast
!             global_area = global_area + area(i,j)
!             l_area(i,j) = area(i,j)
!           enddo
!         enddo
!         call mpp_sum(global_area)
!!        if ( mpp_pe().eq.mpp_root_pe() ) write(*,*) 'Global Area=',global_area
!         g0_sum_initialized = .true.
!      end if

!      gsum = 0.
!      do j=jfirst,jlast
!        do i=ifirst,ilast
!          gsum = gsum + p(i,j)*l_area(i,j)
!        enddo
!      enddo
!      call mpp_sum(gsum)
!
!      if ( mode==1 ) then
!        g0_sum = gsum / global_area
!      else
!        g0_sum = gsum
!      endif
!
! Make it reproducible by truncating to 32-bit precision:
!      g0_sum = real(g0_sum, 4)

! end function g0_sum_old

!  subroutine W_Largescale_init(is, ie, js, je, nwat, ts, time, axes, lat)

! subroutine W_Largescale_init(time, axes)
 ! integer, intent(IN) :: is, ie, js, je
 ! integer, intent(IN) :: nwat
 ! real, INTENT(inout)::    ts(is:ie,js:je)
 ! real, INTENT(IN) :: lat(is:ie,js:je)
!  integer,         intent(in) :: axes(4)
!  type(time_type), intent(in) :: time
!  integer :: unit, ierr, io, i, j
!  real:: total_area

!   ----- read and write namelist -----
!    if ( file_exist('input.nml')) then
!         unit = open_namelist_file ('input.nml')
!         read  (unit, nml=W_Largescale_nml, iostat=io, end=10)
!         ierr = check_nml_error(io,'W_Largescale_nml')
! 10     call close_file (unit)
!    endif


!    id_W_Largescale = register_diag_field(mod_name, 'W_Largescale3D', axes(1:3), Time, 'W Large-scale', 'cm/s')

!    write(mpp_pe()+2000,*) id_W_Largescale

!end subroutine W_Largescale_init




 subroutine W_Largescale_init(axes,time)

 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: time

 integer :: unit, ierr, io

!   ----- read and write namelist -----
!   ----- read and write namelist -----
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=W_Largescale_nml, iostat=io)
    ierr = check_nml_error(io,'W_Largescale_nml')
#else
    if ( file_exist('input.nml')) then
       unit = open_namelist_file( )
       ierr=1; do while (ierr /= 0)
       read  (unit, nml=W_Largescale_nml, iostat=io, end=10)
       ierr = check_nml_error(io,'W_Largescale_nml')
       enddo
  10   call close_file (unit)
    endif
#endif
if (mpp_pe() == mpp_root_pe() ) then
    write(6,*) 'TEST OUTPUT ', tau_WTG, wavenumber, damping, W_Largescale_flag
endif


!    if ( file_exist('input.nml')) then
!         unit = open_namelist_file ('input.nml')
!         read  (unit, nml=W_Largescale_nml, iostat=io, end=10)
!         ierr = check_nml_error(io,'W_Largescale_nml')
! 10     call close_file (unit)
!    endif


!uanber: W_Largescale
    id_W_Largescale = register_diag_field(mod_name, 'W_Largescale3D', axes(1:3), Time, 'W_Largescale3D', 'cm/s')

 end subroutine W_Largescale_init


end module W_Largescale_mod
