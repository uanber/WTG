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

 logical, INTENT(IN):: hydrostatic

 type(time_type), intent(in) :: Time


 real, INTENT(INOUT)::   pt(is:ie,js:je,npz)!, theta(is-ng:ie+ng,js-ng:je+ng,npz) !theta is potential temperature 
 !real, INTENT(INOUT):: delp(is:ie,js:je,npz)
 real, INTENT(INOUT)::    q(is:ie,js:je,npz,nq)  !
 !real, INTENT(INOUT):: delz(is:ie,js:je,npz)
 real,            intent(in),    dimension(is:ie,js:je,npz) :: pfull, zfull
 real, intent(in), dimension(is:ie,js:je,npz+1) :: phalf, zhalf
 real, intent(in), dimension(:,:)   :: area


! Tendencies:
 real, INTENT(INOUT):: t_dt(is:ie,js:je,npz)
 real, INTENT(INOUT):: q_dt(is:ie,js:je,npz,nq)

!================================================!
!uanber: W_Largescale. define variables as local !
!================================================!
 real, dimension(npz) :: W_Largescale, W_dT_dz, W_dq_dz, q_H_dW_dz
 real, dimension(npz) :: d_theta_ref_dz
 real, dimension(is:ie,js:je,npz):: theta, temp
 real, dimension(is:ie,js:je,npz,nq):: q3
 real, dimension(npz):: dW_dz !vertical advection of T, q, dW/dz,horizontal advection of q (respect.)
 real, INTENT(IN), dimension(is:ie, js:je) :: pbltop  ! top of the PBL   
 real :: H_tropopause ! tropopause height in m
 integer:: K_PBL!, W_Largescale_flag  !PBL_K is the first model level above the boundary layer.

 real, parameter:: g= 9.81   !gravity acceleration
! wevecoupling
 real, dimension(is:ie,js:je,npz):: T_v!(is-ng:ie+ng,js-ng:je+ng,npz)
 real, dimension(is:ie,js:je,npz):: T_v_ref!(is-ng:ie+ng,js-ng:je+ng,npz)
 real, dimension(npz):: Tv_ref, Tv_avg

 real, dimension(npz)::a, b, c, e, f, y,RHS, D

 real const

!=================================================!
!     for spectral WTG                            !
!=================================================!
real, dimension(npz) :: theta_W, NN, int_N, int_theta ! integral N and integral theta
real :: N
integer , parameter:: JJ=10  ! number of vertical modes
!J=10
real, dimension(JJ) :: theta_mode, m, tau_mode

!=============================================================!

 real, dimension(is:ie,js:je,npz) :: p3, dz
 real, dimension(is:ie,js:je):: ps !, qs, rho, clouds
! real, dimension(is:ie,js:je):: olr, lwu, lwd, sw_surf, wet_t, net_rad
 real, INTENT(OUT), dimension(is:ie,js:je,npz):: W_Largescale3D
! Flux diag:
! real, dimension(is:ie,js:je):: rflux, qflux
 real, dimension(is:ie,npz):: den   ! density
 real :: tvm, rrg, zvir
 integer  i,j,k, km, iq, k_mp, top, srf ! uanber: top is the height of the tropopause. srf =km-1
 integer  seconds, days
 logical print_diag, used

!============================!
! global average variables   !
!============================!

real, dimension(npz):: z_avg, theta_avg, T_avg, q_avg, p_avg, dz_avg, dp_avg, tdt_avg, theta_diff
real :: pbltop_avg



   zvir = rvgas/rdgas - 1.
   rrg  = rdgas / grav
   km = npz


   do k=km,1,-1  !Linjoing: reversed the loop as in hydrostatic case above
     do j=js,je
        do i=is,ie
                dz(i,j,k) = zhalf(i,j,k) - zhalf(i,j,k+1)               
         enddo
      enddo
   enddo

!---------------------------------------------!
!      calculate potential temperature        !
!---------------------------------------------!

do k=1, km
   do j=js,je
      do i=is,ie
          theta(i,j,k)=pt(i,j,k)*( 1.0E5/pfull(i,j,k) )**(kappa)  !potential temperature (phalf(i,ikm+1,j) is surface pressure)
         temp(i,j,k)=pt(i,j,k)  ! temperature
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! uanber: domain mean averages for zfull, potential temperature theta, temperature pt, and moisture q: !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1, km

   z_avg(k)     = g0_sum(zfull(is,js,k), is, ie, js, je, 1,.true., area)
   theta_avg(k)     = g0_sum(theta(is,js,k), is, ie, js, je, 1,.true., area) 
   theta_diff(k) = theta_avg(k) - theta_ref(k)
   T_avg(k)     = g0_sum(temp(is,js,k), is, ie, js, je, 1,.true., area)
   q_avg(k)     = g0_sum(q(is,js,k,1), is, ie, js, je, 1,.true., area)
   dz_avg(k)     = g0_sum(dz(is,js,k), is, ie, js, je, 1,.true., area)
  
enddo

   pbltop_avg    = g0_sum(pbltop(is,js), is, ie, js, je, 1,.true., area)



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


!=========================!
! initialize W_Largescale !
!=========================!

do k=1, km
   W_Largescale(k)=0
enddo


! WTG flag
!========================================!
!  Weak Temperature Gradient (WTG):      !    
!========================================!


if (W_Largescale_flag == 1) then    ! flag for: 1- WTG


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


  else ! in the boundary layer W_largescale is interpolated from the surface (=0) to the top of hte PBL (=K_PBL)
      W_Largescale(k) =  (W_Largescale(K_PBL) - 0.0)* (z_avg(k) - z_avg(km))/(z_avg(K_PBL)-z_avg(km))
      d_theta_ref_dz(k) =0.0 ! mixing in the boundary layer so no vertical gradient of temperature.
  endif

  
 
 !   write(mpp_pe()+8000,*) k, ' ',theta_ref(k)
     write(mpp_pe()+9000,*) k, ' ',W_Largescale(k)

enddo


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




!uanber: W_Largescale
    id_W_Largescale = register_diag_field(mod_name, 'W_Largescale3D', axes(1:3), Time, 'W_Largescale3D', 'cm/s')

 end subroutine W_Largescale_init


end module W_Largescale_mod
