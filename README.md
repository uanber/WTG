# WTG

 This module is used for parameterizing the large scale vertical
 velocity (W_LargeScale) (or similarly the mass flux) in the GFDL cloud 
 resolving limited-domain model.
 Three methods can be used to do this:

 1- Weak Temperature Gradient (WTG) method: relaxing domain mean virtual
  temperature towards a prescribed profile over a  constant relaxation time
  scale (tau_wtg). This method is implemented using a simple finite difference
  scheme.

 2- Damped Gravity Waves (DGW) method: simulating a single gravity wave doing
  the actual adjustment of the density anomaly (or virtual temperature
  anomaly). Anomaly is difference from a prescribed profile as above. This
  method is implemented using a finite element scheme that solves the linear
  tridiagonal system. Note: this is advection diffusion equation with variable
  source term (RHS) and zero advection velocity. It requires two boundary
  conditions at surface and tropopause W =0.

 3- Spectral Weak Temperature Gradient (SWTG) method: spectral (Fourier)
  decomposition of the virtual temperature anomaly (from the prescribed profile).
  sine functions are chosen as bases for the decomposition as they represent
  wave adjustment in the vertical. This method assumes variable relaxation time
  with hight. 

 REMEMBER: 1 is top layer, km is the lowest layer    


# Weak Temperature Gradient (WTG): 

```
  W*dtheta_ref/dz  = (theta - theta_ref) /tau
 
  In the PBL:
  W = interpolation from above PBL to zero at the surface
```

# Damped Gravity Waves: DGW  

```
 d^2W/d^2Z = const*(T_v - T_v_ref)/T = RHS

 (1/dz) * W_k-1  - (2/dz) * W_k + (1/dz) * W_k+1 = RHS* dz
      W_k-1  -2* W_k +  W_k+1 = RHS* dz*dz =D
 ```
 
 or:
 
 ```
 B*W = D

              | -2  1  0  0 0 . . .  0 |
              |  1 -2  1  0 0 . . .  0 |
 B=           |  0  1 -2  1 0 . . .  0 |
              |  . . . . . .  . . . .  |
              |  0  0  0  0 .. 1 -2  1 |
              |  0  0  0  0 . .0  1 -2 |

 a = -2 ; b=c= 1
```
 Construction of L and U:
``` 
 L = 

      |1  0  0  0 . . . . 0|
      |e2 1  0  0 . . . . 0|
      |0  e3 1  0 . . . . 0|
      |.  .  .  .  .  . . 0|
      |0  0  0  . e_n-1 1 0| 
      |0  0  0  .  .  e_n 1|
 U =
       |f1  b1  0  0  . . .    0 |     
       |0   f2  b2  0  . . .   0 |
       |0   0   f3  b3  . .    0 |
       |.   .   .   .   . . .    |
       |0   0   0   .f_n-1  b_n-1|
       |0   0   0  . .  .    fn  |

 Now: BW=D  ==> (LU)W=D ==> L(UW)=D
      Ly = D ...(1) solve for y
      UW = y ...(2) solve for W
 ```


# Spectral WTG: SWTG         

 A- loop on vertical levels k:
 1- calculate: d_theta_dz 
 2- calculate theta_W = (theta_avg-theta_ref)/d_theta_dz
 3- Calculate Brunt-Vaisala frequency (NN)
 
 B - vertical integrated average of Brunt-Vaisala frequency N ~= 1e-2 s-1
 
 C- loop on vertical modes j (1D arrays):
 1- calculate m_j
 2- calculate spectral theta and vertically integrate
 3- calculate spectral tau

 D- loop over vertical levels k
  - loop over vertical modes j:
  - calculate W_Largescale by summing over dummy j


