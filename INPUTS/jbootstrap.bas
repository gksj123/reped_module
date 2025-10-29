function fluxavh
# This funtion is to evaluate the <h>,<h*h>,<hf>
real s_pin = sign(1.0, dpsi0)
real hsrf(msrf), h2srf(msrf), hfsrf(msrf)
integer l
do l=1, msrf
  dq2 = s_pin * 2.*pi/(mls-1) * zlsd(:,l,4)
	dq2(mls) = 0.
  vprime = sum(dq2)
  
	real ff   = frsrf(l)    # F
	real rr   = rlsd(:,l,1) # R array on flux surface
	real bphi = ff/rr
	real bpol = rlsd(:,l,4)
	real bb   = sqrt(bpol**2 + bphi**2) # |B|
  real bmax = max(bb)
	real hh   = bb/bmax
	hsrf(l)  = sum(dq2 * hh)/vprime     # <h>
	h2srf(l) = sum(dq2 * hh**2 )/vprime # <h*h>
	real temp     = (1.0 - sqrt(1.-hh)* ( 1 + 0.5*hh))/(hh*hh)
	hfsrf(l) = sum(dq2 * temp)/vprime
enddo
return [hsrf, h2srf, hfsrf]
endf

function ftrap
real temp = fluxavh
real hsrf  = temp(:,1)
real h2srf = temp(:,2)
real hfsrf = temp(:,3)

real ftu = 1.0 - h2srf / (hsrf*hsrf) * (1. - sqrt(1. - hsrf)*(1. + 0.5*hsrf))
real ftl = 1.0 - h2srf * hfsrf
real ft  = 0.75*ftu + 0.25*ftl
return ft
endf




function jbootstrap(pres, ne, te, tii, zeff, pbeam, zion, zimp, gradb)

real rr   = rsrf/100.  # <R> in m

real nii = ( ne/zeff - ne * zimp**2) /(zion**2 - zimp**2)
real pth = pres - pbeam
real pe  = ne * te * 1.602

real epsilon = asrf/rsrf   # inverse aspace ratio=a/R 
epsilon(1) = epsilon(2)    # to avoid epsilon = 0, for later using

real ft0 = sqrt(epsilon)
real ft1 = 1.0 - (1.0-epsilon)**2 / sqrt(1.0-epsilon**2) /(1.0+1.46*sqrt(epsilon))
real ft = ftrap()          # call ftrap

#
real lambdae = 31.3 - log( sqrt(ne)/te )
real lambdai = 30.  - log( zeff**3 * sqrt(nii) / tii**1.5)

real nz = 0.58 + 0.74/(0.76 + zeff)
real sspitz = 1.9012e4 * te**1.5 / zeff / nz / lambdae

real nustare = abs(6.921e-18 * qsrf * rr * ne * zeff * lambdae / te**2 / epsilon**1.5)
#!? nustare does not need the abs operator
real nustari = abs(4.90e-18 * qsrf * rr * nii * zeff**4 * lambdai / tii**2 / epsilon**1.5)


# XGC modification for trapped particle fraction (ft)
real delta_psie = rr * sqrt( epsilon * 2.0 * te * 9.11e-31 / 1.602e-19)  #banana width
if (gradb > 0 ) then
  real h_xgc = 1.- 0.2/zeff**4 * exp( - abs( (1. - psibar)*dpsi0 /2.7 /log(epsilon**1.5 *
	             nustare /3.2 +3) /delta_psie )) 
else
  real h_xgc = 1.- 0.6/zeff**4 * exp( - abs( (1. - psibar)*dpsi0 /3.3 /log(epsilon**1.5 *
	             nustare + 2.) /delta_psie))
endif
real ft_xgc = ft * h_xgc  
real beta_xgc = float( (epsilon - (0.44+0i))**0.7)

real az = (-zeff**2 + 5.998*zeff -4.981) / (4.294* zeff**2 - 14.07*zeff + 12.61)
integer l
do l=1,msrf
  if (zeff(l) > 5.0)then
	  az(l) = 0.
	endif
enddo
real zaz = 1/ (zeff**az)
real delta_xgc = 0.55 * zeff**0.2 * ( tanh(3.2 * beta_xgc * (epsilon**1.5 * nustare)**1.4 *
                  zaz) + (1. - exp(-nustare/0.1)) * tanh(2.2 * beta_xgc * epsilon**2.8 *
									nustare**0.1 * zaz))

#
real sqnue = sqrt(nustare)
real f31     = ft / ( 1. +  ( 1.00 - 0.10*ft ) * sqnue + 0.50*( 1.00 - ft ) * nustare/zeff )
real f31_xgc = ft_xgc *( 1. + delta_xgc ) / (
        1. + ( 1.00 - 0.10*ft_xgc ) * sqnue + 0.50*( 1.00 - ft_xgc ) * nustare/zeff )

real f32_ee = ft / (
    1. + 0.26*( 1.00 -      ft ) * sqnue + 0.18*( 1.00 - 0.37*ft ) * nustare/zeff**0.5 )
real f32_xgc_ee = ft_xgc*( 1. + delta_xgc ) / (
    1. + 0.26*( 1.00 -      ft_xgc ) * sqnue + 0.18*( 1.00 - 0.37*ft_xgc ) * nustare/zeff**0.5 )

real f32_ei = ft / (
    1. +      ( 1.00 + 0.60*ft ) * sqnue + 0.85*( 1.00 - 0.37*ft ) * nustare*(1.+zeff) )
real f32_xgc_ei = ft_xgc*( 1. + delta_xgc ) / (
    1. +      ( 1.00 + 0.60*ft_xgc ) * sqnue + 0.85*( 1.00 - 0.37*ft_xgc ) * nustare*(1.+zeff) )

# equation 13b, article 1a
real f33    = ft / (
    1. +      ( 0.55 - 0.10*ft ) * sqnue + 0.45*( 1.00 -      ft ) * nustare/zeff**1.5 )

# equation 16b, article 1a
real f34    = ft /(
    1. +      ( 1.00 - 0.10*ft ) * sqnue + 0.50*( 1.00 - 0.50*ft ) * nustare/zeff      )
# xgc0 modification, article 3
real f34_xgc = ft_xgc*( 1. + delta_xgc ) /(
    1. +      ( 1.00 - 0.10*ft_xgc ) * sqnue + 0.50*( 1.00 - 0.50*ft_xgc ) * nustare/zeff      )

# equation 13a, article 1a
real z1 = 1. / zeff
real f33 = ( ( - 0.23*z1 *f33 + 0.59*z1 )   \
                           *f33 - (1.+0.36*z1) ) \
                           *f33 + 1.0

# neoclassical conductivity
real sneo = f33*sspitz

# equation 14a, article 1a
real z1 = 1./( 1. + zeff )
real l31 = ( ( ( 0.2*z1 *f31 + 0.3*z1 )  \
                          *f31 - 1.9*z1 )  \
                          *f31 + (1. + 1.4*z1) ) \
                          *f31
# xgc0 modification, article 3
real l31_xgc = ( ( ( 0.2*z1 *f31_xgc + 0.3*z1 )  \
                              *f31_xgc - 1.9*z1 )  \
                              *f31_xgc + (1. + 1.4*z1) ) \
                              *f31_xgc

# equation 16a, article 1a
real l34 = ( ( ( 0.2*z1 *f34 + 0.3*z1 )  \
                          *f34 - 1.9*z1 )   \
                          *f34 + (1. + 1.4*z1) ) \
                          *f34
# xgc0 modification, article 3
real l34_xgc = ( ( ( 0.2*z1 *f34_xgc + 0.3*z1 )   \
                              *f34_xgc - 1.9*z1 )  \
                              *f34_xgc + (1. + 1.4*z1) ) \
                              *f34_xgc


# equation 15b, article 1a
real f32_ee = (
    ( 0.05 + 0.62*zeff ) /( zeff * ( 1.+ 0.44*zeff ) ) *
    ( f32_ee - f32_ee**4 )  \
    + 1.00/( 1.+ 0.22*zeff ) *
    ( f32_ee**2 - f32_ee**4 - 1.20*( f32_ee**3 - f32_ee**4 ) ) \
    + 1.2/( 1. + 0.5*zeff ) *
    f32_ee**4  \
    )
# xgc0 modification, article 3
real f32_xgc_ee = (
    ( 0.05 + 0.62*zeff ) /( zeff * ( 1.+ 0.44*zeff ) ) *
    ( f32_xgc_ee - f32_xgc_ee**4 )  \
    + 1.00/( 1.+ 0.22*zeff ) *
    ( f32_xgc_ee**2 - f32_xgc_ee**4 - 1.20*( f32_xgc_ee**3 - f32_xgc_ee**4 ) ) \
    + 1.2/( 1. + 0.5*zeff ) *
    f32_xgc_ee**4 \
    )

# equation 15c, article 1a
real f32_ei = (
    -( 0.56 + 1.93*zeff ) /( zeff * ( 1.+ 0.44*zeff ) ) *
    ( f32_ei - f32_ei**4 ) \
    + 4.95/( 1.+ 2.48*zeff ) *
    ( f32_ei**2 - f32_ei**4 - 0.55*( f32_ei**3 - f32_ei**4 ) ) \
    - 1.2/( 1. + 0.5*zeff ) *
    f32_ei**4 \
    )
# xgc0 modification, article 3
real f32_xgc_ei = (
    -( 0.56 + 1.93*zeff ) /( zeff * ( 1.+ 0.44*zeff ) ) *
    ( f32_xgc_ei - f32_xgc_ei**4 ) \
    + 4.95/( 1.+ 2.48*zeff ) *
    ( f32_xgc_ei**2 - f32_xgc_ei**4 - 0.55*( f32_xgc_ei**3 - f32_xgc_ei**4 ) ) \
    - 1.2/( 1. + 0.5*zeff ) *
    f32_xgc_ei**4 \
    )

# equation 15a, article 1a
real l32 = f32_ei + f32_ee
# xgc0 modification, article 3
real l32_xgc = f32_xgc_ei + f32_xgc_ee

# equation 17a, article 1a
real alf0 = -1.17 * ( 1. - ft )/( 1. - 0.22*ft - 0.19*ft*ft )

# equation 1, article 1b
real sqnui   = sqrt( nustari )
real st = nustari * nustari * ft**6
real alf = (
              ( alf0 + 0.25*(1.-ft**2)*sqnui )/(1.+0.5*sqnui) + 0.315*st  \
             ) / ( 1. + 0.15*st )

real dpth(msrf), dte(msrf), dti(msrf)
call interp(psibar, pth, msrf, -111, -111, psibar, pth, &dpth, msrf)
call interp(psibar, te , msrf, -111, -111, psibar, te , &dte , msrf)
call interp(psibar, tii, msrf, -111, -111, psibar, tii, &dti , msrf)
real dpsi0_sign = sign(dpsi0, jtsrf(1))  # the sign of dpsi0 should be consistent with j
dpth = dpth / dpsi0_sign * 1.e8     # convert dpsi0 to T*m^2/radian, but why /radian?
dte  = dte  / dpsi0_sign * 1.e8
dti  = dti  / dpsi0_sign * 1.e8
#?! why dp/dpsi, not d(lnp)/dpsi as in article? And why psi is in /radian?

# Equation 2, article 1b
real rpe = pe/pth
real gpp     = l31    *dpth/pe
real gp_xgc  = l31_xgc*dpth/pe
real gte     = l32    *dte /te
real gte_xgc = l32_xgc*dte /te
real gti     = l34    *dti /tii*alf*(1.-rpe)/rpe
real gti_xgc = l34_xgc*dti /tii*alf*(1.-rpe)/rpe

# Currents are in A/m**2
# This form is <Jparallel*B>/(f/R0)
real jp      = -1.e-2*rcntr*pe*gpp 
real jp_xgc  = -1.e-2*rcntr*pe*gp_xgc 
real jte     = -1.e-2*rcntr*pe*gte    
real jte_xgc = -1.e-2*rcntr*pe*gte_xgc
real jti     = -1.e-2*rcntr*pe*gti    
real jti_xgc = -1.e-2*rcntr*pe*gti_xgc
real jbs     = jp     + jte     + jti    
real jbs_xgc = jp_xgc + jte_xgc + jti_xgc

real jbs_teqpar     = jbs * 1.e-4 * frsrf(1)/rcntr /bsqrf   # convert to corsica parallel current
real jbs_teqpar_xgc = jbs_xgc * 1.e-4 * frsrf(1)/rcntr /bsqrf
real jbs_teqt       = jbs_teqpar * bsqrf /frsrf /rsqirf     # convert to corsica ohmic current
real jbs_teqt_xgc   = jbs_teqpar_xgc * bsqrf /frsrf /rsqirf

return [jbs, jbs_xgc, nustare, nustari, sspitz, sneo]
endf

