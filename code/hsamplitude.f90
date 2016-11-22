MODULE hsamplitude
  
contains
  real(kind(1.d0)) function Amplitude_triplet(kF,a,m,hbar,theta,phi)
    implicit none
    integer, parameter:: dp = kind(1.d0)
    real(kind=dp):: kF, a, theta, phi
    real(kind=dp):: m, hbar, pi
    real(kind=dp):: const, at
    
    pi = acos(-1.d0)
    const = 4.d0 * hbar**2*kF*a**2/m
    At = const *( Ufunc(theta,phi) - Ufunc( theta, pi + phi) )
    Amplitude_triplet = At
    return
  end function Amplitude_Triplet
  




  real(kind(1.d0)) function Amplitude_singlet(kF,a,m,hbar,theta,phi)
    implicit none
    integer, parameter:: dp = kind(1.d0)
    real(kind=dp):: kF, a, theta, phi
    real(kind=dp):: m, hbar, pi
    real(kind=dp):: const, aa, as
    
    pi = acos(-1.d0)

    const = 16.d0 * pi * hbar**2*a / m
    aa  =  3.d0 - sin(theta/2.d0) * log( abs( (1.d0 + sin(theta/2.d0) ) / ( 1.d0 - sin(theta/2.d0) ) ) )
    aa = aa + Ufunc(theta,phi)/4.d0 + Ufunc( theta, pi + phi)/4.d0 
    As = const * (1.d0 + kF*a/pi * aa )

    Amplitude_singlet = As

    return
  end function Amplitude_Singlet
  
  
  
  real(kind(1.d0)) function Ufunc(theta, phi)
    implicit none
    integer, parameter:: dp = kind(1.d0)
    real(kind=dp):: theta, phi, theta2, phi2
    real(kind=dp):: u1, u2, uu
    
    uu = sin(theta/2.d0)* cos(phi/2.d0)    
    u1 = (1.d0 - uu**2 )/ uu
    u2 = log ( abs ( (1.d0 + uu)/ (1.d0 - uu)) )
    
    Ufunc = u1*u2
    return
  end function Ufunc
  
end MODULE hsamplitude
