function chere,z,omegam,h,Mh,noscatc=noscatc


 Ms=Mh+alog10(h) ;to convert to Msun/h

 ;read-off c-Mh from Benedikt's files (taken from http://www.benediktdiemer.com/data/)
 diro='\\soton.ac.uk\ude\PersonalFiles\Users\fs1y12\mydocuments\WORKING_DIRECTORIES\Articolo\SAMsSizeEnviron\draftrho\'
 ;#   z     nu     M200c  c200c     M500c  c500c    M2500c c2500c      Mvir   cvir     M200m  c200m
 ;readcol,diro+'Concentration_bolshoi_median.txt',zb  ,   nu  ,   M200c,  c200c  ,   M500c,  c500c   , M2500c ,c2500c  ,    Mvir ,  cvir ,    M200m , c200m,/silent
 ;readcol,diro+'Concentration_francesco_median.txt',zb  ,   nu  ,   M200c,  c200c  ,   M500c,  c500c   , M2500c ,c2500c  ,    Mvir ,  cvir ,    M200m , c200m,/silent
 ;save,filename=diro+'BenediktSig0.8Omegam0.3h0.7.sav',zb,M200c,c200c,Mvir,cvir,M200m,c200m,M500c,c500c
 ;stop
 restore,filename=diro+'BenediktSig0.8Omegam0.3h0.7.sav';,zb,M200c,c200c,Mvir,cvir,M200m,c200m,M500c,c500c

 jc=value_locate(zb,z) ; index where z sits in this sequence
 zc=zb(jc) ; closest indexed value? It's actually the lowest closest
 ix=where(zb eq zc) ; vector mask for redshifts of this value

 cb=c200c(ix) ; the value of c200c at the redshift in question
 Mhb=alog10(M200c(ix)) ; halo masses at this value

 logc=alog10(interpol(cb,Mhb,Ms)) ; interpolate input halo mass to cb

 ;from Chae+14:
 ;logc=alog10((5.61/(1.d0+z))*(10.^Ms/1e14)^(-0.13))

 if not keyword_set(noscatc) then scatc=0.16 else scatc=0.
 c=10.^(logc+scatc*randomn(seed,n_elements(Mh)))

 if n_elements(Mh) eq 1 then c=total(c)

 return,c

end


function rvhere,Mvir,z,h,omegam

 OmegaL=1.-Omegam
 Omegaz=Omegam*(1.+z)^3./(Omegam*(1.+z)^3.+OmegaL)
 d=Omegaz-1.
 Deltac=18.*(!pi)^2.+82.*d-39.*d^2.

 ;A1z=31.*(Omegam/Omegaz*Deltac/18./!pi^2.)^(-1./3.)*((1.+z)/7.)^(-1.);*h^(-2./3.)
 ;r=(h/0.7)^(-1.)*A1z*10.^(1./3.*(Mvir-12.)) ;kpc

 H0=100.*h   ;km/s/Mpc
 HH=H0*sqrt(Omegam*(1.+z)^3.+(1.-Omegam))
 G=4.302e-9  ;Mpc/Msun(km/s)^2
 rhoc=3.*HH^2./(8.*!pi*G)
 k=4.*!pi/3.

 ;rhom[z]=rhoc[0]*Omegam*(1+z)^3=(rhoc[0]*E[z])*(Omegam*(1+z)^3/E[z])=rhoc[z]*omegam[z]

 ;r=(10.^Mvir/rhoc/k/Deltac)^(1./3.)*1000. ;kpc

 ;NOTE: the latter definition of halo mass is equivalent to say that the overdensity is Deltavir=(18.*(!pi)^2.+82.*d-39.*d^2.)/(1+d)=(18.*(!pi)^2.+82.*d-39.*d^2.)/Omegaz times the
 ;background density of the Universe at z, i.e., rhi_c(0)*Omegam*(1+z)^3
 ;in fact Mvir=4/3*pi*Deltac/Omegaz*Rv^3*rhoc_c*E(z)*Omegaz=4/3*pi*Deltavir*Rv^3*(rhoc_c*Omegam*(1+z)^3) where
 ;rho_background(z)=rhoc_c*Omegam*(1+z)^3

 ;virial radius defined for halos 200 times the critical density of the Universe at z:
 r=(10.^Mvir/rhoc/k/200.)^(1./3.)*1000. ;kpc

 ;from Chae+14 (about identical to the computation of Rvir when adopting the critical density of the Universe!):
 ;r=199.9*(10.^Mvir/10.^12.)^(1./3.)   ;kpc

 return,r

end


function Vvirhere,Mh,zs,h,omegam

;following Loeb & Peebles 2003 (Mh must be in units of Msun):
Omegaz=Omegam*(1.+zs)^3./(Omegam*(1.+zs)^3.+1.-Omegam)
d=Omegaz-1.
Deltac=18.*(!pi)^2.+82.*d-39.*d^2.
A1z=375.*(Omegam/Omegaz*(Deltac/(18.*(!pi)^2.)))^(1./6.)*((1.+zs)/7.)^0.5
Vh=alog10(A1z)+(1./3.)*(Mh-12.)

;following Barkana & Loeb 2001:
;Vh=alog10(23.4)+1./3*(Mh-8.)+alog10((Omegam/Omegaz*(Deltac/(18.*(!pi)^2.)))^(1./6.)*((1.+zs)/10.)^0.5)

 return,Vh ;km/s

end

