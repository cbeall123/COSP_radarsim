! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, British Crown Copyright, the Met Office
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2015 - D. Swales - Modified for COSPv2.0
! July 2021 - C. Beall - Addition of mass-weighted fall speed calculations for implementation with EAMv1
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_mwfs
  USE COSP_KINDS,     ONLY: wp
  USE math_lib,	      ONLY: gamma
  use mod_cosp_error, ONLY: errorMessage

  implicit none
  
contains
  subroutine calc_mwfs(npoints,nlev,at,pfull,mr_hydrolsr,mr_hydrolss,mr_hydrolsl,mr_hydrolsi,&
		       Niclsr,Niclss,Niclsl,Niclsr,mwfs_lsliq,mwfs_lsice,mwfs_lsrain, mwfs_lssnow)
    INTEGER :: npoints,    &    ! Number of model points in the horizontal
               nlev             ! Number of model levels in columns
 
    INTEGER :: i,j,ilev,ncolprint

    REAL(WP), dimension(npoints,nlev) ::  &
         at,         &    ! Temperature (K)
         pfull,      &    ! Pressure (hPa)
         mr_hydrolsr,   &    ! Mixing ratio for large-scale stratiform rain (kg/kg)
	 mr_hydrolss,	&    ! Mixing ratio large-scale stratiform snow (kg/kg)
	 mr_hydrolsl,	&    ! Mixing ratio large-scale stratiform cloud liquid (kg/kg)
	 mr_hydrolsi,	&    ! Mixing ratio large-scale stratiform cloud ice (kg/kg)
         Niclsr,        &    ! Number concentration large-scale stratiform rain (kg-1)
	 Niclss,        &    ! Number concentration large-scale stratiform snow (kg-1)
	 Niclsl,        &    ! Number concentration large-scale stratiform liquid (kg-1)
	 Niclsi              ! Number concentration large-scale stratiform ice (kg-1)

    !OUTPUTS
    REAL(WP),intent(inout), dimension(npoints,nlev) :: &
         mwfs_lsliq,   &  ! Mass-weighted fall speed large-scale stratiform cloud liquid (m/s)
         mwfs_lsice,   &  ! Mass-weighted fall speed large-scale stratiform cloud ice (m/s)
	 mwfs_lsrain,  &  ! Mass-weighted fall speed large-scale stratiform rain (m/s)
	 mwfs_lssnow      ! Mass-weighted fall speed large-scale stratiform snow (m/s)

    !INTERNAL VARIABLES
    real(wp) :: &
         qsmall,	   & ! Threshold for mixing ratio considered cloud in EAMv1
	 rhow,		   & ! Bulk density liquid
	 shr_const_bolz,   & ! Boltzman constant
	 shr_const_avogad, & ! Avogadro's number 
	 shr_const_mwdair, & ! Molecular weight dry air (kg/kmol)
	 rair,		   & ! (J/kg/K)
	 tmelt,            & ! Melt temperature (K)
	 rhosu,		   & ! Typical air density at 850 mb
	 ar,		   & ! Fall speed air density adjustment
	 br,		   & ! 
	 rhosn,		   & !
	 pi,		   & !
	 cons5,		   & !
	 at,		   & !
	 bs,		   & !
	 ds,		   & !
	 cons8,		   & !
	 cons6,		   & !
	 rhosn,		   & ! Bulk density snow (kg/m3)
	 cs,		   & ! Snow mass-diameter relationship
	 bc,		   & !
	 ac,		   & !
	 ai,bi,di,	   & !
	 rhoi,		   & ! Bulk density ice (kg/m3)
	 cons17,cons1,ci,dcs,lammaxs,lammins,lamminr,lammaxr,lammini,lammaxi
     
     real(wp), dimension(npoints,nlev),allocatable :: &
      	 rhoa,              & ! Density (kg/m3)
	 rhof,              & !
	 arn,		    & !
	 asn,		    & !
	 acn,		    & !
	 ain,		    & !
	 lamr,		    & !
	 lams,     	    & !
	 lamc,		    & !
	 lami,		    & !
	 lammin,	    & !
	 lammax,	    & !
	 pgam 		      ! Spectral width of droplet size distribution (liquid)
	 

    
    ! #######################################################################
    ! Initialize working variables
    ! #######################################################################
    
    ! Initialize mwfs to 0
    mwfs_lsliq(1:npoints,1:nlev)=0.0     
    mwfs_lsice(1:npoints,1:nlev)=0.0     
    mwfs_lsrain(1:npoints,1:nlev)=0.0
    mwfs_lssnow(1:npoints,1:nlev)=0.0

    ! Constants for EAMv1 mwfs calculations
    qsmall = 1E-18
    rhow = 1000._wp
    shr_const_bolz = 1.38065E-23
    shr_const_avogad = 6.02214E26
    shr_const_mwdair = 28.966
    rair = shr_const_avogad*shr_const_bolz/shr_const_mwdair
    tmelt = 273.15
    rhosu = 85000./rair/tmelt
    ar = 841.99667
    br = 0.8
    rhosn = 250._wp
    pi = acos(-1._wp)
    cons5 = gamma(4+br)

    at = 11.72
    bs = 0.41
    ds = 3._wp
    cons8 = gamma(4+bs)
    cons6 = gamma(1+ds)
    cs = rhosn + pi/6

    bc = 2._wp
    ac = 3E7
    ai = 700._wp
    bi = 1._wp
    di = 3._wp
    rhoi = 500._wp
    cons17 = gamma(4+bi)
    cons1 = gamma(1+di)
    ci = rhoi*pi/6

    dcs = 400E-6
    lammaxs = 1._wp/1E-6
    lammins = 1._wp/2000E-6
    lamminr = 1._wp/500E-6
    lammaxr = 1._wp/20E-6

    lammini = 1._wp/2._wp*dcs
    lammaxi = 1._wp/1E-6
    
    allocate(rhoa(npoints,nlev),rhof(npoints,nlev),arn(npoints,nlev),asn(npoints,nlev),lamr(npoints,nlev),&
	     lams(npoints,nlev),lamc(npoints,nlev),pgam(npoints,nlev),lammin(npoints,nlev),lammax(npoints,&
	     nlev),acn(npoints,nlev),ain(npoints,nlev),lami(npoints,nlev))

    ! #######################################################################
    ! Mass-weighted fall speed calculations for large-scale stratiform hydrometeors
    ! #######################################################################

    !Calculate density
    rhoa(1:npoints,1:nlev)=0.
    rhoa(1:npoints,1:nlev) = pfull(1:npoints,1:nlev)/(rair*at(1:npoints,1:nlev))

    !Calculate umr, mass-weighted fall speed rain
    rhof(1:npoints,1:nlev) = (rhosu/rhoa(1:npoints,1:nlev))**0.54
    arn(1:npoints,1:nlev) = ar*rhof(1:npoints,1:nlev)
    asn(1:npoints,1:nlev) = at*rhof(1:npoints,1:nlev)
    acn(1:npoints,1:nlev) = ac*rhof(1:npoints,1:nlev)
    ain(1:npoints,1:nlev) = ai*rhof(1:npoints,1:nlev)

    do k=1,npoints	  !Loop over locations
           do ilev=1,nlev    !Loop over levels
	      if (mr_hydrolsr(k,ilev).gt.qsmall) then
	      	 lamr(k,ilev) = (pi*rhow*Niclsr(k,ilev)/mr_hydrolsr(k,ilev))**(1._wp/3._wp)
		 lamr(k,ilev) = max(lamr(k,ilev),lamminr)
		 lamr(k,ilev) = min(lamr(k,ilev),lammaxr)
		 mwfs_lsrain(k,ilev) = min(arn(k,ilev)*cons5/(6._wp*lamr(k,ilev)**br), 9.1*rhof(k,ilev))
	      endif

	      if (mr_hydrolss(k,ilev).gt.qsmall) then
	      	 lams(k,ilev) = (cons6*cs*Niclss(k,ilev)/mr_hydrolss(k,ilev))**(1._wp/ds)
		 lams(k,ilev) = max(lams(k,ilev),lammins)
		 lams(k,ilev) = min(lams(k,ilev),lammaxs)
		 mwfs_lssnow(k,ilev) = min(asn(k,ilev)*cons8/(6._wp*lams(k,ilev)**bs), 1.2*rhof(k,ilev))
	      endif

	      if (mr_hydrolsl(k,ilev).gt.qsmall) then
	         pgam(k,ilev) = 0.0005714*(Niclsl(k,ilev)/1E6*rhoa(k,ilev))+0.2714
		 pgam(k,ilev) = 1._wp/(pgam(k,ilev)**2._wp) - 1._wp
		 pgam(k,ilev) = max(pgam(k,ilev),2._wp)
		 pgam(k,ilev) = min(pgam(k,ilev),15._wp)
		 lamc(k,ilev) = (pi/6._wp*rhow*Niclsl(k,ilev)*gamma(pgam(k,ilev) + 4._wp) / &
		 (mr_hydrolsl(k,ilev) * gamma(pgam(k,ilev)+1._wp)))**(1._wp/3._wp)
		 lammin(k,ilev) = (pgam(k,ilev)+1._wp)/50E-6
		 lammax(k,ilev) = (pgam(k,ilev)+1._wp)/2E-6
		 lamc(k,ilev) = max(lamc(k,ilev),lammin(k,ilev))
		 lamc(k,ilev) = min(lamc(k,ilev),lammax(k,ilev))
		 mwfs_lsliq(k,ilev) = acn(k,ilev)*gamma(4._wp+bc+pgam(k,ilev))/(lamc(k,ilev)**bc &
		 *gamma(pgam(k,ilev)+4._wp)
	      endif

	      if (mr_hydrolsi(k,ilev).gt.qsmall) then
	      	 lami(k,ilev) = (cons1*ci*Niclsi(k,ilev)/mr_hydrolsi(k,ilev))**(1._wp/di)
		 lami(k,ilev) = max(lami(k,ilev),lammini)
		 lami(k,ilev) = min(lami(k,ilev),lammaxi)
		 mwfs_lsice(k,ilev) = ain(k,ilev)*cons17/(6._wp*lami(k,ilev)**bi)
		 mwfs_lsice(k,ilev) = min(mwfs_lsice(k,ilev),1.2*rhof(k,ilev))
	      endif
	         
	   enddo
    enddo
    
    ! Clean up space
    deallocate(rhoa,rhof,arn,asn,lamr,lams,lamc,pgam,lammin,lammax,acn,ain,lami)
	      	 

  end subroutine calc_mwfs
end module mod_mwfs
