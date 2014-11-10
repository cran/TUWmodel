!      MODULE HBVMODEL
!
!      contains
!
      subroutine hbvmodel(itsteps,nzones,area,param,incon,  
     *                 prec,airt,ep,output)
!     ----- shell for HBV model -------------------------------------

!      USE PARAMETERS 

      PARAMETER(imodincon=4, maxout=20, maxparam=15)

      integer itsteps,nzones,hzbnr

!      real*8 param(maxparam),incon(imodincon),area
      real*8 param(maxparam,nzones),incon(imodincon,nzones),area
	real*8 prec(itsteps),airt(itsteps),ep(itsteps),snowd(itsteps)
	real*8 output

      real*8 csf,ddf,tr,ts,meltt,sweq(itsteps),lprat,lp,fc,beta,
     *       k0,k1,k2,lsuz,cperc,bmax,croute

      real*8 temp,precip,swe,rain,snow,melt
	real*8 etp,dmoist,moist,dq,eta
      real*8 suz,slz,dquh,qg,q0,q1,q2,q
	integer bql,age

      real*8 aa,ddfmin,ddfmax,ddfage

	real*8 sdold,chng
	integer*4 nasim(nzones),nasim2(nzones),nasim3(nzones)

      dimension dquh(itsteps),q(itsteps)
      dimension area(nzones)
!	dimension input(4*nzones,itsteps)  !input(4*maxzones,maxdays) 
	dimension output(nzones,maxout,itsteps)      
	
	logical assimilation
	common /logic/ calibration,assimilation,hzbnr


!     ---- set up parameters ----
!	open(88,file="hbv.log")
!	write(88,*) nzones
      do izone=1,nzones
      csf=param(1,izone)	
      ddf=param(2,izone)	
	tr=param(3,izone)
	ts=param(4,izone)
	meltt=param(5,izone) 

	lprat=param(6,izone)
	FC=param(7,izone)	
	BETA=param(8,izone)

	LP=lprat*fc

	k0=param(9,izone)	
	k1=param(10,izone)	
	k2=param(11,izone)
	lsuz=param(12,izone)
			
	cperc=param(13,izone)
	bmax=param(14,izone)	
	croute=param(15,izone)	
	
!	write(88,*) csf,ddf
	ddfmin=1.0
	ddfmax=6.0
	ddfage=ddf

      output=0.  
	nasim=0
	nasim2=0
	nasim3=0

!      do izone=1,nzones
	age=0
	moist=incon(1,izone) 
      swe=incon(2,izone)   
      suz=incon(3,izone)   
      slz=incon(4,izone) 
	sdold=0. !for data assimilation, initial % of snow cover
	
!	write(88,*) swe, suz, slz, area(izone)

      if(area(izone).gt.0.0) then
  	 do i=1,itsteps
	  q(i)=0.
	  dquh(i)=0.
	 end do
       aa=area(izone)
       do it=1,itsteps
		precip=prec(it)
		temp=airt(it)
		etp=ep(it)
		sd=snowd(it)
	
!	if (it .lt.10) write(88,*) it, precip, temp, etp

	if (temp.lt.-0.1)then
	etp=0.
	endif

       if(precip.lt.-998.00) then
	  output(izone,1,it)=-999.99 
	  output(izone,2,it)=-999.99
	  output(izone,3,it)=-999.99
	  output(izone,4,it)=-999.99
	  output(izone,5,it)=-999.99
	  output(izone,6,it)=-999.99
          output(izone,7,it)=-999.99
          output(izone,8,it)=-999.99
          output(izone,9,it)=-999.99
          output(izone,10,it)=-999.99
        goto  330
       endif

        call snowmod(csf,ddf,tr,ts,meltt,temp,precip,swe,
     *       rain,snow,melt)

	
!Calculation of refreezing process 

!Refreezing of the liquid water is included in the model. It is defined in this case by the diurnal average temperature so that, below the threshold temperature TB,F, refreezing MF (positive value, mm d-1)is calculated according to the following formula: 

 !mf=kf*(tbf-ta)^ef

!where TB,F (°C), refreezing parameter KF (mm °C-1 d-1) and empirical exponent eF have to be calibrated. According to Vehviläinen (1992), average values and standard deviations of the parameters are as follows: TB,F –1.7 and 1.2 (°C), KF 1.5 and 1.4 and eF 0.36 and 0.44. 



!	if(assimilation) then
!

!	if ((swe.gt.0.1).and.(sd.le.0.1)) then
!	nasim3(izone)=nasim3(izone)+1
!	swe=0.
!	endif

!	if ((swe .le. 0.1) .and. (sd.gt.0.1))then !assimilation by snow accumulation
!	if ((swe .le. 0.0) .and. (sd.gt.0.0))then !assimilation by snow accumulation
!	swe=precip
!	if (precip.le.0.1) then
!	print*,hzbnr,swe,sd,sdold,precip,melt,rain,snow
!	stop
!	swe=melt*0.2
!	melt=melt*0.8
!	endif
	
!	if (swe.le.0.1) then
!	swe=melt
!	melt=0.
!	endif

!	if (swe.lt.0.1) swe=0.11

!	nasim(izone)=nasim(izone)+1
!	if (chng.ge.0.)nasim2(izone)=nasim2(izone)+1
	 
!	endif
	


!	end if
	sdold=sd

        call soilmoisture(rain,melt,etp,LP,FC,Beta,dmoist,moist,
     *       dq,eta) 

         call respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
     *	slz,bql,dquh,qg,q0,q1,q2)

        do 300,irf=1,bql
         if((it+irf-1).gt.itsteps) goto 301
         q(it+irf-1)=q(it+irf-1)+dquh(irf)
300     continue
301     continue

!        q(it)=q(it)/24./3.6*area(izone)
!        qg=qg/24./3.6*area(izone)
!        q0=q0/24./3.6*area(izone)
!        q1=q1/24./3.6*area(izone)
!        q2=q2/24./3.6*area(izone)

	 output(izone,1,it)=q(it) 
	 output(izone,2,it)=swe
	 sweq(it)=swe      ! added by al. 
	 output(izone,3,it)=moist
	 output(izone,4,it)=rain
	 output(izone,5,it)=snow
	 output(izone,6,it)=melt
         output(izone,7,it)=q0
         output(izone,8,it)=q1
         output(izone,9,it)=q2
         output(izone,10,it)=eta

!      write(88,'(100f12.4)') rain,snow,melt,swe,moist,q(it)
!	write(88,'(100f12.4)') izone*1.,it*1.,q(it),output(izone,1,it), 
!     *   swe,output(izone,2,it)
!      write(88,*) '#############################'
!	pause
       enddo
330   continue
      else
	 do it=1,itsteps
	  output(izone,1,it)=0.
	  output(izone,2,it)=0. 
	  output(izone,3,it)=0.
	  output(izone,4,it)=0.
	  output(izone,5,it)=0. 
	  output(izone,6,it)=0.
          output(izone,7,it)=0.
          output(izone,8,it)=0.
          output(izone,9,it)=0.
          output(izone,10,it)=0.

	 enddo
      endif
	enddo
	
!	if(assimilation) then
!	write(37,'(100I8)') hzbnr,(nasim(i),i=1,nzones),
!     *	(nasim2(ii),ii=1,nzones),(nasim3(ii),ii=1,nzones)
!	write(37,'(100I8)') hzbnr,(nasim(i),i=1,nzones)
!	print*,hzbnr,nasim
!	pause
!	endif
!	close(88)
      return

      end subroutine hbvmodel





      subroutine respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
     *	slz,bql,dquh,qg,q0,q1,q2)

      parameter(maxday=15000)
      integer bql
      real*8 dq
      real*8 k0,lsuz,k1,k2,cperc
      real*8 bmax,croute
      real*8 suz,slz,suzold,slzold
      real*8 q0,q1,q2,qg
      real*8 dquh(maxday)

      dt=1.

!     dt=1.
!     ----- new ---
      rat=1.0
      suzold=suz+rat*dq
      slzold=slz+(1.-rat)*dq
      slzin=cperc

      if(suzold.lt.0.) suzold=0.
      if(slzold.lt.0.) slzold=0.


!     --- upper storage ---
      if(suzold.gt.lsuz) then
       q0=(suzold-lsuz)/k0*exp(-1./k0)
       if(q0.lt.0.) q0=0.
       if(q0.gt.(suzold-lsuz)) q0=suzold-lsuz
      else
       q0=0.
      endif
      suzold=suzold-q0 
      q1=-slzin+(slzin+suzold/k1)*exp(-1./k1)
      if(q1.lt.0.0) q1=0.
      suz=suzold-q1-slzin
      if(suz.lt.0.) then
       suz=0.
       slzin=suzold
       endif

!     --- lower storage ---    
       q2=slzin-(slzin-slzold/k2)*exp(-1./k2)
       if(q2.lt.0.) q2=0.
       slz=slzold-q2+slzin
        if(slz.lt.0.) then
         slz=0.
         q2=slzold+slzin
        endif
      
      
      qg=q0+q1+q2

!     --- transformation function ---
      if((bmax-croute*qg).gt.1.0) then 
       bq=bmax-croute*qg
       bql=INT(bq)
       sum=0.
       do 400, j=1,bql
        if(j.le.bql/2) then
         dquh(j)=((j-0.5)*4.*qg)/(bql*bql*1.)
        elseif(abs(j-(bql/2.+0.5)).lt.0.1) then
         dquh(j)=((j-0.75) *4.*qg)/(bql*bql*1.)
        else
         dquh(j)=((bql-j+0.5)*4.*qg)/(bql*bql*1.)
        endif
         sum=sum+dquh(j)
400    continue
      else
       bql=1
       dquh(1)=qg
       sum=qg 
      endif
      return 
      
	end subroutine respfunc

      subroutine soilmoisture(rain,melt,etp,LP,FC,Beta,dmoist,moist,
     *	dq,eta)

      real*8 rain,melt,lp,beta,fc,moistold,moist
      real*8 dq,dmoist
      real*8 etp,eta

!     --- soil mositure accounting ----
      moistold=moist
      dq=((moistold/FC)**beta)*(rain+melt)     
      if(dq.gt.(rain+melt)) dq=rain+melt
      dmoist=rain+melt-dq
      if(dmoist.lt.0.) dmoist=0.
      moist=moistold+dmoist
      if(moist.gt.fc) then
        dq=(moist-fc)+dq
        moist=fc
      endif
!     --- calculate evapotranspiration ---
      if(moist.lt.LP) then
       ETA=moist*ETP/LP 
       if(eta.gt.etp) then
               eta=etp    
       end if
	else
       ETA=ETP
      endif
      if(eta.lt.0.) eta=0.
!     --- substract eta of soilmoisture ---
      xx=moist
      moist=moist-eta
      if(moist.lt.0.) then
       eta=xx     
       moist=0.
      endif
      return
      
	end subroutine soilmoisture

      subroutine snowmod(csf,ddf,tr,ts,melttemp,temp,precip,swe,
     *	rain,snow,melt)

       real*8 ddf,csf,tr,ts,melttemp,dsrtemp,dd1
       real*8 temp,precip,rain,snow,sweold,swe,melt     

      if(temp.lt.ts) then
	 snow=precip
      elseif(temp.gt.tr) then
	 snow=0.
      else
       snow=precip*abs(temp-tr)/abs(tr-ts)
	endif
      rain=precip-snow


       melt=(temp-melttemp)*ddf 
       if (melt .lt. 0.) melt=0.
!      --- Bestimme SWE
       sweold=swe
       swe=sweold+csf*snow-melt
       if (swe.lt.0.0001) then
          swe=0.
          melt=sweold+csf*snow
	    if(melt.lt.0.) melt=0.
       end if   
	 return
	end subroutine snowmod

!	END MODULE HBVMODEL

