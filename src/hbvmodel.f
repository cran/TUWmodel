!      MODULE HBVMODEL
!
!      contains
!
      subroutine hbvmodel(itsteps,nzones,area,param,incon,  
     *                  input,output)
!     *                 prec,airt,ep,output)
!     *                 prec,airt,ep,q,sweq)
!     ----- shell for HBV model -------------------------------------

!      USE PARAMETERS 

      PARAMETER(imodincon=4, maxout=20, maxparam=15)

      integer itsteps,nzones,hzbnr

      real*4 param(maxparam),incon(imodincon),area
      real*4 prec(itsteps),airt(itsteps),ep(itsteps),snowd(itsteps)
       real*8 input(4*nzones,itsteps),output(nzones,maxout,itsteps)

      real*4 csf,ddf,tr,ts,meltt,lprat,lp,fc,beta,
     *       k0,k1,k2,lsuz,cperc,bmax,croute

      real*4 temp,precip,swe,rain,snow,melt
      real*4 etp,dmoist,moist,dq,eta
      real*4 suz,slz,dquh,qg,q0,q1,q2,q
      integer bql,age

      real*4 aa,ddfmin,ddfmax,ddfage

      real*4 sdold,chng
      integer*4 nasim(nzones),nasim2(nzones),nasim3(nzones)

      dimension dquh(itsteps),q(itsteps)
      dimension area(nzones)
!      dimension input(4*nzones,itsteps)  !input(4*maxzones,maxdays) 
!      dimension output(nzones,maxout,itsteps)      
      
      logical assimilation
      common /logic/ calibration,assimilation,hzbnr


!	open(66, file='/home/alviglio/hbv.log')
!     ---- set up parameters ----
!	write(66, *) area
!	write(66, *) nzones
      csf=param(1)      
      ddf=param(2)      
      tr=param(3)
      ts=param(4)
      meltt=param(5) 

      lprat=param(6)
      FC=param(7)      
      BETA=param(8)

      LP=lprat*fc

      k0=param(9)      
      k1=param(10)      
      k2=param(11)
      lsuz=param(12)
                  
      cperc=param(13)
      bmax=param(14)      
      croute=param(15)      
      
      ddfmin=1.0
      ddfmax=6.0
      ddfage=ddf

      output=0.  
      nasim=0
      nasim2=0
      nasim3=0

      do izone=1,nzones
      age=0
      moist=incon(1) 
      swe=incon(2)   
      suz=incon(3)   
      slz=incon(4) 
      sdold=0. !for data assimilation, initial % of snow cover
      

      if(area(izone).gt.0.0) then
         do i=1,itsteps
        q(i)=0.
        dquh(i)=0.
       end do
       aa=area(izone)
       do it=1,itsteps
!            precip=prec(it)
            precip=input(((izone-1)*4)+1,it)
!            temp=airt(it)
            temp=input(((izone-1)*4)+2,it)
!            etp=ep(it)
            etp=input(((izone-1)*4)+3,it)
!            sd=snowd(it)
            sd=input(((izone-1)*4)+4,it)
!	write(66, *) izone, it, precip, temp      
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
        goto  330
       endif

        call snowmod(csf,ddf,tr,ts,meltt,temp,precip,swe,
     *       rain,snow,melt)

      
!Calculation of refreezing process 

!Refreezing of the liquid water is included in the model. It is defined in this case by the diurnal average temperature so that, below the threshold temperature TB,F, refreezing MF (positive value, mm d-1)is calculated according to the following formula: 

!mf=kf*(tbf-ta)^ef

!where TB,F (°C), refreezing parameter KF (mm °C-1 d-1) and empirical exponent eF have to be calibrated. According to Vehviläinen (1992), average values and standard deviations of the parameters are as follows: TB,F –1.7 and 1.2 (°C), KF 1.5 and 1.4 and eF 0.36 and 0.44. 



!      if(assimilation) then
!

!      if ((swe.gt.0.1).and.(sd.le.0.1)) then
!      nasim3(izone)=nasim3(izone)+1
!      swe=0.
!      endif

!      if ((swe .le. 0.1) .and. (sd.gt.0.1))then !assimilation by snow accumulation
!      if ((swe .le. 0.0) .and. (sd.gt.0.0))then !assimilation by snow accumulation
!      swe=precip
!      if (precip.le.0.1) then
!      print*,hzbnr,swe,sd,sdold,precip,melt,rain,snow
!      stop
!      swe=melt*0.2
!      melt=melt*0.8
!      endif
      
!      if (swe.le.0.1) then
!      swe=melt
!      melt=0.
!      endif

!      if (swe.lt.0.1) swe=0.11

!      nasim(izone)=nasim(izone)+1
!      if (chng.ge.0.)nasim2(izone)=nasim2(izone)+1
       
!      endif
      


!      end if
      sdold=sd

        call soilmoisture(rain,melt,etp,LP,FC,Beta,dmoist,moist,
     *       dq,eta) 

         call respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
     *      slz,bql,dquh,qg,q0,q1,q2)

        do 300,irf=1,bql
         if((it+irf-1).gt.itsteps) goto 301
         q(it+irf-1)=q(it+irf-1)+dquh(irf)
300     continue
301     continue
        q(it)=q(it)/24./3.6*area(izone)
        qg=qg/24./3.6*area(izone)
        q0=q0/24./3.6*area(izone)
        q1=q1/24./3.6*area(izone)
        q2=q2/24./3.6*area(izone)

       output(izone,1,it)=q(it) 
!	write(66,*) it, precip, q(it), output(izone,1,it)
       output(izone,2,it)=swe
       output(izone,3,it)=moist/FC
       output(izone,4,it)=rain
       output(izone,5,it)=snow
       output(izone,6,it)=melt
       output(izone,7,it)=q0
       output(izone,8,it)=q1
       output(izone,9,it)=q2
       output(izone,10,it)=eta
       output(izone,11,it)=suz
       output(izone,12,it)=slz


!      write(66,'(100f12.4)') rain,snow,melt,swe,moist,q(it)
!      write(88,'(100f12.4)') izone*1.,it*1.,q(it),output(izone,1,it), 
!     *   swe,output(izone,2,it)
!      write(88,*) '#############################'
!      pause
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


       enddo
      endif
      enddo
      
!      if(assimilation) then
!      write(37,'(100I8)') hzbnr,(nasim(i),i=1,nzones),
!     *      (nasim2(ii),ii=1,nzones),(nasim3(ii),ii=1,nzones)
!      write(37,'(100I8)') hzbnr,(nasim(i),i=1,nzones)
!      print*,hzbnr,nasim
!      pause
!      endif

!	close(66)
      return

      end subroutine hbvmodel





      subroutine respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
     *      slz,bql,dquh,qg,q0,q1,q2)

      parameter(maxday=15000)
      integer bql
      real*4 dq
      real*4 k0,lsuz,k1,k2,cperc
      real*4 bmax,croute
      real*4 suz,slz,suzold,slzold
      real*4 q0,q1,q2,qg
      real*4 dquh(maxday)

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
     *      dq,eta)

      real*4 rain,melt,lp,beta,fc,moistold,moist
      real*4 dq,dmoist
      real*4 etp,eta

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
     *      rain,snow,melt)

       real*4 ddf,csf,tr,ts,melttemp,dsrtemp,dd1
       real*4 temp,precip,rain,snow,sweold,swe,melt     

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

!      END MODULE HBVMODEL
