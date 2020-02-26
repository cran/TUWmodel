!      MODULE HBVMODEL
! contains dual soil layers for implementation with
! satellite soil moisture data
!
!      contains
!
      subroutine hbvmodel_dual(itsteps,nzones,area,param,incon,  
     *                 prec,airt,ep,output)
!     ----- shell for HBV model -------------------------------------

!      USE PARAMETERS 

      PARAMETER(imodincon=5, maxout=20, maxparam=18)

      integer itsteps,nzones,hzbnr

      real(8) param(maxparam,nzones),incon(imodincon,nzones),area
      real(8) prec(itsteps,nzones),airt(itsteps,nzones)
      real(8) ep(itsteps,nzones),snowd(itsteps,nzones)
      real(8) output

      real(8) csf,ddf,tr,ts,meltt,sweq(itsteps),lprat,lp,fc,beta,
     *       k0,k1,k2,lsuz,cperc,bmax,croute

      real(8) temp,precip,swe,rain,snow,melt
      real(8) etp,dmoist,moist,dq,eta,sd
      real(8) suz,slz,dquh,qg,q0,q1,q2,q
      integer bql,age

      real(8) aa,ddfmin,ddfmax,ddfage

      real(8) sdold
      real(8) fc_skin,k_moist,f_eta,moist_skin  ! new 3-Nov-2018

!      real(8) chng
!      integer(4) nasim(nzones),nasim2(nzones),nasim3(nzones)

      dimension dquh(itsteps),q(itsteps)
      dimension area(nzones)

      dimension output(nzones,maxout,itsteps)      

      logical assimilation
      common /logic/ calibration,assimilation,hzbnr


!     ---- set up parameters ----
!	open(88,file="hbv.log")
!	write(88,*) nzones

      output=0.  

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

        fc_skin=param(16,izone)   ! new 3-Nov-2018
        f_eta=param(17,izone)     ! new 3-Nov-2018
        k_moist=param(18,izone)   ! new 3-Nov-2018

!	write(88,*) csf,ddf
        ddfmin=1.0
        ddfmax=6.0
        ddfage=ddf


        age=0
        moist=incon(1,izone) 
        swe=incon(2,izone)   
        suz=incon(3,izone)   
        slz=incon(4,izone) 
        moist_skin=incon(5,izone)   ! new 3-Nov-2018
        sdold=0. !for data assimilation, initial % of snow cover

      if(area(izone).gt.0.0) then
         do i=1,itsteps
          q(i)=0.
          dquh(i)=0.
         end do
       aa=area(izone)
       do it=1,itsteps
                precip=prec(it,izone)
                temp=airt(it,izone)
                etp=ep(it,izone)
                sd=snowd(it,izone)

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
          output(izone,11,it)=-999.99
          output(izone,12,it)=-999.99
          output(izone,13,it)=-999.99
        goto  330
       endif

        call snowmod(csf,ddf,tr,ts,meltt,temp,precip,swe,
     *       rain,snow,melt)


!Calculation of refreezing process 

!Refreezing of the liquid water is included in the model. It is defined in this case by the diurnal average temperature so that, below the threshold temperature TB,F, refreezing MF (positive value, mm d-1)is calculated according to the following formula: 

!mf=kf*(tbf-ta)^ef

!where TB,F (°C), refreezing parameter KF (mm °C-1 d-1) and empirical exponent eF have to be calibrated. According to Vehviläinen (1992), average values and standard deviations of the parameters are as follows: TB,F –1.7 and 1.2 (°C), KF 1.5 and 1.4 and eF 0.36 and 0.44. 


!        call soilmoisture(rain,melt,etp,LP,FC,Beta,dmoist,moist,
!     *       dq,eta) 

       call soilmoisture_skin(rain,melt,etp,LP,FC,Beta,dmoist,moist,   ! new 3-Nov-2018 
     *       dq,eta,FC_skin,f_eta,k_moist,moist_skin)                  ! new 3-Nov-2018


        call respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
     *       slz,bql,dquh,qg,q0,q1,q2)

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
         output(izone,11,it)=suz    ! added by al.
         output(izone,12,it)=slz    ! added by al.
         output(izone,13,it)=moist_skin   ! new 3-Nov-2018

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
          output(izone,11,it)=0.
          output(izone,12,it)=0.
          output(izone,13,it)=0.   ! new 3-Nov-2018
        enddo
      endif
        enddo

!	close(88)
      return

      end subroutine hbvmodel_dual





!----- already in hbvmodel.f -----!           subroutine respfunc(dq,k0,lsuz,k1,k2,cperc,bmax,croute,suz,
!----- already in hbvmodel.f -----!          *           slz,bql,dquh,qg,q0,q1,q2)
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           parameter(maxday=15000)
!----- already in hbvmodel.f -----!           integer bql
!----- already in hbvmodel.f -----!           real(8) dq
!----- already in hbvmodel.f -----!           real(8) k0,lsuz,k1,k2,cperc
!----- already in hbvmodel.f -----!           real(8) bmax,croute,bq
!----- already in hbvmodel.f -----!           real(8) suz,slz,suzold,slzold,slzin
!----- already in hbvmodel.f -----!           real(8) q0,q1,q2,qg,sum
!----- already in hbvmodel.f -----!           real(8) dquh(maxday)
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           dt=1.
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     !     dt=1.
!----- already in hbvmodel.f -----!     !     ----- new ---
!----- already in hbvmodel.f -----!           rat=1.0
!----- already in hbvmodel.f -----!           suzold=suz+rat*dq
!----- already in hbvmodel.f -----!           slzold=slz+(1.-rat)*dq
!----- already in hbvmodel.f -----!           slzin=cperc
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           if(suzold.lt.0.) suzold=0.
!----- already in hbvmodel.f -----!           if(slzold.lt.0.) slzold=0.
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     !     --- upper storage ---
!----- already in hbvmodel.f -----!           if(suzold.gt.lsuz) then
!----- already in hbvmodel.f -----!            q0=(suzold-lsuz)/k0*exp(-1./k0)
!----- already in hbvmodel.f -----!            if(q0.lt.0.) q0=0.
!----- already in hbvmodel.f -----!            if(q0.gt.(suzold-lsuz)) q0=suzold-lsuz
!----- already in hbvmodel.f -----!           else
!----- already in hbvmodel.f -----!            q0=0.
!----- already in hbvmodel.f -----!           endif
!----- already in hbvmodel.f -----!           suzold=suzold-q0 
!----- already in hbvmodel.f -----!           q1=-slzin+(slzin+suzold/k1)*exp(-1./k1)
!----- already in hbvmodel.f -----!           if(q1.lt.0.0) q1=0.
!----- already in hbvmodel.f -----!           suz=suzold-q1-slzin
!----- already in hbvmodel.f -----!           if(suz.lt.0.) then
!----- already in hbvmodel.f -----!            suz=0.
!----- already in hbvmodel.f -----!            slzin=suzold
!----- already in hbvmodel.f -----!            endif
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     !     --- lower storage ---    
!----- already in hbvmodel.f -----!            q2=slzin-(slzin-slzold/k2)*exp(-1./k2)
!----- already in hbvmodel.f -----!            if(q2.lt.0.) q2=0.
!----- already in hbvmodel.f -----!            slz=slzold-q2+slzin
!----- already in hbvmodel.f -----!             if(slz.lt.0.) then
!----- already in hbvmodel.f -----!              slz=0.
!----- already in hbvmodel.f -----!              q2=slzold+slzin
!----- already in hbvmodel.f -----!             endif
!----- already in hbvmodel.f -----!           
!----- already in hbvmodel.f -----!           
!----- already in hbvmodel.f -----!           qg=q0+q1+q2
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     !     --- transformation function ---
!----- already in hbvmodel.f -----!           if((bmax-croute*qg).gt.1.0) then 
!----- already in hbvmodel.f -----!            bq=bmax-croute*qg
!----- already in hbvmodel.f -----!            bql=INT(bq)
!----- already in hbvmodel.f -----!            sum=0.
!----- already in hbvmodel.f -----!            do 400, j=1,bql
!----- already in hbvmodel.f -----!             if(j.le.bql/2) then
!----- already in hbvmodel.f -----!              dquh(j)=((j-0.5)*4.*qg)/(bql*bql*1.)
!----- already in hbvmodel.f -----!             elseif(abs(j-(bql/2.+0.5)).lt.0.1) then
!----- already in hbvmodel.f -----!              dquh(j)=((j-0.75) *4.*qg)/(bql*bql*1.)
!----- already in hbvmodel.f -----!             else
!----- already in hbvmodel.f -----!              dquh(j)=((bql-j+0.5)*4.*qg)/(bql*bql*1.)
!----- already in hbvmodel.f -----!             endif
!----- already in hbvmodel.f -----!              sum=sum+dquh(j)
!----- already in hbvmodel.f -----!     400    continue
!----- already in hbvmodel.f -----!           else
!----- already in hbvmodel.f -----!            bql=1
!----- already in hbvmodel.f -----!            dquh(1)=qg
!----- already in hbvmodel.f -----!            sum=qg 
!----- already in hbvmodel.f -----!           endif
!----- already in hbvmodel.f -----!           return 
!----- already in hbvmodel.f -----!           
!----- already in hbvmodel.f -----!           end subroutine respfunc
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           subroutine soilmoisture(rain,melt,etp,LP,FC,Beta,dmoist,moist,
!----- already in hbvmodel.f -----!          *           dq,eta)
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           real(8) rain,melt,lp,beta,fc,moistold,moist,xx
!----- already in hbvmodel.f -----!           real(8) dq,dmoist
!----- already in hbvmodel.f -----!           real(8) etp,eta
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     !     --- soil mositure accounting ----
!----- already in hbvmodel.f -----!           moistold=moist
!----- already in hbvmodel.f -----!           dq=((moistold/FC)**beta)*(rain+melt)     
!----- already in hbvmodel.f -----!           if(dq.gt.(rain+melt)) dq=rain+melt
!----- already in hbvmodel.f -----!           dmoist=rain+melt-dq
!----- already in hbvmodel.f -----!           if(dmoist.lt.0.) dmoist=0.
!----- already in hbvmodel.f -----!           moist=moistold+dmoist
!----- already in hbvmodel.f -----!           if(moist.gt.fc) then
!----- already in hbvmodel.f -----!             dq=(moist-fc)+dq
!----- already in hbvmodel.f -----!             moist=fc
!----- already in hbvmodel.f -----!           endif
!----- already in hbvmodel.f -----!     !     --- calculate evapotranspiration ---
!----- already in hbvmodel.f -----!           if(moist.lt.LP) then
!----- already in hbvmodel.f -----!            ETA=moist*ETP/LP 
!----- already in hbvmodel.f -----!            if(eta.gt.etp) then
!----- already in hbvmodel.f -----!                    eta=etp    
!----- already in hbvmodel.f -----!            end if
!----- already in hbvmodel.f -----!             else
!----- already in hbvmodel.f -----!            ETA=ETP
!----- already in hbvmodel.f -----!           endif
!----- already in hbvmodel.f -----!           if(eta.lt.0.) eta=0.
!----- already in hbvmodel.f -----!     !     --- substract eta of soilmoisture ---
!----- already in hbvmodel.f -----!           xx=moist
!----- already in hbvmodel.f -----!           moist=moist-eta
!----- already in hbvmodel.f -----!           if(moist.lt.0.) then
!----- already in hbvmodel.f -----!            eta=xx     
!----- already in hbvmodel.f -----!            moist=0.
!----- already in hbvmodel.f -----!           endif
!----- already in hbvmodel.f -----!           return
!----- already in hbvmodel.f -----!           
!----- already in hbvmodel.f -----!           end subroutine soilmoisture


      ! New 3-Nov-2018
      subroutine soilmoisture_skin(rain,melt,etp,LP,FC,Beta,dmoist,
     *  moist,dq,eta, FC_skin, f_eta,k_moist, moist_skin)

      real(8) rain,melt,lp,beta,fc,moistold,moist
      real(8) dq,dmoist
      real(8) etp,eta
      real(8) moistold_skin, moist_skin, FC_skin, dQ_skin, by, bypass
      real(8) eta_skin, f_eta, Q_moist,moist_diff,k_moist

!     --- soil mositure accounting ----
      moistold=moist
      moistold_skin=moist_skin
      
      moist_skin=moistold_skin+rain+melt
      if (moist_skin .gt. FC_skin) then
             dQ_skin=moist_skin-FC_skin
             moist_skin=FC_skin
      else
             dQ_skin=0.0
      endif

      dq=((moistold/FC)**beta)*dQ_skin

!		  if (by.gt.0.) then
!	   if (moistold/FC.gt.0.4.and.moistold/FC.lt.0.9) then
!	    bypass=by*dQ_skin
!		if (bypass.gt.0.05) bypass=0.05
!	    else
!	    bypass=0.
!	   endif
!	  else 
!	    bypass=0.
!	  endif

  
      if(dq.gt.dQ_skin) dq=dQ_skin
  
!      dmoist=dQ_skin-dq-bypass
      dmoist=dQ_skin-dq

      if(dmoist.lt.0.) dmoist=0.
      moist=moistold+dmoist
      if(moist.gt.FC) then
        dq=(moist-FC)+dq
        moist=FC
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
      eta_skin=eta*f_eta
      eta=eta-eta_skin
      if (eta_skin .gt. moist_skin) then
           eta=eta+(eta_skin-moist_skin)
           eta_skin=moist_skin
      endif


!     --- substract eta of soilmoisture ---
      xx=moist
      moist=moist-eta
      moist_skin=moist_skin-eta_skin

      if(moist.lt.0.) then
         eta=xx     
         moist=0.
      endif

      moist_diff=moist_skin/FC_skin-moist/FC
      Q_moist=moist_diff*k_moist


      if (Q_moist .ge. 0.) then
          if (Q_moist .le. moist_skin) then
               moist=moist+abs(Q_moist)
               if (moist .gt. FC) then
                     moist_skin=moist_skin-abs(Q_moist)+(moist-FC)
                     moist=FC
               else
                     moist_skin=moist_skin-abs(Q_moist)
               endif
          else
               moist=moist+moist_skin
               if (moist .gt. FC) then
                     moist_skin=(moist-FC)
                     moist=FC
               else
                     moist_skin=0.
               endif
          endif
      else
          if (abs(Q_moist) .le. moist) then
               moist_skin=moist_skin+abs(Q_moist)
               if (moist_skin .gt. FC_skin) then
                     moist=moist-abs(Q_moist)+(moist_skin-FC_skin)
                     moist_skin=FC_skin
               else
                     moist=moist-abs(Q_moist)
               endif
          else
               moist_skin=moist_skin+moist
               if (moist_skin .gt. FC_skin) then
                     moist=(moist_skin-FC_skin)
                     moist_skin=FC_skin
               else
                     moist=0.
               endif
          endif
      endif  

      return
      
      end subroutine soilmoisture_skin





!----- already in hbvmodel.f -----!           subroutine snowmod(csf,ddf,tr,ts,melttemp,temp,precip,swe,
!----- already in hbvmodel.f -----!          *                   rain,snow,melt)
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!            real(8) ddf,csf,tr,ts,melttemp
!----- already in hbvmodel.f -----!     !       real(8) dsrtemp,dd1
!----- already in hbvmodel.f -----!            real(8) temp,precip,rain,snow,sweold,swe,melt     
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!           if(temp.lt.ts) then
!----- already in hbvmodel.f -----!              snow=precip
!----- already in hbvmodel.f -----!           elseif(temp.gt.tr) then
!----- already in hbvmodel.f -----!              snow=0.
!----- already in hbvmodel.f -----!           else
!----- already in hbvmodel.f -----!            snow=precip*abs(temp-tr)/abs(tr-ts)
!----- already in hbvmodel.f -----!             endif
!----- already in hbvmodel.f -----!           rain=precip-snow
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!     
!----- already in hbvmodel.f -----!            melt=(temp-melttemp)*ddf 
!----- already in hbvmodel.f -----!            if (melt .lt. 0.) melt=0.
!----- already in hbvmodel.f -----!     !      --- Bestimme SWE
!----- already in hbvmodel.f -----!            sweold=swe
!----- already in hbvmodel.f -----!            swe=sweold+csf*snow-melt
!----- already in hbvmodel.f -----!            if (swe.lt.0.0001) then
!----- already in hbvmodel.f -----!               swe=0.
!----- already in hbvmodel.f -----!               melt=sweold+csf*snow
!----- already in hbvmodel.f -----!                if(melt.lt.0.) melt=0.
!----- already in hbvmodel.f -----!            end if   
!----- already in hbvmodel.f -----!            return
!----- already in hbvmodel.f -----!             end subroutine snowmod

!	END MODULE HBVMODEL
