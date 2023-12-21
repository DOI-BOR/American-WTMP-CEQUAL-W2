!  programs calcs temp stats
      
      parameter(nstart=1)
      parameter(ih=800000,ns=1)
      REAL JDAY
      
      real daymod(ih),pmod(ih,ns)
	real temp1(39),nh4n,nox

      common /stats/daymod,pmod,ncomp

      call sample_data_read
      
      OPEN(UNIT=10,FILE='two_34.csv',status='old')
	do i=1,3
	read(10,*)
	end do

      NITER=1

!c****START OF READ LOOP
      
200   CONTINUE

!c reading output file...                 
         read(10,*,end=502)JDAY,temp
     
5024  format(100g10.0)

     if(temp > 0.0)then  ! throwing out zero values

!c storing wq predictions
        daymod(niter)=jday

        
	  pmod(niter,1)=temp

        niter=niter+1
     end if

      GO TO 200     !!!!****** END OF READ LOOP

502   CONTINUE

!c calling subroutine calculating stats      
      ncomp=niter-1
      call wq_statistics

            STOP
            END


!******************************************************
!c  subroutine reads data files

      subroutine sample_data_read

      parameter(ih=800000,ns=1)

      real wlqdata(ih,ns),daydatacnt(ih)
	  real no2no3, nh3

      common /sampdata /daydatacnt,wlqdata,ndata
     

          open(25,file='T_11446500_AmerR_FairOaks.csv',status='old')	  

	 ! read(25,*)	
      !read(25,*)
      read(25,*)

        nit=1

 
23      read(25,*,end=200)daydatacnt(nit),temp
           wlqdata(nit,1)=temp


          nit=nit+1
        go to 23
200     continue

        ndata=nit-1

        close(25)

      return
      end

!******************************************************************
!c  subroutine calculates wq error statistics

      subroutine wq_statistics     

      parameter(ih=800000,ns=1)

      real wlqdata(ih,ns),daymod(ih),pmod(ih,ns)
      real merr(ns),daydatacnt(ih),abserr(ns)
	dimension std(ns)
    integer npt(ns)

      common /stats/daymod,pmod,ncomp
	common /sampdata /daydatacnt,wlqdata,ndata

	totmerr=0.0
	totstd=0.0
	nsite=0

	do is=1,ns

       taerr=0.0
       tsqe=0.0
       tmerr=0.0
       npt(is)=0
       
	
	  idata=ndata
	
       do id=1,idata
         
           daydata=daydatacnt(id)
		   if(wlqdata(id,is).le.-99.0)go to 345
	   

         do im=2,ncomp
           if(daydata.le.daymod(1))go to 345
           if(daydata.gt.daymod(ncomp))go to 345
           if(daymod(im).gt.daydata)then
             tdif=daymod(im)-daymod(im-1)
             w1=(daymod(im)-daydata)/tdif
             w2=(daydata-daymod(im-1))/tdif
             doint=w2*pmod(im,is)+w1*pmod(im-1,is)
             go to 25
           end if
         end do

25       npt(is)=npt(is)+1
         tmerr=tmerr+(doint-wlqdata(id,is))
         taerr=taerr+abs(doint-wlqdata(id,is))
         tsqe=tsqe+(doint-wlqdata(id,is))**2
         
345      continue
       end do

       if(npt(is).ne.0)then
         merr(is)=tmerr/real(npt(is))
         abserr(is)=taerr/real(npt(is))
         std(is)=sqrt(tsqe/real(npt(is)))
	   totmerr=totmerr+merr(is)
	   totstd=totstd+std(is)
	   nsite=nsite+1
	  else
         merr(is)=-999.0
         abserr(is)=-999.0
         std(is)=-999.0
	   
	  end if

	end do

!      1 : temperature
!	    2 : conductivity
!      3 : po4p
!	   4 : nh4n
!	   5 : nox
!	   6 : dox
!	  7 : ph
      

	open(99,file='stats_T_OUTFLOW.csv',status='unknown')
	write(99,344)
!344   format(2x,"Temp_AME_Error",2x,"Temp_RMS_error",  &
!			 2x,"cond_AME_Error",2x,"cond_RMS_error",  &
!			 2x,"po4p_AME_Error",2x,"po4p_RMS_error",  &
!			 2x,"nh4n_AME_Error",2x,"nh4n_RMS_error",  &
!             2x," nox_AME_Error",2x," nox_RMS_error",  &
!			 2x," dox_AME_Error",2x," dox_RMS_error",  &
!			 2x,"  pH_AME_Error",2x,"  pH_RMS_error")
344   format("Count",",","ME_Error",",","AME_Error",",","RMS_error",",")
!			 2x,"cond_AME_Error",2x,"cond_RMS_error",  &
!			 2x,"po4p_AME_Error",2x,"po4p_RMS_error",  &
!			 2x,"nh4n_AME_Error",2x,"nh4n_RMS_error",  &
!             2x," nox_AME_Error",2x," nox_RMS_error",  &
!			 2x," dox_AME_Error",2x," dox_RMS_error",  &
!			 2x,"  pH_AME_Error",2x,"  pH_RMS_error")
      	write(99,'(7(i8,",",f14.4,",",f14.4,",",f14.4,","))')(npt(is),merr(is),abserr(is),std(is), is=1,ns)
	
      

      return
      end



