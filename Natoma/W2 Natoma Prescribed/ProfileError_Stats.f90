integer :: iend,run,run_number,proNum34,n_TEMP,n_TEMP34
real :: error_me,error_ame,error_rms,seg
character*18 infile
character*8, dimension (4) :: title
real, dimension (100000) :: JDAY,TEMP34,depth,junk,depth2
real, dimension (100000) :: dseg,djday,ddata,ddepth,ddate

real :: me_TEMP,ame_TEMP,rms_TEMP,me_ave_Temp,ame_ave_Temp,rms_ave_Temp
real :: me_TEMP34,ame_TEMP34,rms_TEMP34,me_ave_TEMP34,ame_ave_TEMP34,rms_ave_TEMP34

real, dimension (100000) :: seg2,jday2,elev2,me2,ame2,rms2
character*2 rnum

open(40,file='NatomaProfileData.txt',status='old')
!open(41,file='TempData_Seg13.txt',status='old')
!open(42,file='TempData_Seg14.txt',status='old')
!open(43,file='TempData_Seg18.txt',status='old')


run=1

      write(*,*)'Run number:'
      read(*,*)run_number

write(rnum,'(i2)')run_number
rnum=adjustl(trim(rnum))
open(150,file='prof_error_all'//rnum//'.txt',status='unknown')

open(200,file='spr_1.opt',status='old')

!	read(10,*)
!12	if (run<=7) then
!		read(10,'(a15,4x,f5.1,f3.0)',end=200)infile,depth,seg
!	else if (run<=13) then
!		read(10,'(a16,3x,f5.1,f3.0)',end=200)infile,depth,seg
!	else
!		read(10,'(a17,2x,f5.1,f4.0)',end=200)infile,depth,seg
!	end if
!
!	open(200,file=infile,status='old')

n_TEMP34=0
n_TEMP=0

    read(200,*)

do i=1,100000
	read(200,'(41x,f7.3,f10.3,15x,f5.2)',end=450)JDAY(i),Depth(i),TEMP34(i)
	
    if (TEMP34(i)>40) then
        TEMP34(i)=-99
        else
    end if

end do

450 continue

iend=i-1

		READ(40,*)
		READ(40,*)
		
		do k=1,100000

		read(40,*,end=300)djday(k),junk(k),ddepth(k),ddata(k),ddate(k)
		dseg(k)=34
			do i=1,iend
					if (JDAY(i)<=(djday(k)+0.6) .and. JDAY(i)>=(djday(k)-0.6)) then
					
						if (ddepth(k)<=(depth(i)+0.5) .and. ddepth(k)>=(depth(i)-0.5)) then
						if(TEMP34(i)>0.0) then
						    if(ddepth(k)>depth(i))then
						      if(i==1 .or. TEMP34(i+1)==-99.00)then
						        t=TEMP34(i)
						        else
							t=TEMP34(i)-((TEMP34(i)-TEMP34(i+1))/(depth(i)-depth(i+1)))*(depth(i)-ddepth(k))
							    end if
								error_me=(t-ddata(k))
								error_ame=abs(t-ddata(k))
								error_rms=(t-ddata(k))**2
								write(150,'(f10.0,2x,f10.3,2x,f10.1,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)')dseg(k),djday(k),ddepth(k),error_me,error_ame,error_rms,t,ddata(k)
								exit
								
							    else
						    if(i==1)then
						            t=TEMP34(i)
						        else if (TEMP34(i-1)==-99.00) then
						            t=TEMP34(i)
						        else
							        t=TEMP34(i)+((TEMP34(i)-TEMP34(i-1))/(depth(i)-depth(i-1)))*(depth(i)-ddepth(k))
							end if
								error_me=(t-ddata(k))
								error_ame=abs(t-ddata(k))
								error_rms=(t-ddata(k))**2
								write(150,'(f10.0,2x,f10.3,2x,f10.1,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)')dseg(k),djday(k),ddepth(k),error_me,error_ame,error_rms,t,ddata(k)
								exit
						    end if
						else
						end if
						end if	
					end if
             end do
         end do
         
300 continue
           

rewind(150)

proNum34=1
me_Temp=0
ame_Temp=0
rms_Temp=0
n_Temp=0
me_Temp34=0
ame_Temp34=0
rms_Temp34=0
n_Temp34=0


do i=1,100000
	read(150,'(f10.0,f12.3,f12.1,f12.4,f12.4,f12.4)',end=555)seg2(i),jday2(i),depth2(i),me2(i),ame2(i),rms2(i)
		
				me_Temp=me_Temp+me2(i)
				ame_Temp=ame_Temp+ame2(i)
				rms_Temp=rms_Temp+rms2(i)
				n_Temp=n_Temp+1
				
				if (seg2(i)==34) then
				    me_TEMP34=me_TEMP34+me2(i)
				    ame_TEMP34=ame_TEMP34+ame2(i)
				    rms_TEMP34=rms_TEMP34+rms2(i)
				    n_TEMP34=n_TEMP34+1
				        
				        if (i==1) then
				                go to 999
				        else
				            if (jday2(i)/=jday2(i-1))then
				                if (jday2(i)<=jday2(i-1)+0.6 .and. jday2(i)>=jday2(i-1)-0.6)then
				                    go to 999
				                else
				                    proNum34=proNum34+1
				                end if
				            else
				            end if
				        end if
				
999                     continue                       
                end if
                end do

555 continue

	me_ave_Temp=me_Temp/n_Temp
	ame_ave_Temp=ame_Temp/n_Temp
	rms_ave_Temp=sqrt(rms_Temp/n_Temp)

	me_ave_TEMP34=me_TEMP34/n_TEMP34
	ame_ave_TEMP34=ame_TEMP34/n_TEMP34
	rms_ave_TEMP34=sqrt(rms_TEMP34/n_TEMP34)
	




open(99,file='ProfFinalStats'//rnum//'.txt',status='unknown')

write(99,*)'Overall Model Data Comparison Statistics'
write(99,'(5x,3a14,5x,a6,5x,a12,2x,i2)')'ME','AME','RMS','Num','Run number:',run_number
write(99,'(a10,3f14.4,i)')'Temp',me_ave_Temp,ame_ave_Temp,rms_ave_Temp,n_Temp

write(99,*)
write(99,*)
write(99,*)'Model Statistics by Segment'
write(99,'(10x,3a14,6x,a6,2x,12a)')'ME','AME','RMS','Num','# Profiles'
write(99,'(a10,3f14.4,i,i)')'Segment 34',me_ave_TEMP34,ame_ave_TEMP34,rms_ave_TEMP34,n_TEMP34,proNum34


end
