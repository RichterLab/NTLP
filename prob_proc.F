      
      subroutine prob_proc(nprb,iplocd)
      use pars
      use fields
      use con_stats
   
      character(len=500) :: fname
!     Locate proc that is responsible for izs and ize
      if(izs .le. iplocd .and. ize .ge. iplocd)then
!     Currently set probe location to y=0 (arbitrary for now)
        if(iys .eq. 1) then
          write(*,*),'proc,zz,izs,iplocd,ize:',
     +            myid,zz(iplocd),izs,iplocd,ize
          if(myid==0)then
            write(*,*)'wbc: ',w(1,1,0)
          end if
          probe_loc = myid !assign processor id to global variable
          probe_id  = iplocd !
          write(fname,"(A,I0.3,A)")'probe',probe_loc,'.out'
          open(nprb,file=fname,form='unformatted',
     +            status='replace')
        end if
      end if
      
      
      
      
      end
