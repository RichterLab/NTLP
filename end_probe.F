
       subroutine end_probe(p1,p2)
       use pars
       use fields
       use con_stats
       integer :: p1, p2
!      p1 represents probe_id (izs<=probe_id<=ize)
!      p2 is the do loop for p3=1:3
         call probe_locator(probe_z_locations(p2))
         if(iploc .eq. p1)then
           write(*,*)'proc, p1, zz, izs, iploc,ize: ',
     +           myid,p1,zz(iploc),izs,iploc,ize
           write(*,*)'id,w(0),w(1),w(2),w(3),w(4),w(5)',
     +           myid,w(1,1,0),w(1,1,1),w(1,1,2),w(1,1,3),
     +           w(1,1,4),w(1,1,5)
         end if

       end
