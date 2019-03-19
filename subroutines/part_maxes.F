      
cc     This subroutine shifts the wind speeds and subsequently
cc     the concentration flux of aerosols based on Monahan's whitecap
cc     data (not yet for concentration flux)
c      subroutine shifting_conditions
c      use pars
c      use fields
c      use con_stats
c
c      
c      if(myid==0)then
c        write(*,*)'Old ug:'
c        do iz=1,nnz
c          write(*,*) ug(iz)
c        enddo
c      end if
c      do iz=1,nnz
c         ug(iz)   = ugcont_new
c         vg(iz)   = vgcont_new 
c      enddo 
c      if(myid==0)then
c        write(*,*)'New ug:'
c        do iz=1,nnz
c          write(*,*) ug(iz)
c        enddo
c      end if
c      end
c      
      subroutine part_maxes
      use pars
      use con_stats
      use particles 

      real :: components(3),magnitude

      part => first_particle
      do while (associated(part))
         components=(part%vp(1:3)+part%u_sub(1:3))
         magnitude=sqrt(components(1)**2+
     +             components(2)**2+components(3)**2)
         if(magnitude .ge. part%max_speed)then
           part%max_speed = magnitude
         end if

         if(part%xp(3) .ge. part%z_highest)then
           part%z_highest = part%xp(3)
         end if
         part => part%next
      end do



      end
