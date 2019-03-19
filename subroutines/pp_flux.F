
c      !APR - Airborne particle reader
c      !Know when to take the reading
c      subroutine APR_update
c      use pars
c      use fields
c      use con_data
c      use con_stats
c      use fftwk
c      use particles
c#if defined(SWAP)
c      use module_byteswap
c#endif
c      include 'mpif.h'
c      character cgrid*4
c!     pseudocode:
c!     make sure to update and record locations equivalently as 
c!     the other particle recording information so things are in "phase"
c!     update aircraft particle location
c!     get particle count based off a specific "range"
c!     
c
c
c
c      end






c      subroutine root_gather
c      use pars
c      use fields
c      use con_data
c      use con_stats
c      use fftwk
c      use particles
c#if defined(SWAP)
c      use module_byteswap
c#endif
c      include 'mpif.h'
c
c      call probe_locator(probe_z_locations(1))
c      if(probe_id .eq. iploc)then
c      !Proceed to send u v w partcount information to root
c        call mpi_sendrecv(
c        
c      end if
c
c      call probe_locator(probe_z_locations(2))
c      if(probe_id .eq. iploc)then
c        call mpi_sendrecv(
c      end if
c
c      call probe_locator(probe_z_locations(3))
c      if(probe_id .eq. iploc)then
c        call mpi_sendrecv(
c      end if
c
c
c      end
c


       subroutine pp_flux(dxdy,type_f,iz,jpt,ipt)
       use particles
       use pars
       use con_data
       use con_stats
       implicit none
       include 'mpif.h'
       integer :: dxdy, x_p, y_p
       integer :: type_f !type of flux: 1 = up, 2 = down, 3 = total
       integer :: iz, jpt, ipt
       integer :: per_plane
       real    :: x_d, y_d
       real    :: xleft, xright
       real    :: yleft, yright

       per_plane = 4 !How many planes per one dimension?
                     !currently assuming a perfect square here.
       y_p = dxdy/per_plane + 1 !Figure out which plane located
       x_p = dxdy

       x_d = xl/per_plane
       y_d = yl/per_plane
       if(mod(dxdy,per_plane) .eq. 0)then
         y_p = dxdy/per_plane
       end if
       if(x_p .gt. per_plane)then
         x_p = dxdy-(y_p-1)*per_plane
       end if

       xleft  = x_d*(x_p-1)
       xright = x_d*(x_p)
       yleft  = y_d*(y_p-1)
       yright = y_d*(y_p)

       if(type_f .eq. 1)then !+1 up flux
         if(part%xp(1) .ge. xleft .and. part%xp(1) .lt. xright)then
           if(part%xp(2) .ge. yleft .and. part%xp(2) .lt. yright)then
             select case (dxdy)
               case(1)
                 pp1flux_t(iz,jpt,ipt) = pp1flux_t(iz,jpt,ipt) + 1.0
               case(2)
                 pp2flux_t(iz,jpt,ipt) = pp2flux_t(iz,jpt,ipt) + 1.0
               case(3)
                 pp3flux_t(iz,jpt,ipt) = pp3flux_t(iz,jpt,ipt) + 1.0
               case(4)
                 pp4flux_t(iz,jpt,ipt) = pp4flux_t(iz,jpt,ipt) + 1.0
               case(5)
                 pp5flux_t(iz,jpt,ipt) = pp5flux_t(iz,jpt,ipt) + 1.0
               case(6)
                 pp6flux_t(iz,jpt,ipt) = pp6flux_t(iz,jpt,ipt) + 1.0
               case(7)
                 pp7flux_t(iz,jpt,ipt) = pp7flux_t(iz,jpt,ipt) + 1.0
               case(8)
                 pp8flux_t(iz,jpt,ipt) = pp8flux_t(iz,jpt,ipt) + 1.0
               case(9)
                 pp9flux_t(iz,jpt,ipt) = pp9flux_t(iz,jpt,ipt) + 1.0
               case(10)
                 pp10flux_t(iz,jpt,ipt) = pp10flux_t(iz,jpt,ipt) + 1.0
               case(11)
                 pp11flux_t(iz,jpt,ipt) = pp11flux_t(iz,jpt,ipt) + 1.0
               case(12)
                 pp12flux_t(iz,jpt,ipt) = pp12flux_t(iz,jpt,ipt) + 1.0
               case(13)
                 pp13flux_t(iz,jpt,ipt) = pp13flux_t(iz,jpt,ipt) + 1.0
               case(14)
                 pp14flux_t(iz,jpt,ipt) = pp14flux_t(iz,jpt,ipt) + 1.0
               case(15)
                 pp15flux_t(iz,jpt,ipt) = pp15flux_t(iz,jpt,ipt) + 1.0
               case(16)
                 pp16flux_t(iz,jpt,ipt) = pp16flux_t(iz,jpt,ipt) + 1.0
!               case(17)
!                 pp17flux_t(iz,jpt,ipt) = pp17flux_t(iz,jpt,ipt) + 1.0
!               case(18)
!                 pp18flux_t(iz,jpt,ipt) = pp18flux_t(iz,jpt,ipt) + 1.0
!               case(19)
!                 pp19flux_t(iz,jpt,ipt) = pp19flux_t(iz,jpt,ipt) + 1.0
!               case(20)
!                 pp20flux_t(iz,jpt,ipt) = pp20flux_t(iz,jpt,ipt) + 1.0
!               case(21)
!                 pp21flux_t(iz,jpt,ipt) = pp21flux_t(iz,jpt,ipt) + 1.0
!               case(22)
!                 pp22flux_t(iz,jpt,ipt) = pp22flux_t(iz,jpt,ipt) + 1.0
!               case(23)
!                 pp23flux_t(iz,jpt,ipt) = pp23flux_t(iz,jpt,ipt) + 1.0
!               case(24)
!                 pp24flux_t(iz,jpt,ipt) = pp24flux_t(iz,jpt,ipt) + 1.0
!               case(25)
!                 pp25flux_t(iz,jpt,ipt) = pp25flux_t(iz,jpt,ipt) + 1.0
!               case(26)
!                 pp26flux_t(iz,jpt,ipt) = pp26flux_t(iz,jpt,ipt) + 1.0
!               case(27)
!                 pp27flux_t(iz,jpt,ipt) = pp27flux_t(iz,jpt,ipt) + 1.0
!               case(28)
!                 pp28flux_t(iz,jpt,ipt) = pp28flux_t(iz,jpt,ipt) + 1.0
!               case(29)
!                 pp29flux_t(iz,jpt,ipt) = pp29flux_t(iz,jpt,ipt) + 1.0
!               case(30)
!                 pp30flux_t(iz,jpt,ipt) = pp30flux_t(iz,jpt,ipt) + 1.0
!               case(31)
!                 pp31flux_t(iz,jpt,ipt) = pp31flux_t(iz,jpt,ipt) + 1.0
!               case(32)
!                 pp32flux_t(iz,jpt,ipt) = pp32flux_t(iz,jpt,ipt) + 1.0
!               case(33)
!                 pp33flux_t(iz,jpt,ipt) = pp33flux_t(iz,jpt,ipt) + 1.0
!               case(34)
!                 pp34flux_t(iz,jpt,ipt) = pp34flux_t(iz,jpt,ipt) + 1.0
!               case(35)
!                 pp35flux_t(iz,jpt,ipt) = pp35flux_t(iz,jpt,ipt) + 1.0
!               case(36)
!                 pp36flux_t(iz,jpt,ipt) = pp36flux_t(iz,jpt,ipt) + 1.0
!               case(37)
!                 pp37flux_t(iz,jpt,ipt) = pp37flux_t(iz,jpt,ipt) + 1.0
!               case(38)
!                 pp38flux_t(iz,jpt,ipt) = pp38flux_t(iz,jpt,ipt) + 1.0
!               case(39)
!                 pp39flux_t(iz,jpt,ipt) = pp39flux_t(iz,jpt,ipt) + 1.0
!               case(40)
!                 pp40flux_t(iz,jpt,ipt) = pp40flux_t(iz,jpt,ipt) + 1.0
!               case(41)
!                 pp41flux_t(iz,jpt,ipt) = pp41flux_t(iz,jpt,ipt) + 1.0
!               case(42)
!                 pp42flux_t(iz,jpt,ipt) = pp42flux_t(iz,jpt,ipt) + 1.0
!               case(43)
!                 pp43flux_t(iz,jpt,ipt) = pp43flux_t(iz,jpt,ipt) + 1.0
!               case(44)
!                 pp44flux_t(iz,jpt,ipt) = pp44flux_t(iz,jpt,ipt) + 1.0
!               case(45)
!                 pp45flux_t(iz,jpt,ipt) = pp45flux_t(iz,jpt,ipt) + 1.0
!               case(46)
!                 pp46flux_t(iz,jpt,ipt) = pp46flux_t(iz,jpt,ipt) + 1.0
!               case(47)
!                 pp47flux_t(iz,jpt,ipt) = pp47flux_t(iz,jpt,ipt) + 1.0
!               case(48)
!                 pp48flux_t(iz,jpt,ipt) = pp48flux_t(iz,jpt,ipt) + 1.0
!               case(49)
!                 pp49flux_t(iz,jpt,ipt) = pp49flux_t(iz,jpt,ipt) + 1.0
!               case(50)
!                 pp50flux_t(iz,jpt,ipt) = pp50flux_t(iz,jpt,ipt) + 1.0
!               case(51)
!                 pp51flux_t(iz,jpt,ipt) = pp51flux_t(iz,jpt,ipt) + 1.0
!               case(52)
!                 pp52flux_t(iz,jpt,ipt) = pp52flux_t(iz,jpt,ipt) + 1.0
!               case(53)
!                 pp53flux_t(iz,jpt,ipt) = pp53flux_t(iz,jpt,ipt) + 1.0
!               case(54)
!                 pp54flux_t(iz,jpt,ipt) = pp54flux_t(iz,jpt,ipt) + 1.0
!               case(55)
!                 pp55flux_t(iz,jpt,ipt) = pp55flux_t(iz,jpt,ipt) + 1.0
!               case(56)
!                 pp56flux_t(iz,jpt,ipt) = pp56flux_t(iz,jpt,ipt) + 1.0
!               case(57)
!                 pp57flux_t(iz,jpt,ipt) = pp57flux_t(iz,jpt,ipt) + 1.0
!               case(58)
!                 pp58flux_t(iz,jpt,ipt) = pp58flux_t(iz,jpt,ipt) + 1.0
!               case(59)
!                 pp59flux_t(iz,jpt,ipt) = pp59flux_t(iz,jpt,ipt) + 1.0
!               case(60)
!                 pp60flux_t(iz,jpt,ipt) = pp60flux_t(iz,jpt,ipt) + 1.0
!               case(61)
!                 pp61flux_t(iz,jpt,ipt) = pp61flux_t(iz,jpt,ipt) + 1.0
!               case(62)
!                 pp62flux_t(iz,jpt,ipt) = pp62flux_t(iz,jpt,ipt) + 1.0
!               case(63)
!                 pp63flux_t(iz,jpt,ipt) = pp63flux_t(iz,jpt,ipt) + 1.0
!               case(64)
!                 pp64flux_t(iz,jpt,ipt) = pp64flux_t(iz,jpt,ipt) + 1.0

             end select
           end if
         end if
       end if

       if(type_f .eq. 2)then !-1 down flux
         if(part%xp(1) .ge. xleft .and. part%xp(1) .lt. xright)then
           if(part%xp(2) .ge. yleft .and. part%xp(2) .lt. yright)then
             select case (dxdy)
               case(1)
                 pp1flux_t(iz,jpt,ipt) = pp1flux_t(iz,jpt,ipt) - 1.0
               case(2)
                 pp2flux_t(iz,jpt,ipt) = pp2flux_t(iz,jpt,ipt) - 1.0
               case(3)
                 pp3flux_t(iz,jpt,ipt) = pp3flux_t(iz,jpt,ipt) - 1.0
               case(4)
                 pp4flux_t(iz,jpt,ipt) = pp4flux_t(iz,jpt,ipt) - 1.0
               case(5)
                 pp5flux_t(iz,jpt,ipt) = pp5flux_t(iz,jpt,ipt) - 1.0
               case(6)
                 pp6flux_t(iz,jpt,ipt) = pp6flux_t(iz,jpt,ipt) - 1.0
               case(7)
                 pp7flux_t(iz,jpt,ipt) = pp7flux_t(iz,jpt,ipt) - 1.0
               case(8)
                 pp8flux_t(iz,jpt,ipt) = pp8flux_t(iz,jpt,ipt) - 1.0
               case(9)
                 pp9flux_t(iz,jpt,ipt) = pp9flux_t(iz,jpt,ipt) - 1.0
               case(10)
                 pp10flux_t(iz,jpt,ipt) = pp10flux_t(iz,jpt,ipt) - 1.0
               case(11)
                 pp11flux_t(iz,jpt,ipt) = pp11flux_t(iz,jpt,ipt) - 1.0
               case(12)
                 pp12flux_t(iz,jpt,ipt) = pp12flux_t(iz,jpt,ipt) - 1.0
               case(13)
                 pp13flux_t(iz,jpt,ipt) = pp13flux_t(iz,jpt,ipt) - 1.0
               case(14)
                 pp14flux_t(iz,jpt,ipt) = pp14flux_t(iz,jpt,ipt) - 1.0
               case(15)
                 pp15flux_t(iz,jpt,ipt) = pp15flux_t(iz,jpt,ipt) - 1.0
               case(16)
                 pp16flux_t(iz,jpt,ipt) = pp16flux_t(iz,jpt,ipt) - 1.0
!               case(17)
!                 pp17flux_t(iz,jpt,ipt) = pp17flux_t(iz,jpt,ipt) - 1.0
!               case(18)
!                 pp18flux_t(iz,jpt,ipt) = pp18flux_t(iz,jpt,ipt) - 1.0
!               case(19)
!                 pp19flux_t(iz,jpt,ipt) = pp19flux_t(iz,jpt,ipt) - 1.0
!               case(20)
!                 pp20flux_t(iz,jpt,ipt) = pp20flux_t(iz,jpt,ipt) - 1.0
!               case(21)
!                 pp21flux_t(iz,jpt,ipt) = pp21flux_t(iz,jpt,ipt) - 1.0
!               case(22)
!                 pp22flux_t(iz,jpt,ipt) = pp22flux_t(iz,jpt,ipt) - 1.0
!               case(23)
!                 pp23flux_t(iz,jpt,ipt) = pp23flux_t(iz,jpt,ipt) - 1.0
!               case(24)
!                 pp24flux_t(iz,jpt,ipt) = pp24flux_t(iz,jpt,ipt) - 1.0
!               case(25)
!                 pp25flux_t(iz,jpt,ipt) = pp25flux_t(iz,jpt,ipt) - 1.0
!               case(26)
!                 pp26flux_t(iz,jpt,ipt) = pp26flux_t(iz,jpt,ipt) - 1.0
!               case(27)
!                 pp27flux_t(iz,jpt,ipt) = pp27flux_t(iz,jpt,ipt) - 1.0
!               case(28)
!                 pp28flux_t(iz,jpt,ipt) = pp28flux_t(iz,jpt,ipt) - 1.0
!               case(29)
!                 pp29flux_t(iz,jpt,ipt) = pp29flux_t(iz,jpt,ipt) - 1.0
!               case(30)
!                 pp30flux_t(iz,jpt,ipt) = pp30flux_t(iz,jpt,ipt) - 1.0
!               case(31)
!                 pp31flux_t(iz,jpt,ipt) = pp31flux_t(iz,jpt,ipt) - 1.0
!               case(32)
!                 pp32flux_t(iz,jpt,ipt) = pp32flux_t(iz,jpt,ipt) - 1.0
!               case(33)
!                 pp33flux_t(iz,jpt,ipt) = pp33flux_t(iz,jpt,ipt) - 1.0
!               case(34)
!                 pp34flux_t(iz,jpt,ipt) = pp34flux_t(iz,jpt,ipt) - 1.0
!               case(35)
!                 pp35flux_t(iz,jpt,ipt) = pp35flux_t(iz,jpt,ipt) - 1.0
!               case(36)
!                 pp36flux_t(iz,jpt,ipt) = pp36flux_t(iz,jpt,ipt) - 1.0
!               case(37)
!                 pp37flux_t(iz,jpt,ipt) = pp37flux_t(iz,jpt,ipt) - 1.0
!               case(38)
!                 pp38flux_t(iz,jpt,ipt) = pp38flux_t(iz,jpt,ipt) - 1.0
!               case(39)
!                 pp39flux_t(iz,jpt,ipt) = pp39flux_t(iz,jpt,ipt) - 1.0
!               case(40)
!                 pp40flux_t(iz,jpt,ipt) = pp40flux_t(iz,jpt,ipt) - 1.0
!               case(41)
!                 pp41flux_t(iz,jpt,ipt) = pp41flux_t(iz,jpt,ipt) - 1.0
!               case(42)
!                 pp42flux_t(iz,jpt,ipt) = pp42flux_t(iz,jpt,ipt) - 1.0
!               case(43)
!                 pp43flux_t(iz,jpt,ipt) = pp43flux_t(iz,jpt,ipt) - 1.0
!               case(44)
!                 pp44flux_t(iz,jpt,ipt) = pp44flux_t(iz,jpt,ipt) - 1.0
!               case(45)
!                 pp45flux_t(iz,jpt,ipt) = pp45flux_t(iz,jpt,ipt) - 1.0
!               case(46)
!                 pp46flux_t(iz,jpt,ipt) = pp46flux_t(iz,jpt,ipt) - 1.0
!               case(47)
!                 pp47flux_t(iz,jpt,ipt) = pp47flux_t(iz,jpt,ipt) - 1.0
!               case(48)
!                 pp48flux_t(iz,jpt,ipt) = pp48flux_t(iz,jpt,ipt) - 1.0
!               case(49)
!                 pp49flux_t(iz,jpt,ipt) = pp49flux_t(iz,jpt,ipt) - 1.0
!               case(50)
!                 pp50flux_t(iz,jpt,ipt) = pp50flux_t(iz,jpt,ipt) - 1.0
!               case(51)
!                 pp51flux_t(iz,jpt,ipt) = pp51flux_t(iz,jpt,ipt) - 1.0
!               case(52)
!                 pp52flux_t(iz,jpt,ipt) = pp52flux_t(iz,jpt,ipt) - 1.0
!               case(53)
!                 pp53flux_t(iz,jpt,ipt) = pp53flux_t(iz,jpt,ipt) - 1.0
!               case(54)
!                 pp54flux_t(iz,jpt,ipt) = pp54flux_t(iz,jpt,ipt) - 1.0
!               case(55)
!                 pp55flux_t(iz,jpt,ipt) = pp55flux_t(iz,jpt,ipt) - 1.0
!               case(56)
!                 pp56flux_t(iz,jpt,ipt) = pp56flux_t(iz,jpt,ipt) - 1.0
!               case(57)
!                 pp57flux_t(iz,jpt,ipt) = pp57flux_t(iz,jpt,ipt) - 1.0
!               case(58)
!                 pp58flux_t(iz,jpt,ipt) = pp58flux_t(iz,jpt,ipt) - 1.0
!               case(59)
!                 pp59flux_t(iz,jpt,ipt) = pp59flux_t(iz,jpt,ipt) - 1.0
!               case(60)
!                 pp60flux_t(iz,jpt,ipt) = pp60flux_t(iz,jpt,ipt) - 1.0
!               case(61)
!                 pp61flux_t(iz,jpt,ipt) = pp61flux_t(iz,jpt,ipt) - 1.0
!               case(62)
!                 pp62flux_t(iz,jpt,ipt) = pp62flux_t(iz,jpt,ipt) - 1.0
!               case(63)
!                 pp63flux_t(iz,jpt,ipt) = pp63flux_t(iz,jpt,ipt) - 1.0
!               case(64)
!                 pp64flux_t(iz,jpt,ipt) = pp64flux_t(iz,jpt,ipt) - 1.0

             end select
           end if
         end if
       end if

       
       end 
