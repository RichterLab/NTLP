      subroutine get_dt
c
c ---------- routine computes max time step for given cfl number
c            from max's found previously
c
      use pars
      use con_data
      use con_stats
c

      ucfl = umax
      vcfl = vmax
      wcfl = wmax
      ucflm = ucfl
      vcflm = vcfl
      wcflm = wcfl
      vel_max = wcflm
      vel_max = amax1(ucflm,vel_max)
      vel_max = amax1(vcflm,vel_max)

      if(vel_max .le. 0.0) then
          write(6,6000) ucflm, vcflm,wcflm, vel_max
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +           ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +           ' vel_max = ',e15.5,/,
     +           ' infinite time step !!!')
          stop
      endif
c
c ---------------- choose fixed or variable time step
c
c       write(*,*) 'dt_con', dt_con
      if(ifix_dt .ne. 0) then
c
c ------------- if used, change to fit your problem
c
        dt_new = 0.5
        dt_new = dt_con
      else
c
c ------------------- new estimate of best time step
c                     from cfl constraint
c
      dt_new  = cfl/vel_max
      dt_new = amin1(dt_new, 5.0)
c     dt_new = amin1(dt_new, 10.0)
      endif

c ---------------- compare against viscous stability limit
c
      if(vismax*dt_new .gt. 0.5) then
         dt_cfl = dt_new
         dt_new = 0.5/vismax
         if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +      '   viscous time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' vismax = ',e15.6)
         endif
      endif
c
c -------- for safety if restart set timestep = saved timestep in
c          restart file
c
      if(it_step .eq. iti .and. iti .ne. 0) then
        dt_new = dt
      endif

c
      return
      end
