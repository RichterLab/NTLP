      subroutine set_sav(it,istart)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ------------------- increment time if not first time through
c
         time=time+dt
      endif
c
      it=it+1
c
      dt    = dt_new
      dtg   = dt
      mnout = (mod(it,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it,itape).eq.0)
      micut = (mod(it,itcut).eq.0)
      mviz  = (mod(it,i_viz).eq.0)
      mphist = (mod(it,iphist).eq.0)
      
      if(ihst .lt. 0) then
         mhis = .false.
      else
         mhis = (mod(it,ihst).eq.0 .and. it .ge. it_his)
      endif
c     mtrans - transpose for the particle quantities from z to x
       mtrans = (mod(it+1,ihst).eq.0)

      if (i_viz .lt. 0) then
         msave_v = .false.
      else
         msave_v = (mod(it,i_viz).eq.0 .and. it .ge. it_viz)
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(it .ge. itstr .and. mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c ---------- decide whether viz fields are saved
c
      !msave_v = .false.
      !if(it .ge. itstr .and. mviz .and. i_viz .gt. 0) then
      !   msave_v = .true.
      !   if(ionce .eq. 0) then
      !      ionce = 1
      !      call open_viz
      !   endif
      !endif
      if((i_viz .gt. 0) .and. (it .ge. it_viz_nxt)) then
        call viz_output_filename(it)
        it_viz_nxt = it_viz_nxt + itape
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it .ge. it_his_nxt)) then
         call open_his(it)
         it_his_nxt = it_his_nxt + itape
      endif
c
c ------decide whether the particle histogrm is saved
c
      part_save = .false.
      if(it .ge. itstr .and. mphist) then   
         part_save = .true.
      endif
 
c
c ------probe reading specifications (3 different locations)
c

      mprobe = (mod(it,it_probe).eq.0 .and. it .ge. itstr)

      if((it .ge. itstr) .and. (it .ge. it_pro_nxt)) then
c        call open_probe(it)
        it_pro_nxt = it_pro_nxt + itape
      endif

      return
      end
