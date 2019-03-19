      function ranf()
      data inc /1/
      save inc, ix, ia, m, fm
      if(inc.eq.1) then
        inc = 2
        m = 2**20
        fm = float(m)
        ix = 566387
        ia = 2**10 + 3
      endif
      ix = mod(ia*ix,m)
      fx = float(ix)
      ranf = fx/fm
      return
      end
