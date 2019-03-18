      subroutine range(n1,n2,nprocs,irank,ista,iend)
c
c ---------- the ibm range finder to balance load
c
      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if(iwork2 .gt. irank) iend = iend + 1
c
      return
      end
