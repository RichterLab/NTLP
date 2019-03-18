      subroutine tridv(b,d,a,r,n,j1,j2)
c
c --- tridiagonal matrix solver with multiple vectors
c     (note j and i loops are reversed from cray version)
c
c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c              j1:j2  range of input vectors
c
c --- output:  r   solution vector
c
      real b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)
c
      if(n .le. 1 ) then
         do j=j1,j2
            r(1,j) = r(1,j)/d(1,j)
         enddo
         go to 999
      endif
      do j=j1,j2
         d(1,j) = 1.0/d(1,j)
      enddo
      do j=j1,j2
      do i=2,n
         fac = b(i,j)*d(i-1,j)
         d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
         r(i,j) = r(i,j) - fac*r(i-1,j)
      enddo
      enddo
      do j=j1,j2
         r(n,j) = r(n,j)*d(n,j)
      enddo
      do j=j1,j2
      do i=n-1,1,-1
         r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      enddo
      enddo
  999 continue
c
      return
      end
