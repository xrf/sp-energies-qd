      program yalla
      implicit none
      integer :: i, j, n
      real(kind=8) :: hw, a, b, c, d, dmc
      real(kind=8) :: elog1, elog2, elog3, corre, hf
      real(kind=8), allocatable, dimension(:) :: r, x1, x2, x3

      read(5,*)n, hw, dmc
      allocate(r(n),x1(n),x2(n), x3(n))
      do i=1,n
         read(5,*) j, a, b, c,d
         r(i) = j
         x1(i) = b*hw
         x2(i) = c*hw
         x3(i) = d*hw
      enddo 
      hf = a*hw
      corre = abs(a*hw-dmc)
!      dmc = x3(n)

      do i=1,n
         a= (((x1(i)-hf)-corre))/corre
         b = (((x2(i)-hf)-corre))/corre
         c = (((x3(i)-hf)-corre))/corre
         elog1 = 100-abs(a)*100
         elog2 = 100-(abs( b ))*100
         elog3 =  100-(abs( c ))*100
         write(6,'(4E20.10)') r(i), elog1, elog2, elog3
      enddo 


      end
