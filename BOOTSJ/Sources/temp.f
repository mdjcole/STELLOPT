
      function temp (x, T)
      use stel_kinds
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real(rprec) :: x, temp, T(0:10)
C-----------------------------------------------
      temp = (T(0) + x*(T(1) + x*(T(2) + x*(T(3) + x*(T(4) +
     1     x*(T(5) + x*(T(6) + x*(T(7) + x*(T(8) + x*(T(9) + x*T(10)
     2   ))))))))))
      end function temp
