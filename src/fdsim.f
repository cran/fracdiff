      subroutine fdsim( n, ip, iq, ar, ma, d, rmu, y, s,
     *                  flmin, flmax, epmin, epmax)

      implicit double precision (a-h,o-z)

c  generates a random time series for use with fracdf
c
c  Input :
c
c  n      integer  length of the time series
c  ip     integer  number of autoregressive parameters
c  ar     real    (ip) autoregressive parameters
c  ma     real    (iq) moving average parameters
c  d      real     fractional differencing parameters
c  rmu    real     time series mean
c  y      real    (n+iq) 1st n : normalized random numbers
c  s      real    (n+iq) workspace
c
c  Output :
c
c  s      real   (n) the generated time series

c-----------------------------------------------------------------------------
c
c        Simulates a series of length n from an ARIMA (p,d,q) model
c        with fractional d (0 < d < 0.5). 
c
c-----------------------------------------------------------------------------

      integer            n, ip, iq
c     real               ar(ip), ma(iq), rmu, d
      double precision   ar(*), ma(*), rmu, d

      double precision   g0, vk, amk, sum, dk1, dk1d, dj, temp
c     real               y(n+iq), s(n+iq)
      double precision   y(*), s(*)

      double precision   flmin, flmax, epmin, epmax

      double precision   dgamr, dgamma

      external           dgamr, dgamma

      integer            k, j, i

      double precision   FLTMIN, FLTMAX, EPSMIN, EPSMAX
      common /MACHFD/    FLTMIN, FLTMAX, EPSMIN, EPSMAX
      save   /MACHFD/

      integer            IGAMMA, JGAMMA
      common /GAMMFD/    IGAMMA, JGAMMA
      save   /GAMMFD/

      real              zero, one, two
      parameter        (zero = 0.0, one = 1.0, two = 2.0)

*--------------------------------------------------------------------------

        IGAMMA = 0
        JGAMMA = 0

        FLTMIN  = flmin
        FLTMAX  = flmax
        EPSMIN  = epmin
        EPSMAX  = epmax
c
c	 Calculate g0

        temp = real(dgamr(dble(one-d)))
        if (IGAMMA .ne. 0) then
          do i = 1, n
            s(i) = zero
          end do
          return
        end if

        g0   = real(dgamma(dble(one-two*d)))*(temp*temp)
        if (IGAMMA .ne. 0) then
          do i = 1, n
            s(i) = zero
          end do
          return
        end if
c
c	 Generate y(1)
c
	y(1) = y(1)*sqrt(g0)
c
c	 Generate y(2) and initialise vk,phi(j)
c
	temp  = d / (one-d)
	vk    = g0*(one-(temp*temp))

	amk   = temp*y(1)
        s(1)  = temp
	y(2)  = amk + y(2)*sqrt(vk)
c
c	 Generate y(3),...,y(n+iq)
c
	do k = 3, n + iq
          dk1  = real(k) - one
          dk1d = dk1 - d
c
c	 Update the phi(j) using the recursion formula on W498
c
          do j = 1, k-2
            dj   = dk1 - real(j) 
            s(j) = s(j)*(dk1*(dj-d)/(dk1d*dj))
          end do

       	  temp   = d / dk1d
          s(k-1) = temp
c
c	 Update vk
c
	  vk = vk * (one-(temp*temp))
c
c	 Form amk
c
	  amk = zero
	  do j = 1, k-1
	    amk = amk + s(j)*y(k-j)
          end do
c
c	 Generate y(k)
c
	  y(k) = amk + y(k)*sqrt(vk)

	end do
c
c	 We now have an ARIMA (0,d,0) realisation of length n+iq in 
c	 y(k),k=1,n+iq. We now run this through an inverse ARMA(p,q)
c	 filter to get the final output in x(k),k=1,n.
c

	do k = 1, n 

	  sum = zero

          do i = 1, ip
	    if (k .le. i) go to 10
	    sum = sum + ar(i)*s(k-i)
          end do

10        continue

          do j = 1, iq
	    sum = sum-ma(j)*y(k+iq-j)
          end do

   	  s(k) = sum + y(k+iq)

        end do

        if (rmu .ne. zero) then
          do i = 1, n
            s(i) = s(i) + rmu
          end do
        end if

       return
       end
