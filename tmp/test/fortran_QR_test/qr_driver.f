      program main

      implicit none

      double precision A( 2000, 2000 ), t( 2000 ), work( 2000*104 )
      double precision Aold( 2000, 2000 )
      double precision b( 2000 ), x( 2000 )

      double precision extern rand, fla_clock_f, dnrm2

      double precision time, bnrm

      integer m, i, j, lda, lwork, info

      print *, "input m:"
      read *, m

      print *, "m = ", m

C     Create a random m x m matrix A and copy to Aold

      do j=1,m
         do i=1,m
            A( i, j )    = rand()
            Aold( i, j ) = A( i, j )
         enddo
      enddo

C     Create a random right-hand-side vector b and copy to x

      do i=1,m
         b( i ) = rand()
         x( i ) = b( i )
      enddo

      lda = 2000
      lwork = 2000*104

      time = fla_clock_f()

C     Factor the matrix

      call dgeqrf( m, m, A, lda, t, work, lwork, info )

C     Apply Q to x (=b)

      call dormqr( "L", "T", m, 1, m, A, lda, t, x, lda, work, lwork, 
     &             info )

C     Perform backward substitution

      call dtrsv( "U", "N", "N", m, a, lda, x, 1 )

      time = fla_clock_f() - time

C     Check residual

      call dgemv( "N", m, m, -1.0d0, Aold, lda, x, 1, 1.0d0, b, 1 )

C     print residual

      bnrm = dnrm2( m, b, 1 );
      print *, "norm of residual = ", bnrm
      print *, "time             = ", time
      print *, "MFLOPS           = ", 
     &     4.0 / 3.0 * m * m * m / time * 1e-6

      stop
      end
