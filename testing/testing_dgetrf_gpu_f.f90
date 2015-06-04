!
!   -- MAGMA (version 1.4.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      December 2013
!
!  @generated d Tue Dec 17 13:18:56 2013
!
      program testing_dgetrf_gpu_f

      use magma

      external cublas_init, cublas_set_matrix, cublas_get_matrix
      external cublas_shutdown, cublas_alloc
      external dlange, dgemm, dgesv, dlamch

      double precision dlange, dlamch
      integer cublas_alloc

      double precision              :: rnumber(2), Anorm, Bnorm, Rnorm, Xnorm 
      double precision, allocatable :: work(:)
      double precision, allocatable       :: h_A(:), h_B(:), h_X(:)
      magma_devptr_t                :: devptrA, devptrB
      integer,    allocatable       :: ipiv(:)

      double precision                    :: zone, mzone
      integer                       :: i, n, info, stat, lda, ldda
      integer                       :: size_of_elt, nrhs
      real(kind=8)                  :: flops, t, tstart, tend

      PARAMETER          ( nrhs = 1, zone = 1., mzone = -1. )
      
      call cublas_init()

      n   = 2048
      lda  = n
      ldda = ((n+31)/32)*32
      size_of_elt = sizeof_double
 
!------ Allocate CPU memory
      allocate(h_A(lda*n))
      allocate(h_B(lda*nrhs))
      allocate(h_X(lda*nrhs))
      allocate(work(n))
      allocate(ipiv(n))

!------ Allocate GPU memory
      stat = cublas_alloc(ldda*n, size_of_elt, devPtrA)
      if (stat .ne. 0) then
         write(*,*) "device memory allocation failed"
         stop
      endif

      stat = cublas_alloc(ldda*nrhs, size_of_elt, devPtrB)
      if (stat .ne. 0) then
         write(*,*) "device memory allocation failed"
         stop
      endif

!---- Initializa the matrix
      do i=1,lda*n
         call random_number(rnumber)
         h_A(i) = rnumber(1)
      end do

      do i=1,lda*nrhs
        call random_number(rnumber)
        h_B(i) = rnumber(1)
      end do
      h_X(:) = h_B(:)

!---- devPtrA = h_A
      call cublas_set_matrix(n, n, size_of_elt, h_A, lda, devptrA, ldda)

!---- devPtrB = h_B
      call cublas_set_matrix(n, nrhs, size_of_elt, h_B, lda, devptrB, ldda)

!---- Call magma LU ----------------
      call magma_wtime_f(tstart)
      call magmaf_dgetrf_gpu(n, n, devptrA, ldda, ipiv, info)
      call magma_wtime_f(tend)

      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if

!---- Call magma solve -------------
      call magmaf_dgetrs_gpu('n', n, nrhs, devptrA, ldda, ipiv, devptrB, ldda, info)

      if ( info .ne. 0 )  then
         write(*,*) "Info : ", info
      end if

!---- h_X = devptrB
      call cublas_get_matrix (n, nrhs, size_of_elt, devptrB, ldda, h_X, lda)

!---- Compare the two results ------
      Anorm = dlange('I', n, n,    h_A, lda, work)
      Bnorm = dlange('I', n, nrhs, h_B, lda, work)
      Xnorm = dlange('I', n, nrhs, h_X, lda, work)
      call dgemm('n', 'n', n,  nrhs, n, zone, h_A, lda, h_X, lda, mzone, h_B, lda)
      Rnorm = dlange('I', n, nrhs, h_B, lda, work)

      write(*,*)
      write(*,*  ) 'Solving A x = b using LU factorization:'
      write(*,105) '  || A || = ', Anorm
      write(*,105) '  || b || = ', Bnorm
      write(*,105) '  || x || = ', Xnorm
      write(*,105) '  || b - A x || = ', Rnorm

      flops = 2. * n * n * n / 3.
      t = tend - tstart

      write(*,*)   '  Gflops  = ',  flops / t / 1e9
      write(*,*)

      Rnorm = Rnorm / ( (Anorm*Xnorm+Bnorm) * n * dlamch('E') )

      if ( Rnorm > 60. ) then
         write(*,105) '  Solution is suspicious, ', Rnorm
      else
         write(*,105) '  Solution is CORRECT' 
      end if

!---- Free CPU memory
      deallocate(h_A, h_X, h_B, work, ipiv)

!---- Free GPU memory
      call cublas_free(devPtrA)
      call cublas_free(devPtrB)
      call cublas_shutdown()

 105  format((a35,es10.3))

      end
