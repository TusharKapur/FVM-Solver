  !1-D nozzle inviscid
program oneDnozzle
implicit none

       double precision, parameter :: m_dot=1.0
       double precision, parameter :: rho=1.0
       integer, parameter :: nx=1000, n = nx-2
     double precision, parameter :: L=2.00
     double precision, parameter :: dx=L/(nx-1)
     double precision ar(2*nx), au(nx-1,3), a(nx-2,3), u(nx-1), Di(nx), F(nx), &
      S(nx-1), d(nx-1), po, p(nx), pcr(nx), ust(nx-1), Fstw(nx-2), Fste(nx-2), Dst(nx-2), b(n), &
      c(n,n), pst(nx), m_res(nx-1), p_res(nx-1), ur, diag(n), lowdiag(n-1), updiag(n-1), x(n), &
      uru, urp

     double precision :: ua=2.0
     integer :: i, j, nt, time, nrhs, info, lda, ldb
     integer, dimension(n) :: ipiv

     !double precision, dimension (:,:), allocatable :: pmat
     !allocate ( pmat(n,n) )

     ! Under Relaxation Factor
     ur = 0.7
     time = 1
     uru=0.7; urp=0.7

     
     
     
     arallocationloop: do i = 2, 2*nx

                    ar(i) = 0.5 - i*(0.5-0.1)*dx/(2*L)
                   ! ar(2*i-1) = 0.6 - i*(0.5-0.1)*dx/2.0

     end do arallocationloop
    ar(1) = 0.5
    ar(2*nx-1) = 0.1
     
     

     ! U Initial conditions
     usetloop: do i = 1,nx-1

            ust(i) = m_dot/(rho*ar(2*i))
            
     end do usetloop

     ! P Initial conditions
     po = 10.0
     pst(1) = po


     psetloop: do i = 2,nx

            pst(i) = pst(1) - (i-1)*(10-7.5)*dx/0.5

     end do psetloop
     u = ust
     p = pst
!     print* , ust


timeloop: do nt = 1,1000

     ! U constants
     F(1) = rho*ust(1)*ar(1)
     
     if(time==1) then
            F(nx) = 1.0
     else
            F(nx) = rho*ar(2*nx-2)*ust(nx-1)
     end if

     fdsloop: do i = 2,nx-1

            F(i) = rho*(ust(i)+ust(i-1))*ar(2*i-1)/2
            Di(i) = 0  ! 0.0000181*nx/L  ! Tau/(x1-x2)
            S(i) = (pst(i)-pst(i+1))*ar(2*i)

     end do fdsloop

     auloop: do i = 2,nx-1

            au(i,1) = Di(i) + max(F(i),0.0)
            au(i,3) = Di(i+1) + max(-F(i+1),0.0)
            au(i,2) = au(i,1) + au(i,3) + F(i+1) - F(i)

     end do auloop

     ! Left BCS
     ua = ust(1)*ar(2)/ar(1)
     F(1) = rho*ua*ar(1)
     au(1,:) = 0.0
     au(1,2) = F(2) + F(1)*((ar(2)/ar(1))**2)/2
     S(1) = (pst(1)-pst(2))*ar(2) + F(1)*(ar(2)/ar(1))*ust(1)

     ! Right BCS
     !pst(nx) = 0.0


     d(1) = ar(2)*uru/au(1,2)
     dloop: do i =2,nx-1
            d(i) = ar(2*i)*uru/au(i,2)
     end do dloop


     !if (time==1) then
     !Solving the U Momentum Equation
     ust(1) = uru*(au(1,3)*ust(2) + S(1) + (1-uru)*au(1,2)*ust(1)/uru)/au(1,2)
     uloop1: do i = 2,nx-1
            
            ust(i) = uru*(au(i,1)*ust(i-1) + au(i,3)*ust(i+1) + S(i) + &
              (1-uru)*au(i,2)*ust(i)/uru)/au(i,2)

     end do uloop1
     !end if

     

     ! P constants

     aloop: do i=1,nx-2

          a(i,1) = rho*d(i)*ar(2*i)
          a(i,3) = rho*d(i+1)*ar(2*i+2)
          a(i,2) = a(i,1) + a(i,3)
          Fstw(i) = rho*ust(i)*ar(2*i)
          Fste(i) = rho*ust(i+1)*ar(2*i+2)
          b(i) = Fstw(i) - Fste(i)

     end do aloop

     ! Assembling Pressure Correction Matrix (Tri-diagonal)
     pcr(1) = 0.0; pcr(nx) = 0.0
     !pmat = 0
     !pmat(1,1) = a(1,2)
     !pmat(2,1) = -a(1,3)
     !pmat(nx-3,nx-2) = -a(nx-2,1)
     !pmat(nx-2,nx-2) = a(nx-2,2)
     !prloop1: do j=2,nx-3
     !
     !       pmat(j-1,j) = -a(j,1)
     !       pmat(j,j) = a(j,2)
     !       pmat(j+1,j) = -a(j,3)
     !
     !end do prloop1

     !do i =1,n
     ! do j = 1, i
     !       pmat(i,j) = pmat(j,i)
     ! end do
     !end do
     !c=0


    !do i = 1,n
    !  diag(i) = pmat(i,i)
    !  updiag(i) = pmat(i+1,i)
    !end do
    !lowdiag = updiag
    !diag(n) = pmat(n,n)
 
     diag = 0
     updiag = 0
     lowdiag = 0
     diag(n) = a(nx-2,2)
     diag(1:n-1) = a(1:n-1,2)
     updiag(1:n-1) = -a(1:n-1,3)
     !prloop1: do j=1,n-1
     
     !       diag(j) = a(j,2)
     !       updiag(j) = -a(j,3)
     
     !end do prloop1
     lowdiag = updiag



    nrhs = 1
    lda = n
    ldb = n
    x = b
    call dgtsv(n,nrhs,lowdiag,diag,updiag,x,ldb,info)

    pcr = 0
    pcr(2:n+1) = x(1:n)
    !do i = 1,n
    !    pcr(i+1) = x(i)
    !enddo

     

     ! Pressure Correction
     p = (1-urp)*pst + urp*(pst + pcr)
     ! Velocity Correction
     u(1:nx-1) = (1-uru)*ust(1:nx-1) + uru*(ust(1:nx-1) + d(1:nx-1)*(pcr(1:nx-1)-pcr(2:nx)))
     !uloop2: do i=1,nx-1

     !       u(i) = (1-ur)*ust(i) + ur*(ust(i) + d(i)*(pcr(i)-pcr(i+1)))

     !end do uloop2

     ! Update inlet Pressure
     !p(1) = po - 0.5*rho*((u(1)*ar(2)/ar(1))**2)

     ! Residuals
     p_res(1:nx-2) = F(2:nx-1) - F(1:nx-2)
     m_res(2:nx-2) = au(2:nx-2,2)*u(2:nx-2) - au(2:nx-2,1)*ust(1:nx-3) - &
     au(2:nx-2,3)*ust(3:nx-1) - S(2:nx-2)
     !resloop: do i=1,nx-2
      
     !       p_res(i) = F(i+1) - F(i)
     !       m_res(i) = au(i,2)*u(i) - au(i,1)*ust(i-1) - au(i,3)*ust(i+1) - S(i)

     !enddo resloop

     ust = u ! Allocating u* = u
     pst = p ! Allocating p* = p


     print*,time, F(4), m_res(nx/2), p_res(nx/2)
     !print*, p

  time = time + 1

end do timeloop

! output data into a file 
  open(1, file = 'data1.dat', status = 'old')  
  do i=1,nx-1
  write(1,*) u(i)
  end do
!deallocate (pmat)
!print*, Fste, Fstw

end program oneDnozzle
