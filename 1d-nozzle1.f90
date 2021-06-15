    !1-D nozzle inviscid
program oneDnozzle
implicit none

       double precision, parameter :: m_dot=1.0
       double precision, parameter :: rho=1.0
       integer, parameter :: nx=1000, n = nx-2
       double precision, parameter :: L=2.00
       double precision, parameter :: dx=L/(nx-1)
       double precision ar(2*nx), x(nx), au(nx-1,3), a(nx-2,3), u(nx-1), Di(nx), F(nx), &
        S(nx-1), d(nx-1), po, p(nx), pcr(nx), ust(nx-1), Fstw(nx-2), Fste(nx-2), Dst(nx-2), b(nx-2), &
        c(nx-2,nx-2), pst(nx), m_res(nx-1), p_res(nx-1), ur

       double precision pmat(nx-2,nx-2)
       double precision :: ua=2.0
       integer :: i, j, nt, time

       ! Under Relaxation Factor
       ur = 0.5
       time = 1

       
       
       
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


timeloop: do nt = 1,200

       ! U constants
       F(1) = rho*ust(1)*ar(1)
       
       if(time==1) then
                    F(nx) = 1.0
       else
                    F(nx) = rho*ar(2*nx-2)*ust(nx-1)
       end if

       fdsloop: do i = 2,nx-1

                    F(i) = rho*(ust(i)+ust(i-1))*ar(2*i-1)/2
                    Di(i) = 0 !0.0000181*nx/L  ! Tau/(x1-x2)
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
       pst(nx) = 0.0


       d(1) = ar(2)*ur/au(1,2)
       dloop: do i =2,nx-1
                    d(i) = ar(2*i)*ur/au(i,2)
       end do dloop


       !if (time==1) then
       !Solving the U Momentum Equation
       ust(1) = ur*(au(1,3)*ust(2) + S(1) + (1-ur)*au(1,2)*ust(1)/ur)/au(1,2)
       uloop1: do i = 2,nx-1
                    
                    ust(i) = ur*(au(i,1)*ust(i-1) + au(i,3)*ust(i+1) + S(i) + &
                        (1-ur)*au(i,2)*ust(i)/ur)/au(i,2)

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
       pmat = 0
       pmat(1,1) = a(1,2)
       pmat(2,1) = -a(1,3)
       pmat(nx-3,nx-2) = -a(nx-2,1)
       pmat(nx-2,nx-2) = a(nx-2,2)
       prloop1: do j=2,nx-3

                    pmat(j-1,j) = -a(j,1)
                    pmat(j,j) = a(j,2)
                    pmat(j+1,j) = -a(j,3)

       end do prloop1

       do i =1,nx-2
        do j = 1, i
                    pmat(i,j) = pmat(j,i)
        end do
       end do
       c=0
       call solve(pmat,c,n)

       !write (*,202)
       !do i = 1,n
       !write (*,201)  (c(i,j),j=1,n)
       !end do
       !200 format (' Computing Inverse mpmattrix ',/,/, &
       !' Mpmattrix pmat')
       !201 format (6f12.6)
       !202 format (/,' Inverse mpmattrix pmat^{-1}')
       
       ! Solve {p} = [A]^(-1) {b}
       pcr = 0
       do i = 1,n
        do j = 1,n
                    pcr(i+1) = pcr(i+1) + b(j)*c(j,i)
        end do
       end do

       ! Pressure Correction
       p = pst + ur*pcr
       ! Velocity Correction
       uloop2: do i=1,nx-1

                    u(i) = (1-ur)*ust(i) + ur*(ust(i) + d(i)*(pcr(i)-pcr(i+1)))

       end do uloop2

       ! Update inlet Pressure
       !p(1) = po - 0.5*rho*((u(1)*ar(2)/ar(1))**2)

       ! Residuals
       resloop: do i=1,nx-1
        
                    p_res(i) = F(i+1) - F(i)
                    m_res(i) = au(i,2)*u(i) - au(i,1)*ust(i-1) - au(i,3)*ust(i+1) - S(i)

       enddo resloop

       ust = u ! Allocating u* = u
       pst = p ! Allocating p* = p


       print*, F(4)

    ! output data into a file 
    open(1, file = 'data1.dat', status = 'old')
    do i=1,nx-1
    write(1,*) p(i)
    end do

    time = time + 1

end do timeloop



















       contains
       subroutine solve(pmat,c,n)
        implicit none 
integer n
double precision pmat(n,n), c(n,n)
double precision L(n,n), U(n,n), r(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and r
! Fortran 90/95 aloows such operations on matrices
L(:,:)=0.0
U(:,:)=0.0
r(:)=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=pmat(i,k)/pmat(k,k)
      L(i,k) = coeff
      do j=k+1,n
         pmat(i,j) = pmat(i,j)-coeff*pmat(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of a
do j=1,n
  do i=1,j
    U(i,j) = pmat(i,j)
  end do
end do

! Step 3: compute columns of the solve matrix C
do k=1,n
  r(k)=1.0
  d(1) = r(1)
! Step 3pmat: Solve Ld=r using the forward surstitution
  do i=2,n
    d(i)=r(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  r(k)=0.0
end do
       end subroutine solve


end program oneDnozzle
