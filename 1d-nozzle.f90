	!1-D nozzle inviscid
program oneDnozzle
implicit none

       double precision, parameter :: m_dot=1.0
       double precision, parameter :: rho=1.0
       integer, parameter :: nx=5, n = nx-2
	   double precision, parameter :: L=2.00
	   double precision, parameter :: dx=L/(nx-1)
	   double precision ar(2*nx), x(nx), au(nx-1,nx-2), a(nx-2,nx-2), u(nx-1), Di(nx), F(nx), &
	    S(nx-1), d(nx-1), po, p(nx), pcr(nx), ust(nx-1), Fstw(nx-2), Fste(nx-2), Dst(nx-2), b(nx-2), &
	    c(nx-2,nx-2), pst(nx), m_res(nx-1), p_res(nx-1), ur

	   double precision pmat(nx-2,nx-2)
	   double precision :: ua=2.0
	   integer :: i, j, nt

	   ! Under Relaxation Factor
	   ur = 0.2


	   ! Dimension of the Pressure Correction Matrix
	   
	   
	   allocationloop: do i = 1, nx

                    x(i) = (i-1)*dx
                    ar(2*i) = 0.55 - i*(0.45-0.35)*dx/0.5
                    ar(2*i-1) = 0.6 - i*(0.45-0.35)*dx/0.5

	   end do allocationloop
	   
	   
	   F(1) = rho*u(1)*ar(1)
	   F(nx) = 1.0

	   ! U Initial conditions
	   usetloop: do i = 1,nx-1

	   				u(i) = m_dot/(rho*ar(2*i))
	   				
	   end do usetloop

	   ! P Initial conditions
	   po = 10.0
	   p(1) = po

	   psetloop: do i = 2,nx

	   				p(i) = p(i) - (i-1)*(10-7.5)*dx/0.5

	   end do psetloop

timeloop: do nt = 1,1

	   ! U constants
	   fdsloop: do i = 2,nx-1

	   				F(i) = rho*(u(i)+u(i-1))*ar(2*i-1)/2
	   				Di(i) = 0
	   				S(i) = (p(i)-p(i+1))*ar(2*i)

	   end do fdsloop

	   auloop: do i = 2,nx-1
	   				au(i,1) = Di(i) + max(F(i),0.0)
	   				au(i,3) = Di(i+1) + max(-F(i+1),0.0)
	   				au(i,2) = au(i,1) + au(i,3) + F(i+1) - F(i)

	   end do auloop

	   ! Left BCS
	   ua = u(1)*ar(2)/ar(1)
	   F(1) = rho*ua*ar(1)
	   au(1,:) = 0.0
	   au(1,2) = F(2) + F(1)*((ar(2)/ar(1))**2)/2
	   S(1) = (p(1)-p(2))*ar(2) + F(1)*(ar(2)/ar(1))*u(1)

	   ! Right BCS
	   p(nx) = 0.0

	   !Solving U
	   d(1) = ar(2)*ur/au(1,2)
	   u(1) = ur*(au(1,3)*u(2) + S(1) - (1-ur)*au(1,2)*u(1)/ur)/au(1,2)
	   uloop1: do i = 2,nx-1
	   				
	   				d(i) = ar(2*i)*ur/au(i,2)
	   				u(i) = ur*(au(i,1)*u(i-1) + au(i,3)*u(i+1) + S(i) - (1-ur)*au(i,2)*u(i)/ur)/au(i,2)

	   end do uloop1

	   ust(:) = u(:) ! Allocating u* = u
	   pst(:) = p(:) ! Allocating p* = p

	   ! P constants

	   aloop: do i=1,nx-2

					a(i,1) = rho*d(i)*ar(2*i)
					a(i,3) = rho*d(i+1)*ar(2*i+2)
					a(i,2) = rho*a(i,1) + a(i,3)
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
	   ploop1: do j=2,nx-3

	   				pmat(1,j) = -a(j,1)
	   				pmat(2,j) = a(j,2)
	   				pmat(3,j) = -a(j,3)

	   end do ploop1


	   call solve(pmat,c,n)
	   
	   ! Solve [A]{p}={b}
	   do i = 1,n
	   	do j = 1,n
	   				pcr(i+1) = pcr(i+1) + b(j)*c(j,i)
	   	end do
	   end do

	   ! Pressure Correction
	   p = (1-ur)*p + ur*(pst + pcr)
	   ! Velocity Correction
	   uloop2: do i=1,nx-2

	   				u(i) = (1-ur)*u(i) + ur*(ust(i) + d(i)*(pcr(i)-pcr(i+1)))

	   end do uloop2
	   u(nx-1) = (1-ur)*u(nx-1) + ur*(ust(nx-1) + d(nx-1)*(pcr(nx-1)-0.0))

	   ! Update inlet Pressure
	   p(1) = (1-ur)*p(1) + ur*(po - 0.5*rho*((u(1)*ar(2)/ar(1))**2))

	   ! Residuals
	   resloop: do i=1,nx-1
	   	
	   				m_res(i) = rho*u(i)*ar(2*i)
	   				p_res(i) = au(i,2)*u(i) - au(i,1)*u(i-1) - au(i,3)*u(i+1) - S(i)

	   enddo resloop


	   print*, pcr

	! output data into a file 
	open(1, file = 'data1.dat', status = 'old')  
	do i=1,nx-1
	write(1,*) p(i)
	end do

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
