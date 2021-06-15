!1-D nozzle inviscid
program twoDpipe
implicit none

! EWNS
! nx: row (vert), ny: column (hor)

double precision, parameter :: m_dot=1.0
double precision, parameter :: rho=1.0, tau = 0.0000181
integer, parameter :: nx=5, ny=5, nu = (nx)*(ny-1), nv = (nx-1)*(ny-2), np = (nx)*(ny)
double precision, parameter :: L=2.0, H = 2.0
double precision, parameter :: dx = L/(nx-1), dy =  H/(ny-1)
double precision arx(nx,2*ny), ary(2*nx+1,ny), au(nx,ny-1,5), av(nx+1,ny,5), &
u(nx,ny+1), v(nx+1,ny), ust(nx,ny+1), vst(nx+1,ny), Fxu(nx,ny), Fyu(nx-1,ny-1), &
Fxv(nx,ny+2), Fyv(nx+2,ny), &
Dix(nx,ny), Diy(nx,ny+1), Sx(nx+1,ny+1), Sy(nx+1,ny+1), du(nx,ny+1), dv(nx+1,ny), &
po, p(nx,ny), pcr(nx,ny), a(nx,ny,5), &
Fstx(nx,ny), Fsty(nx,ny), Fsts(nx,ny), bx(nx,ny), by(ny,ny), &
pst(nx,ny), m_res(nx-1,ny-1), p_res(nx-1,ny-1), ur, ua, &
diagp(np), lowdiagp(np-1), updiagp(np-1), xp(np), bp(np), &
diagu(nu), lowdiagu(nu-1), updiagu(nu-1), lowdiagu2(nu-2), updiagu2(nu-2), xu(nu), bu(nu), &
diagv(nv), lowdiagv(nv-1), updiagv(nv-1), xv(nv), bv(nv), &
umat(nu,nu), c(nu, nu)
integer :: i, j, k, nt, time, nrhs, info, lda, ldb, n
integer, dimension(nu) :: ipiv
! Under Relaxation Factor
ur = 1.0
time = 1
arx = dy*1.0
ary = dx*1.0
! U, V Initial conditions
ua = 2.0
ust = ua
vst = 0.0
! P Initial conditions
po = 10.0

!psetloop: do i = 1,ny
!pst(:,i) = po - (i-1)*(10-7.5)*dx/0.5
!end do psetloop
pst = 0.0

u = ust
p = pst

timeloop: do nt = 1,1

! U, V Constants

Fxu(1:nx,ny) = rho*arx(1,1)*ust(1:nx,ny)
Fxu(1:nx,1:ny-1) = rho*(ust(1:nx,1:ny-1)+ust(1:nx,2:ny))*arx(1,1)/2
Fyu(1:nx-1,1:ny-1) = rho*(vst(2:nx,1:ny-1)+vst(2:nx,1:ny-1))*ary(1,1)/2

!Fxv(1:nx-1,:) = rho*(ust(1:nx-1,:)+ust(2:nx,:))*arx(1,1)/2
!Fyv(1:nx,:) = rho*(vst(1:nx,:)+vst(2:nx+1,:))*ary(1,1)/2

Dix(:,:) = tau*arx(1,1)*nx/L  ! Tau.A/(x1-x2)
Diy(:,:) = tau*ary(1,1)*nx/L
Sx(1:nx,1:ny-1) = (pst(:,1:ny-1)-pst(:,2:ny))*arx(1,1)
Sx(1,1:ny-1) = Sx(1,1:ny-1) - tau*ary(1,1)*2/dx
Sx(nx,1:ny-1) = Sx(nx,1:ny-1) - tau*ary(1,1)*2/dx
Sy(1:nx-1,1:ny) = (pst(1:nx-1,:)-pst(2:nx,:))*ary(1,1)
! IC: eg. au(:,ny+1,:) = 0.0 

au(1:nx,1:ny-2,2) = Dix(1:nx,2:ny-1) + max(-Fxu(1:nx,2:ny-1),0.0)
au(1:nx,2:ny-1,3) = Dix(1:nx,2:ny-1) + max(Fxu(1:nx,2:ny-1),0.0)
au(1:nx-1,1:ny-1,4) = Diy(1:nx-1,1:ny-1) + max(-Fyu(1:nx-1,1:ny-1),0.0)
au(2:nx,1:ny-1,5) = Diy(2:nx,1:ny-1) + max(Fyu(1:nx-1,1:ny-1),0.0)
au(2:nx-1,1:ny-1,1) = au(2:nx-1,1:ny-1,2) + au(2:nx-1,1:ny-1,3) + au(2:nx-1,1:ny-1,4) + &
au(2:nx-1,1:ny-1,5) + Fxu(2:nx-1,2:ny) - Fxu(2:nx-1,1:ny-1) + Fyu(2:nx-1,1:ny-1) - &
Fyu(1:nx-2,1:ny-1)
au(1,1:ny-1,1) = au(1,1:ny-1,2) + au(1,1:ny-1,3) + au(1,1:ny-1,4) + &
au(1,1:ny-1,5) + Fxu(1,2:ny) - Fxu(1,1:ny-1)! + Fyu(1,1:ny-1)
au(nx,1:ny-1,1) = au(nx,1:ny-1,2) + au(nx,1:ny-1,3) + au(nx,1:ny-1,4) + &
au(nx,1:ny-1,5) + Fxu(nx,2:ny) - Fxu(nx,1:ny-1)! - Fyu(nx-1,1:ny-1)
au(1:nx,1,1) = au(1:nx,1,2) + au(1:nx,1,3) + au(1:nx,1,4) + au(1:nx,1,5) + Fxu(1:nx,2) - Fxu(1:nx,1) &
+ Fxu(1:nx,1)

bu = 0.0
umat = 0.0

do i = 1,nx
	diagu(((i-1)*(ny-1)+1):(i)*(ny-1)) = au(i,1:ny-1,1)
	if(i>0) then
		lowdiagu(((i-1)*(ny-1)):(i)*(ny-1)) = -au(i,2:ny-1,3)
	end if
	if(i<ny) then
		updiagu(((i-1)*(ny-1)+1):(i)*(ny-1)-1) = -au(i,1:ny-2,2)
	end if
	if(i>2) then
		lowdiagu2(((i-1)*(ny-1)):(i)*(ny-1)) = -au(i,1:ny-1,5)
	end if
	if(i<ny-2) then
		updiagu2(((i-1)*(ny-1)+1):(i)*(ny-1)) = -au(i,1:ny-1,4)
	end if
	bu(((i-1)*(ny-1)+1):(i)*(ny-1)) = Sx(i,1:ny-1)
	bu(((i-1)*(ny-1)+1)) = Fxu(i,1)*ua ! S term addition
enddo

do i = 1,nu
	umat(i,i) = diagu(i)
enddo
do i = 2,nu
	umat(i-1,i) = updiagu(i-1)
	umat(i,i-1) = lowdiagu(i)
end do
do i = 3,nu
	umat(i-2,i) = updiagu2(i-2)
	umat(i,i-2) = lowdiagu2(i-2)
enddo

k=1
umat = 0.0
do i=1,nx
	do j=1,ny-1
		umat(k,k)=au(i,j,1)
		if(k>1) then
		umat(k,k-1)=-au(i,j,3)
		end if
		if(k<nu) then
		umat(k,k+1)=-au(i,j,2)
		end if
		if(k>ny-1) then
		umat(k,k-ny+1)=-au(i,j,5)
		end if
		if(k<nu-ny+2) then
		umat(k,k+ny-1)=-au(i,j,4)
		end if
		k=k+1
	enddo
enddo
!ust=ua
nrhs = 1
n = nu

call dgesv(nu, nrhs, umat, nu, ipiv, bu, nu, info)

!call solve(umat,c,n)
!xu = 0
!       do i = 1,nu
!        do j = 1,nu
!                    xu(i) = xu(i) + bu(j)*c(j,i)
!        end do
!       end do
!bu = xu
do i = 1,nx
	ust(i,2:ny) = bu(((i-1)*(ny-1)+1):(i)*(ny-1))
enddo
print*, umat(1,:)

!write (*,202)
!do i = 1,nu
!write (*,201)  (c(i,j),j=1,n)
!end do
!200 format (' Computing Inverse mumattrix ',/,/, &
!' Mumattrix umat')
!201 format (6f12.6)
!02 format (/,' Inverse mumattrix umat^{-1}')

!diagu = 0.0
!lowdiagu = 0.0
!updiagu = 0.0
!bu = 0.0
!do i = 1,ny-1
!	diagu(((i-1)*(nx)+1):(i)*(nx)) = au(1:nx,i,1)
!	if(i>1) then
!		lowdiagu(((i-1)*(nx)):(i)*(nx)) = -au(1:nx,i,3)
!	end if
!	if(i<ny-1) then
!		updiagu(((i-1)*(nx)+1):(i)*(nx)) = -au(1:nx,i,2)
!	end if

!	bu(((i-1)*(nx)+1):(i)*(nx)) = Sx(1:nx,i)
!	bu(((i-1)*(nx)+2):(i)*(nx)) = bu(((i-1)*(nx)+2):(i)*(nx)) + au(2:nx,i,4)*ust(2:nx,i)
!	bu(((i-1)*(nx)+1):((i)*(nx)-1)) = bu(((i-1)*(nx)+1):((i)*(nx)-1)) + au(1:nx-1,i,5)*ust(1:nx-1,i)
!enddo
!call dgtsv(n,1,lowdiagu,diagu,updiagu,bu,n,info)

!do i = 1,ny-1
!	ust(1:nx,i+1) = bu(((i-1)*(nx)+1):(i)*(nx))
!enddo


!diagu = 0.0
!lowdiagu = 0.0
!updiagu = 0.0
!bu = 0.0
!do i = 1,ny-1
!	diagu(((i-1)*(nx)+1):(i)*(nx)) = au(1:nx,i,1)
!	if(i>1) then
!		lowdiagu(((i-1)*(nx)):(i)*(nx)) = -au(1:nx,i,4)
!		bu(((i-1)*(nx)+1):(i)*(nx)) = au(1:nx,i,3)*ust(1:nx,i-1)
!	end if
!	if(i<ny-1) then
!		updiagu(((i-1)*(nx)+1):(i)*(nx)) = -au(1:nx,i,5)
!		bu(((i-1)*(nx)+1):(i)*(nx)) = au(1:nx,i,2)*ust(1:nx,i+1)
!	end if

!	bu(((i-1)*(nx)+1):(i)*(nx)) = bu(((i-1)*(nx)+1):(i)*(nx)) + Sx(1:nx,i)

!enddo
!call dgtsv(n,1,lowdiagu,diagu,updiagu,bu,n,info)

!do i = 1,ny-1
!	ust(1:nx,i+1) = bu(((i-1)*(nx)+1):(i)*(nx))
!enddo

!av(2:nx,:,2) = Dix(1:nx-1,:) + max(Fxv(1:nx-1,:),0.0)
!av(1:nx-1,:,3) = Dix(2:nx,:) + max(-Fxv(2:nx,:),0.0)
!av(:,2:ny,4) = Diy(:,1:ny-1) + max(Fyv(:,1:ny-1),0.0)
!av(:,1:ny-1,5) = Diy(:,2:ny) + max(-Fyv(:,2:ny),0.0)
!av(2:nx-1,2:ny-1,1) = av(2:nx-1,2:ny-1,2) + av(2:nx-1,2:ny-1,3) + av(2:nx-1,2:ny-1,4) + &
!av(2:nx-1,2:ny-1,5) + Fxv(2:nx-1,2:ny-1) - Fxv(2:nx-1,1:ny-2) + Fyv(2:nx-1,2:ny-1) - &
!Fyv(1:nx-2,2:ny-1)






enddo timeloop


!write (*,202)
!do i = 1,n
!write (*,201)  (pst(i,j),j=1,n)
!end do
!200 format (' Computing Inverse mumattrix ',/,/, &
!' Mumattrix umat')
!201 format (6f12.6)
!202 format (/,' Inverse mumattrix umat^{-1}')


end program twoDpipe
