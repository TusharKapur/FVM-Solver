!1-D nozzle inviscid   #############  U INLET!  #######################33
program cavity
implicit none

! EWNS
! nx: row (vert), ny: column (hor)

double precision, parameter :: m_dot=1.0
double precision, parameter :: rho=1.0, tau = 0.001794
integer, parameter :: nx=10, ny=10, nu = (nx)*(ny-1), nv = (nx-1)*(ny), np = (nx)*(ny)
double precision, parameter :: L=2.0, H = 2.0
double precision, parameter :: dx = H/(ny-1), dy = L/(nx-1)
double precision arx(nx,2*ny), ary(2*nx+1,ny), au(nx,ny-1,5), av(nx-1,ny,5), &
u(nx,ny+1), v(nx+1,ny), ust(nx,ny+1), vst(nx+1,ny), Fxu(nx,ny), Fyu(nx+1,ny-1), &
Dix, Diy, Fxv(nx-1,ny+1), Fyv(nx,ny), ustar(nx,ny+1), vstar(nx+1,ny), &
Sx(nx,ny-1), Sy(nx-1,ny), du(nx,ny-1), dv(nx-1,ny), &
po, p(nx,ny), pcr(nx,ny), a(nx,ny,5), &
pst(nx,ny), m_res(nx-1,ny-1), p_res(nx-1,ny-1), uru, urp, ua, &
xp(np), bd(nx,ny), bp(np), xu(nu), bu(nu), xv(nv), bv(nv), &
umat(nu,nu), vmat(nv,nv), pmat(np,np), c(np,np)
integer :: i, j, k, nt, time, nrhs, info, lda, ldb, n, iterp
integer, dimension(nu) :: ipivu
integer, dimension(nv) :: ipivv
integer, dimension(np) :: ipivp
uru = 0.7; urp = 0.5
time = 1
arx = dx*1.0
ary = dy*1.0
! U, V Initial conditions
ua = 2.0
ust = 0.0; vst = 0.0

! P Initial conditions
po = 0.0; pst=0.0

u = ust
p = pst

iterloop: do nt = 1,1000

vst(1,:)=0.0; vst(nx+1,:)=0.0; ust(:,1)=0.0; ust(:,ny+1)=0.0

Fxu(1:nx,1:ny) = rho*(ust(1:nx,1:ny)+ust(1:nx,2:ny+1))*arx(1,1)/2
Fyu(2:nx,1:ny-1) = rho*(vst(2:nx,1:ny-1)+vst(2:nx,2:ny))*ary(1,1)/2

Dix = tau*arx(1,1)*nx/L  ! Tau.A/(x1-x2)
Diy = tau*ary(1,1)*nx/L

Sx(1:nx,1:ny-1) = (pst(1:nx,1:ny-1)-pst(1:nx,2:ny))*arx(1,1)

au(:,:,:) = 0.0
au(1:nx,1:ny-2,2) = Dix + max(-Fxu(1:nx,2:ny-1),0.0)
au(1:nx,2:ny-1,3) = Dix + max(Fxu(1:nx,2:ny-1),0.0)
au(1:nx-1,1:ny-1,4) = Diy + max(-Fyu(2:nx,1:ny-1),0.0)
au(2:nx,1:ny-1,5) = Diy + max(Fyu(2:nx,1:ny-1),0.0)

au(1:nx,1:ny-1,1) = au(1:nx,1:ny-1,2) + au(1:nx,1:ny-1,3) + au(1:nx,1:ny-1,4) + &
au(1:nx,1:ny-1,5) + Fxu(1:nx,2:ny) - Fxu(1:nx,1:ny-1) + Fyu(2:nx+1,1:ny-1) - &
Fyu(1:nx,1:ny-1)

au(1,:,1)=au(1,:,1)+tau*ary(1,1)*2/dx; au(nx,:,1)=au(nx,:,1)+tau*ary(1,1)*2/dx

bu = 0.0

do i = 1,nx
	bu(((i-1)*(ny-1)+1):(i)*(ny-1)) = Sx(i,1:ny-1) + (1-uru)*au(i,1:ny-1,1)*ust(i,2:ny)/uru
!	bu(i*(ny-1)) = bu(i*(ny-1)) + ua*tau*ary(1,1)*2/dx
enddo
bu(1:ny) = bu(1:ny) + ua*tau*ary(1,1)*2/dx

umat = 0.0

k=1
do i=1,nx
	do j=1,ny-1
		umat(k,k) = au(i,j,1)/uru
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

nrhs = 1
call dgesv(nu, nrhs, umat, nu, ipivu, bu, nu, info)
do i = 1,nx
	ust(i,2:ny) = bu(((i-1)*(ny-1)+1):(i)*(ny-1))
enddo


Fxv(1:nx-1,2:ny) = rho*(ust(1:nx-1,2:ny)+ust(2:nx,2:ny))*arx(1,1)/2
Fyv(1:nx,1:ny) = rho*(vst(1:nx,1:ny)+vst(2:nx+1,1:ny))*ary(1,1)/2

Sy(1:nx-1,1:ny) = (pst(1:nx-1,1:ny)-pst(2:nx,1:ny))*ary(1,1)

! v from 1 to nx+1; Set a from 2 to nx
av(:,:,:) = 0.0
av(1:nx-1,1:ny-1,2) = Dix + max(-Fxv(1:nx-1,2:ny),0.0)
av(1:nx-1,2:ny,3) = Dix + max(Fxv(1:nx-1,2:ny),0.0)
av(1:nx-2,1:ny,4) = Diy + max(-Fyv(2:nx-1,1:ny),0.0)
av(2:nx-1,1:ny,5) = Diy + max(Fyv(2:nx-1,1:ny),0.0)

av(1:nx-1,1:ny,1) = av(1:nx-1,1:ny,2) + av(1:nx-1,1:ny,3) + av(1:nx-1,1:ny,4) + &
av(1:nx-1,1:ny,5) + Fxv(1:nx-1,2:ny+1) - Fxv(1:nx-1,1:ny) + Fyv(2:nx,1:ny) - &
Fyv(1:nx-1,1:ny)

av(:,1,1)=av(:,1,1)+tau*arx(1,1)*2/dy; av(:,ny-1,1)=av(:,ny-1,1)+tau*arx(1,1)*2/dy


bv = 0.0
do i=1,nx-1
	bv(((i-1)*(ny)+1):(i)*(ny)) = Sy(i,1:ny) + (1-uru)*av(i,1:ny,1)*vst(i+1,1:ny)/uru
enddo

k=1
vmat=0.0
do i=1,nx-1
	do j=1,ny
		vmat(k,k) = av(i,j,1)/uru
		if(k>1) then
		vmat(k,k-1)=-av(i,j,3)
		end if
		if(k<nv) then
		vmat(k,k+1)=-av(i,j,2)
		end if
		if(k>ny) then
		vmat(k,k-ny)=-av(i,j,5)
		end if
		if(k<nv-ny+1) then
		vmat(k,k+ny)=-av(i,j,4)
		end if
		k=k+1
	enddo
enddo

info=0; nrhs=1
call dgesv(nv, nrhs, vmat, nv, ipivv, bv, nv, info)


do i=1,nx-1
	vst(i+1,1:ny) = bv(((i-1)*(ny)+1):(i)*(ny))
enddo


ploop: do iterp=1,1

du=0.0; dv=0.0; a(:,:,:)=0.0; bd=0.0; bp=0.0

du(1:nx,1:ny-1) = arx(1,1)*uru/au(1:nx,1:ny-1,1)
dv(1:nx-1,1:ny) = ary(1,1)*uru/av(1:nx-1,1:ny,1)

a(1:nx,1:ny-1,2) = rho*arx(1,1)*du(1:nx,1:ny-1)
a(1:nx,2:ny,3) = rho*arx(1,1)*du(1:nx,1:ny-1)
a(1:nx-1,1:ny,4) = rho*ary(1,1)*dv(1:nx-1,1:ny)
a(2:nx,1:ny,5) = rho*ary(1,1)*dv(1:nx-1,1:ny)
a(1:nx,1:ny,1) = a(1:nx,1:ny,2) + a(1:nx,1:ny,3) + a(1:nx,1:ny,4) + a(1:nx,1:ny,5)

a(1,1,1) = a(1,1,1) + (10**5)

bd(1:nx,1:ny) = rho*arx(1,1)*(ust(1:nx,1:ny)-ust(1:nx,2:ny+1)) + &
rho*ary(1,1)*(vst(1:nx,1:ny)-vst(2:nx+1,1:ny))

do i=1,nx
	bp(((i-1)*(ny)+1):(i)*(ny)) = bd(i,1:ny)
enddo
bp(1) = bp(1) + (10**5)*po

k=1; pmat=0.0
do i=1,nx
	do j=1,ny
		pmat(k,k) = a(i,j,1)
		if(k>1) then
		pmat(k,k-1)=-a(i,j,3)
		end if
		if(k<np) then
		pmat(k,k+1)=-a(i,j,2)
		end if
		if(k>ny) then
		pmat(k,k-ny)=-a(i,j,5)
		end if
		if(k<np-ny+1) then
		pmat(k,k+ny)=-a(i,j,4)
		end if
		k=k+1
	enddo
enddo


!print*, bd(1:nx,1:ny-2)
!do i=1,np-ny+2
!	print*, pmat(i+ny-2,i)
!enddo

!write (*,202)
!do i = 1,np
!write (*,201)  (pmat(i,j),j=1,np)
!end do
!200 format (' Computing Inverse mpmattrix ',/,/, &
!' Mpmattrix pmat')
!201 format (6f12.6)
!202 format (/,'pmat')

!print*, bp

info=0; nrhs=1
call dgesv(np, nrhs, pmat, np, ipivp, bp, np, info)

!bp(:)=bp(:)-bp(np)
pcr = 0.0
do i=1,nx
	pcr(i,1:ny) = bp(((i-1)*(ny)+1):(i)*(ny))
enddo
!pcr(1,1) = 0.0

p = (1-urp)*pst + urp*(pcr+pst)

u(1:nx,2:ny) = (1-uru)*ust(1:nx,2:ny) + uru*(ust(1:nx,2:ny) + &
du(1:nx,1:ny-1)*(pcr(1:nx,1:ny-1)-pcr(1:nx,2:ny)))

v(2:nx,1:ny) = (1-uru)*vst(2:nx,1:ny) + uru*(vst(2:nx,1:ny) + &
dv(1:nx-1,1:ny)*(pcr(1:nx-1,1:ny)-pcr(2:nx,1:ny)))

!print*, u
enddo ploop
!print*, p

print*, ust
!print*, vst(1:nx,:)
ust=u; vst=v; pst=p

!print*, dv

enddo iterloop

open(1, file = 'data-u.dat', status = 'new')  
do i=1,nx+1
	write(1,*) vst(i,1:ny)   
end do 

open(1, file = 'data-v.dat', status = 'new')  
do i=1,ny+1
	write(1,*) ust(i,1:nx)
end do 


end program cavity
