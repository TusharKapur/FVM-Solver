!1-D nozzle inviscid
program twoDpipe
implicit none

! EWNS
! nx: row (vert), ny: column (hor)

double precision, parameter :: m_dot=1.0
double precision, parameter :: tau = 0.0000181
integer, parameter :: nx=5, ny=10, nu = (nx)*(ny-1), nv = (nx-1)*(ny-2), np = (nx)*(ny-2)
double precision, parameter :: L=2.0, H = 2.0
double precision, parameter :: dx = L/(nx-1), dy =  H/(ny-1)
double precision arx(nx,2*ny), ary(2*nx+1,ny), au(nx,ny-1,5), av(nx-1,ny-2,5), &
u(nx,ny-1), v(nx+1,ny), ust(nx,ny-1), vst(nx+1,ny), Fxu(nx,ny), Fyu(nx-1,ny-1), &
Fxv(nx-1,ny-1), Fyv(nx,ny), ustar(nx,ny+1), vstar(nx+1,ny), &
Dix(nx,ny-1), Diy(nx-1,ny), Sx(nx,ny-1), Sy(nx-1,ny-2), du(nx,ny-1), dv(nx-1,ny-2), &
po, p(nx,ny), pcr(nx,ny), a(nx,ny-2,5), &
Fstx(nx,ny), Fsty(nx,ny), Fsts(nx,ny), bx(nx,ny), by(ny,ny), &
pst(nx,ny), m_res(nx-1,ny-1), p_res(nx-1,ny-1), uru, urp, ua, &
xp(np), bd(nx,ny-2), bp(np), xu(nu), bu(nu), xv(nv), bv(nv), &
abu(7,nu), &
umat(nu,nu), vmat(nv,nv), pmat(np,np), c(np,np), rho(nx,ny)
integer :: i, j, k, nt, time, nrhs, info, lda, ldb, n, iterp
integer, dimension(nu) :: ipivu
integer, dimension(nv) :: ipivv
integer, dimension(np) :: ipivp
rho(:,:) = 1.0
uru = 0.9; urp = 0.8
time = 1
arx = dy*1.0
ary = dx*1.0
! U, V Initial conditions
ua = 2.0
ust = ua
vst = 0.0

! P Initial conditions
po = 10.0
pst(:,1)=po
!psetloop: do i = 2,nx
!	pst(:,i) = pst(:,1) - (i-1)*(10-7.5)*dx/0.5
!end do psetloop
pst(:,2)=9.0; pst(:,ny-1)=1.0
!pst(2,:)=pst(1,:)-0.1
!pst(4,:)=pst(5,:)-0.1
!pst(3,:)=pst(2,:)-0.1
!pst = po

u = ust
p = pst

iterloop: do nt = 1,1
vst(1,:)=0.0; vst(nx+1,:)=0.0; vst(:,1) = 0.0
!vst=0.0
! U, V Constants

Fxu(1:nx,2:ny-1) = ((rho(1:nx,))*ust(1:nx,1:ny-2)/2+ust(1:nx,2:ny-1)/2)*arx(1,1)/2
Fxu(1:nx,ny) = rho*arx(1,1)*ust(1:nx,ny-1)
Fxu(1:nx,1) = rho*arx(1,1)*ust(1:nx,1)

Fyu(1:nx-1,1:ny-1) = rho*(vst(2:nx,1:ny-1)+vst(2:nx,2:ny))*ary(1,1)/2

Dix(:,:) = tau*arx(1,1)*nx/L  ! Tau.A/(x1-x2)
Diy(:,:) = tau*ary(1,1)*nx/L

Sx(1:nx,1:ny-1) = (pst(1:nx,1:ny-1)-pst(1:nx,2:ny))*arx(1,1)
Sx(1,1:ny-1) = Sx(1,1:ny-1) - tau*ary(1,1)*2/dx
Sx(nx,1:ny-1) = Sx(nx,1:ny-1) - tau*ary(1,1)*2/dx


au(:,:,:) = 0.0
au(1:nx,2:ny-2,2) = Dix(1:nx,2:ny-2) + max(-Fxu(1:nx,3:ny-1),0.0)
au(1:nx,2:ny-1,3) = Dix(1:nx,2:ny-1) + max(Fxu(1:nx,2:ny-1),0.0)
au(1:nx-1,1:ny-1,4) = Diy(1:nx-1,1:ny-1) + max(-Fyu(1:nx-1,1:ny-1),0.0)
au(2:nx,1:ny-1,5) = Diy(1:nx-1,1:ny-1) + max(Fyu(1:nx-1,1:ny-1),0.0)

au(2:nx-1,1:ny-1,1) = (au(2:nx-1,1:ny-1,2) + au(2:nx-1,1:ny-1,3) + au(2:nx-1,1:ny-1,4) + &
au(2:nx-1,1:ny-1,5) + Fxu(2:nx-1,2:ny) - Fxu(2:nx-1,1:ny-1) + Fyu(2:nx-1,1:ny-1) - &
Fyu(1:nx-2,1:ny-1))/uru

au(1,1:ny-1,1) = (au(1,1:ny-1,2) + au(1,1:ny-1,3) + au(1,1:ny-1,4) + &
au(1,1:ny-1,5) + Fxu(1,2:ny) - Fxu(1,1:ny-1))/uru! + Fyu(1,1:ny-1)/uru
au(nx,1:ny-1,1) = (au(nx,1:ny-1,2) + au(nx,1:ny-1,3) + au(nx,1:ny-1,4) + &
au(nx,1:ny-1,5) + Fxu(nx,2:ny) - Fxu(nx,1:ny-1))/uru! - Fyu(nx-1,1:ny-1)/uru

au(:,1,:)=0.0
au(1:nx,1,1) = (Fxu(1:nx,2) + 0.5*Fxu(1:nx,1))/uru


bu = 0.0

do i = 1,nx
	bu(((i-1)*(ny-1)+1):(i)*(ny-1)) = Sx(i,1:ny-1) + (1-uru)*au(i,1:ny-1,1)*ust(i,1:ny-1)/uru
	bu(((i-1)*(ny-1)+1)) = Sx(i,1) + Fxu(i,1)*ust(i,1) + (1-uru)*au(i,1,1)*ust(i,1)/uru ! S term addition
enddo

umat = 0.0
abu = 0.0

k=1
do i=1,nx
	do j=1,ny-1
		umat(k,k) = au(i,j,1)
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

!		abu(5,k) = au(i,j,1)
!		abu(4,k+1) = -au(i,j,2)
!		abu(3,k+2) = -au(i,j,4)
!		abu(6,k) = -au(i,j,3)
!		abu(7,k) = -au(i,j,5)
		k=k+1
	enddo
enddo
!ust=ua
nrhs = 1

!call dgbsv(nu, 2, 2, nrhs, abu, 7, ipivu, bu, nu, info)
call dgesv(nu, nrhs, umat, nu, ipivu, bu, nu, info)
do i = 1,nx
	ust(i,1:ny-1) = bu(((i-1)*(ny-1)+1):(i)*(ny-1))
enddo

!print*, bu
vst(:,ny) = vst(:,ny-1); vst(1,:)=0.0; vst(nx+1,:)=0.0

Fxv(1:nx-1,1:ny-1) = rho*(ust(1:nx-1,1:ny-1)+ust(2:nx,1:ny-1))*arx(1,1)/2
Fyv(1:nx,1:ny) = rho*(vst(1:nx,1:ny)+vst(2:nx+1,1:ny))*ary(1,1)/2

Diy(:,:) = tau*ary(1,1)*nx/L

Sy(1:nx-1,1:ny-2) = (pst(1:nx-1,2:ny-1)-pst(2:nx,2:ny-1))*ary(1,1)

! v from 1 to nx+1; Set a from 2 to nx
av(:,:,:) = 0.0
av(1:nx-1,1:ny-3,2) = Dix(1:nx-1,2:ny-2) + max(-Fxv(1:nx-1,2:ny-2),0.0)
av(1:nx-1,2:ny-2,3) = Dix(1:nx-1,2:ny-2) + max(Fxv(1:nx-1,2:ny-2),0.0)
av(1:nx-2,1:ny-2,4) = Diy(1:nx-2,2:ny-1) + max(-Fyv(2:nx-1,2:ny-1),0.0)
av(2:nx-1,1:ny-2,5) = Diy(2:nx-1,2:ny-1) + max(Fyv(2:nx-1,2:ny-1),0.0)

av(1:nx-1,1:ny-2,1) = (av(1:nx-1,1:ny-2,2) + av(1:nx-1,1:ny-2,3) + av(1:nx-1,1:ny-2,4) + &
av(1:nx-1,1:ny-2,5) + Fxv(1:nx-1,2:ny-1) - Fxv(1:nx-1,1:ny-2) + Fyv(2:nx,2:ny-1) - &
Fyv(1:nx-1,2:ny-1))/uru

!av(1:nx-1,ny,1) = (av(1:nx-1,ny,2) + av(1:nx-1,ny,3) + av(1:nx-1,ny,4) + av(1:nx-1,ny,5) + &
!Fyv(2:nx,ny) - Fyv(1:nx-1,ny))/uru !- Fxv(1:nx-1,ny-1)/uru


bv = 0.0
do i=1,nx-1
	bv(((i-1)*(ny-2)+1):(i)*(ny-2)) = Sy(i,1:ny-2) + (1-uru)*av(i,1:ny-2,1)*vst(i+1,2:ny-1)/uru
enddo

k=1
vmat=0.0
do i=1,nx-1
	do j=1,ny-2
		vmat(k,k) = av(i,j,1)
		if(k>1) then
		vmat(k,k-1)=-av(i,j,3)
		end if
		if(k<nv) then
		vmat(k,k+1)=-av(i,j,2)
		end if
		if(k>ny-2) then
		vmat(k,k-ny+2)=-av(i,j,5)
		end if
		if(k<nv-ny+3) then
		vmat(k,k+ny-2)=-av(i,j,4)
		end if
		k=k+1
	enddo
enddo

info=0; nrhs=1
call dgesv(nv, nrhs, vmat, nv, ipivv, bv, nv, info)


do i=1,nx-1
	vst(i+1,2:ny-1) = bv(((i-1)*(ny-2)+1):(i)*(ny-2))
enddo
vst(:,ny) = vst(:,ny-1); vst(1,:)=0.0; vst(nx+1,:)=0.0
!vst=0.0
!vst(2:3,:)=0.1
!vst(3:5,:)=-0.1
!print*, vst

ploop: do iterp=1,1

du=0.0; dv=0.0; a(:,:,:)=0.0; bd=0.0; bp=0.0

du(1:nx,1:ny-1) = arx(1,1)*uru/au(1:nx,1:ny-1,1)
dv(1:nx-1,1:ny-2) = ary(1,1)*uru/av(1:nx-1,1:ny-2,1)

a(1:nx,1:ny-3,2) = rho*arx(1,1)*du(1:nx,2:ny-2)
a(1:nx,2:ny-2,3) = rho*arx(1,1)*du(1:nx,2:ny-2)
a(1:nx-1,1:ny-2,4) = rho*ary(1,1)*dv(1:nx-1,1:ny-2)
a(2:nx,1:ny-2,5) = rho*ary(1,1)*dv(1:nx-1,1:ny-2)
a(1:nx,1:ny-2,1) = a(1:nx,1:ny-2,2) + a(1:nx,1:ny-2,3) + a(1:nx,1:ny-2,4) + a(1:nx,1:ny-2,5)

!a(1:nx,1,1) = a(1:nx,1,1) + po*rho*arx(1,1)*du(1:nx,1)

bd(1:nx,1:ny-2) = rho*arx(1,1)*(ust(1:nx,1:ny-2)-ust(1:nx,2:ny-1)) + &
rho*ary(1,1)*(vst(1:nx,2:ny-1)-vst(2:nx+1,2:ny-1))

do i=1,nx
	bp(((i-1)*(ny-2)+1):(i)*(ny-2)) = bd(i,1:ny-2)
enddo

!print*, dv(1:nx-1,1:ny-1)
k=1; pmat=0.0
do i=1,nx
	do j=1,ny-2
		pmat(k,k) = a(i,j,1)
		if(k>1) then
		pmat(k,k-1)=-a(i,j,3)
		end if
		if(k<np) then
		pmat(k,k+1)=-a(i,j,2)
		end if
		if(k>ny-2) then
		pmat(k,k-ny+2)=-a(i,j,5)
		end if
		if(k<np-ny+3) then
		pmat(k,k+ny-2)=-a(i,j,4)
		end if
		k=k+1
	enddo
enddo

do i=1,nx
	bp(((i-1)*(ny-2)+1)) = -0.0
	pmat(((i-1)*(ny-2)+1),:)=0.0
	pmat(((i-1)*(ny-2)+1),((i-1)*(ny-2)+1))=1.0

	bp(((i)*(ny-2))) = -0.0
	pmat(((i)*(ny-2)),:)=0.0
	pmat(((i)*(ny-2)),((i)*(ny-2)))=1.0
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
!202 format (/,' Inverse mpmattrix pmat^{-1}')


info=0; nrhs=1
call dgesv(np, nrhs, pmat, np, ipivp, bp, np, info)

!bp(:)=bp(:)-bp(np)
pcr = 0.0
do i=1,nx
	pcr(i,2:ny-1) = bp(((i-1)*(ny-2)+1):(i)*(ny-2))
enddo

p = (1-urp)*pst + urp*(pcr+pst)
p(:,1)=po; p(:,ny)=0.0
p(1,:)=p(2,:); p(nx-1,:)=p(nx,:)
!p(:,1) = p(:,1) - 0.5*rho*(u(:,1)**2)

u(:,1:ny-1) = (1-uru)*ust(:,1:ny-1) + uru*(ust(:,1:ny-1) + &
du(:,1:ny-1)*(pcr(:,1:ny-1)-pcr(:,2:ny)))
v(2:nx,2:ny-1) = (1-uru)*vst(2:nx,2:ny-1) + uru*(vst(2:nx,2:ny-1) + &
dv(1:nx-1,1:ny-2)*(pcr(1:nx-1,2:ny-1)-pcr(2:nx,2:ny-1)))

print*, ust
enddo ploop
!print*, ust


ust=u; vst=v; pst=p
!print*, vst(1:nx,:)
!print*, ust

enddo iterloop



end program twoDpipe
