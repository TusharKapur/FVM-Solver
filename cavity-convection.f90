!1-D nozzle inviscid   #############  U INLET!  #######################33
program cavity
implicit none

! EWNS
! nx: row (vert), ny: column (hor)

double precision, parameter :: m_dot=1.0
double precision, parameter :: tau = 0.00002088
integer, parameter :: nx=10, ny=10, nu = (nx)*(ny-1), nv = (nx-1)*(ny), np = (nx)*(ny)
double precision, parameter :: L=2.0, H = 2.0
double precision, parameter :: dx = H/(ny-1), dy = L/(nx-1)
double precision arx(nx,2*ny), ary(2*nx+1,ny), au(nx,ny-1,5), av(nx-1,ny,5), &
u(nx,ny+1), v(nx+1,ny), ust(nx,ny+1), vst(nx+1,ny), Fxu(nx,ny), Fyu(nx+1,ny-1), &
Dix, Diy, Fxv(nx-1,ny+1), Fyv(nx,ny), ustar(nx,ny+1), vstar(nx+1,ny), &
Sx(nx,ny-1), Sy(nx-1,ny), du(nx,ny-1), dv(nx-1,ny), &
po, p_abs, p(nx,ny), pcr(nx,ny), a(nx,ny,5), &
pst(nx,ny), m_res, p_res, uru, urp, ua, &
bd(nx,ny), bp(nx,ny), bu(nx,ny-1), bv(nx-1,ny), &
rho(nx+2,ny+2), kt, cp, at(nx,ny,5), T(nx,ny), T1, T2, T0, urt, bt(nx,ny), Pr, R
integer :: i, j, k, nt, time, nrhs, info, n, iteru, iterv, iterp, itert
integer, dimension(nu) :: ipivu
integer, dimension(nv) :: ipivv
integer, dimension(np) :: ipivp

double precision, dimension (:), allocatable :: diag, updiag, lowdiag, x
! Note: AIR properties, tdma for t, gas eq.
! 70 - 80 - 90K
T1=0.000001; T2=-0.000001; T0=0; kt=0.03023; cp = 1010; Pr = 0.701; R = 8314/28.97
rho = 1.0; T = T0
iteru = 10; iterv = 10; iterp = 10; itert = 10
uru = 0.5; urp = 0.5; urt = 1.0
time = 1
arx = dx*1.0
ary = dy*1.0
! U, V Initial conditions
ua = 0.02
ust = 0.0; vst = 0.0

! P Initial conditions
po = 0.0; pst=0.0; p_abs=101325

u = ust
p = pst

iterloop: do nt = 1,10000
vst(1,:)=0.0; vst(nx+1,:)=0.0; ust(:,1)=0.0; ust(:,ny+1)=0.0
rho(1,:)=rho(2,:); rho(nx+2,:)=rho(nx+1,:)
rho(:,1)=rho(:,2); rho(:,ny+2)=rho(:,ny+1)

Fxu(1:nx,1:ny) = ((rho(2:nx+1,1:ny)+rho(2:nx+1,2:ny+1))*ust(1:nx,1:ny) + &
	(rho(2:nx+1,2:ny+1)+rho(2:nx+1,3:ny+2))*ust(1:nx,2:ny+1))*arx(1,1)/4
Fyu(2:nx,1:ny-1) = ((rho(2:nx,2:ny)+rho(3:nx+1,2:ny))*vst(2:nx,1:ny-1) + &
	(rho(2:nx,3:ny+1)+rho(3:nx+1,3:ny+1))*vst(2:nx,2:ny))*ary(1,1)/4

Dix = tau*arx(1,1)*ny/H  ! Tau.A/(x1-x2)
Diy = tau*ary(1,1)*nx/L

Sx(1:nx,1:ny-1) = (pst(1:nx,1:ny-1)-pst(1:nx,2:ny))*arx(1,1)

au(:,:,:) = 0.0
au(1:nx,1:ny-1,2) = Dix + max(-Fxu(1:nx,2:ny),0.0)
au(1:nx,1:ny-1,3) = Dix + max(Fxu(1:nx,1:ny-1),0.0)
au(1:nx-1,1:ny-1,4) = Diy + max(-Fyu(2:nx,1:ny-1),0.0)
au(2:nx,1:ny-1,5) = Diy + max(Fyu(2:nx,1:ny-1),0.0)

au(1:nx,1:ny-1,1) = au(1:nx,1:ny-1,2) + au(1:nx,1:ny-1,3) + au(1:nx,1:ny-1,4) + &
au(1:nx,1:ny-1,5) + Fxu(1:nx,2:ny) - Fxu(1:nx,1:ny-1) + Fyu(2:nx+1,1:ny-1) - &
Fyu(1:nx,1:ny-1)

au(1:nx,ny-1,2) = 0.0; au(1:nx,1,3) = 0.0

!au(1,:,1)=au(1,:,1)+tau*ary(1,1)*2/dx; au(nx,:,1)=au(nx,:,1)+tau*ary(1,1)*2/dx

bu = 0.0
bu(1:nx,1:ny-1) = Sx(1:nx,1:ny-1) + (1-uru)*au(1:nx,1:ny-1,1)*ust(1:nx,2:ny)/uru
!bu(1,1:ny-1) = bu(1,1:ny-1) + ua*tau*ary(1,1)*2/dx

uloop: do k = 1,iteru

allocate (diag(nx), updiag(nx-1), lowdiag(nx-1), x(nx))
do i = 1,ny-1
	info = 0; nrhs = 1
	lowdiag(1:nx-1) = -au(2:nx,i,5)
	updiag(1:nx-1) = -au(1:nx-1,i,4)
	diag(1:nx) = au(1:nx,i,1)/uru
	x(1:nx) = bu(1:nx,i) + au(1:nx,i,2)*ust(1:nx,i+2) + au(1:nx,i,3)*ust(1:nx,i)
	!x(1:nx) = Sx(1:nx,i) + (1-uru)*au(1:nx,i,1)*ust(1:nx,i+1)/uru &
	! + au(1:nx,i,2)*ust(1:nx,i+2) + au(1:nx,i,3)*ust(1:nx,i)
	!x(1) = x(1) + ua*tau*ary(1,1)*2/dx

	call dgtsv(nx,nrhs,lowdiag,diag,updiag,x,nx,info)
	ust(1:nx,i+1) = x(1:nx)
enddo
deallocate (diag, updiag, lowdiag, x)

allocate (diag(ny-1), updiag(ny-2), lowdiag(ny-2), x(ny-1))
do i = 1,nx
	info = 0; nrhs = 1
	lowdiag(1:ny-2) = -au(i,2:ny-1,3)
	updiag(1:ny-2) = -au(i,1:ny-2,2)
	diag(1:ny-1) = au(i,1:ny-1,1)/uru
	if(i==1) then
		x(1:ny-1) = bu(i,1:ny-1) + au(i,1:ny-1,4)*ust(i+1,2:ny)
	else if(i==nx) then
		x(1:ny-1) = bu(i,1:ny-1) + au(i,1:ny-1,5)*ust(i-1,2:ny)
	else
		x(1:ny-1) = bu(i,1:ny-1) + au(i,1:ny-1,4)*ust(i+1,2:ny) + au(i,1:ny-1,5)*ust(i-1,2:ny)
	end if

	call dgtsv(ny-1,nrhs,lowdiag,diag,updiag,x,ny-1,info)
	ust(i,2:ny) = x(1:ny-1)
enddo
deallocate (diag, updiag, lowdiag, x)

enddo uloop

!print*, ust

Fxv(1:nx-1,2:ny) = ((rho(2:nx,2:ny)+rho(2:nx,3:ny+1))*ust(1:nx-1,2:ny) + &
	(rho(3:nx+1,2:ny)+rho(3:nx+1,3:ny+1))*ust(2:nx,2:ny))*arx(1,1)/4
Fyv(1:nx,1:ny) = ((rho(1:nx,2:ny+1)+rho(2:nx+1,2:ny+1))*vst(1:nx,1:ny) + &
	(rho(2:nx+1,2:ny+1)+rho(3:nx+2,2:ny+1))*vst(2:nx+1,1:ny))*ary(1,1)/4

Sy(1:nx-1,1:ny) = (pst(1:nx-1,1:ny)-pst(2:nx,1:ny))*ary(1,1) !&
!	- 0.5*(rho(2:nx,2:ny+1)+rho(3:nx+1,2:ny+1))*9.81*dx*dy
!if(nt>100) then
	Sy(1:nx-1,1:ny) = Sy(1:nx-1,1:ny) - (rho(2:nx,2:ny+1)-rho(3:nx+1,2:ny+1))*9.81*dx*dy
!end if
! v from 1 to nx+1; Set a from 2 to nx
av(:,:,:) = 0.0
av(1:nx-1,1:ny-1,2) = Dix + max(-Fxv(1:nx-1,2:ny),0.0)
av(1:nx-1,2:ny,3) = Dix + max(Fxv(1:nx-1,2:ny),0.0)
av(1:nx-1,1:ny,4) = Diy + max(-Fyv(2:nx,1:ny),0.0)
av(1:nx-1,1:ny,5) = Diy + max(Fyv(1:nx-1,1:ny),0.0)

av(1:nx-1,1:ny,1) = av(1:nx-1,1:ny,2) + av(1:nx-1,1:ny,3) + av(1:nx-1,1:ny,4) + &
av(1:nx-1,1:ny,5) + Fxv(1:nx-1,2:ny+1) - Fxv(1:nx-1,1:ny) + Fyv(2:nx,1:ny) - &
Fyv(1:nx-1,1:ny)

av(:,1,1)=av(:,1,1)+tau*arx(1,1)*2/dy; av(:,ny-1,1)=av(:,ny-1,1)+tau*arx(1,1)*2/dy

av(nx-1,1:ny,4) = 0.0; av(1,1:ny,5) = 0.0

bv = 0.0
bv(1:nx-1,1:ny) = Sy(1:nx-1,1:ny) + (1-uru)*av(1:nx-1,1:ny,1)*vst(2:nx,1:ny)/uru

vloop: do k = 1,iterv

allocate (diag(nx-1), updiag(nx-2), lowdiag(nx-2), x(nx-1))
do i = 1,ny
	info = 0; nrhs = 1
	lowdiag(1:nx-2) = -av(2:nx-1,i,5)
	updiag(1:nx-2) = -av(1:nx-2,i,4)
	diag(1:nx-1) = av(1:nx-1,i,1)/uru
	if(i==1) then
		x(1:nx-1) = bv(1:nx-1,i) + av(1:nx-1,i,2)*vst(2:nx,i+1)
	else if(i==ny) then
		x(1:nx-1) = bv(1:nx-1,i) + av(1:nx-1,i,3)*vst(2:nx,i-1)
	else
		x(1:nx-1) = bv(1:nx-1,i) + av(1:nx-1,i,2)*vst(2:nx,i+1) + av(1:nx-1,i,3)*vst(2:nx,i-1)
	end if

	call dgtsv(nx-1,nrhs,lowdiag,diag,updiag,x,nx-1,info)
	vst(2:nx,i) = x(1:nx-1)
enddo
deallocate (diag, updiag, lowdiag, x)

allocate (diag(ny), updiag(ny-1), lowdiag(ny-1), x(ny))
do i = 1,nx-1
	info = 0; nrhs = 1
	lowdiag(1:ny-1) = -av(i,2:ny,3)
	updiag(1:ny-1) = -av(i,1:ny-1,2)
	diag(1:ny) = av(i,1:ny,1)/uru
	x(1:ny) = bv(i,1:ny) + av(i,1:ny,4)*vst(i+2,1:ny) + av(i,1:ny,5)*vst(i,1:ny)
	
	call dgtsv(ny,nrhs,lowdiag,diag,updiag,x,ny,info)
	vst(i+1,1:ny) = x(1:ny)
enddo
deallocate (diag, updiag, lowdiag, x)

enddo vloop

!print*, vst

du=0.0; dv=0.0; a(:,:,:)=0.0; bd=0.0; bp=0.0

du(1:nx,1:ny-1) = arx(1,1)*uru/au(1:nx,1:ny-1,1)
dv(1:nx-1,1:ny) = ary(1,1)*uru/av(1:nx-1,1:ny,1)

a(1:nx,1:ny-1,2) = rho(2:nx+1,2:ny)*arx(1,1)*du(1:nx,1:ny-1)
a(1:nx,2:ny,3) = rho(2:nx+1,3:ny+1)*arx(1,1)*du(1:nx,1:ny-1)
a(1:nx-1,1:ny,4) = rho(2:nx,2:ny+1)*ary(1,1)*dv(1:nx-1,1:ny)
a(2:nx,1:ny,5) = rho(3:nx+1,2:ny+1)*ary(1,1)*dv(1:nx-1,1:ny)
a(1:nx,1:ny,1) = a(1:nx,1:ny,2) + a(1:nx,1:ny,3) + a(1:nx,1:ny,4) + a(1:nx,1:ny,5)

!a(1,1,1) = a(1,1,1) + (10**30)

bp(1:nx,1:ny) = rho(2:nx+1,2:ny+1)*arx(1,1)*(ust(1:nx,1:ny)-ust(1:nx,2:ny+1)) + &
	rho(2:nx+1,2:ny+1)*ary(1,1)*(vst(1:nx,1:ny)-vst(2:nx+1,1:ny))

!bp(1,1) = bp(1,1) + (10**30)*po

ploop: do k = 1,iterp

allocate (diag(nx), updiag(nx-1), lowdiag(nx-1), x(nx))
do i = 1,ny
	info = 0; nrhs = 1
	lowdiag(1:nx-1) = -a(2:nx,i,5)
	updiag(1:nx-1) = -a(1:nx-1,i,4)
	diag(1:nx) = a(1:nx,i,1)
	if(i==1) then
		x(1:nx) = bp(1:nx,i) + a(1:nx,i,2)*pcr(1:nx,i+1)
	else if(i==ny) then
		x(1:nx) = bp(1:nx,i) + a(1:nx,i,3)*pcr(1:nx,i-1)
	else
		x(1:nx) = bp(1:nx,i) + a(1:nx,i,2)*pcr(1:nx,i+1) + a(1:nx,i,3)*pcr(1:nx,i-1)
	end if

	call dgtsv(nx,nrhs,lowdiag,diag,updiag,x,nx,info)
	pcr(1:nx,i) = x(1:nx)
end do
deallocate (diag, updiag, lowdiag, x)

allocate (diag(ny), updiag(ny-1), lowdiag(ny-1), x(ny))
do i = 1,nx
	info = 0; nrhs = 1
	lowdiag(1:ny-1) = -a(i,2:ny,3)
	updiag(1:ny-1) = -a(i,1:ny-1,2)
	diag(1:ny) = a(i,1:ny,1)
	if(i==1) then
		x(1:ny) = bp(i,1:ny) + a(i,1:ny,4)*pcr(i+1,1:ny)
	else if(i==nx) then
		x(1:ny) = bp(i,1:ny) + a(i,1:ny,5)*pcr(i-1,1:ny)
	else
		x(1:ny) = bp(i,1:ny) + a(i,1:ny,4)*pcr(i+1,1:ny) + a(i,1:ny,5)*pcr(i-1,1:ny)
	end if

	call dgtsv(ny,nrhs,lowdiag,diag,updiag,x,ny,info)
	pcr(i,1:ny) = x(1:ny)
enddo
deallocate (diag, updiag, lowdiag, x)

enddo ploop

!print*, ust(:,:)

p = (1-urp)*pst + urp*(pcr+pst)

u(1:nx,2:ny) = (1-uru)*ust(1:nx,2:ny) + uru*(ust(1:nx,2:ny) + &
du(1:nx,1:ny-1)*(pcr(1:nx,1:ny-1)-pcr(1:nx,2:ny)))

v(2:nx,1:ny) = (1-uru)*vst(2:nx,1:ny) + uru*(vst(2:nx,1:ny) + &
dv(1:nx-1,1:ny)*(pcr(1:nx-1,1:ny)-pcr(2:nx,1:ny)))

!print*, u

!print*, p

!print*, ust
!print*, vst(1:nx,:)
ust=u; vst=v; pst=p

!print*, dv
p_res=0.0; m_res=0.0
do i = 1,nx
	do j = 1,ny
		p_res = p_res + Fxu(i,j)
	enddo
enddo
!do i = 1,nx-1
!	do j = 1,ny
!		m_res = m_res + au(i+1,j,2)*ust(i+1,j+2) + au(i+1,j,3)*ust(i+1,j) + &
!		 au(i+1,j,4)*ust(i+2,j+1) + au(i+1,j,5)*ust(i,j+1) + Sx(i+1,)
!	enddo
!enddo
!print*, p_res

if(modulo(nt,10)==0) then

at(:,:,:) = 0.0; bt = 0.0
at(1:nx,1:ny,2) = kt*arx(1,1)*ny/(rho(2:nx+1,2:ny+1)*cp*H) + &
 max(-rho(2:nx+1,2:ny+1)*ust(1:nx,2:ny+1)*arx(1,1),0.0)
at(1:nx,1:ny,3) = kt*arx(1,1)*ny/(rho(2:nx+1,2:ny+1)*cp*H) + &
 max(rho(2:nx+1,2:ny+1)*ust(1:nx,1:ny)*arx(1,1),0.0)
at(1:nx,1:ny,4) = kt*ary(1,1)*nx/(rho(2:nx+1,2:ny+1)*cp*L) + &
 max(-rho(2:nx+1,2:ny+1)*vst(2:nx+1,1:ny)*ary(1,1),0.0)
at(1:nx,1:ny,5) = kt*ary(1,1)*nx/(rho(2:nx+1,2:ny+1)*cp*L) + &
 max(rho(2:nx+1,2:ny+1)*vst(1:nx,1:ny)*ary(1,1),0.0)
at(:,:,1) = at(:,:,2) + at(:,:,3) + at(:,:,4) + at(:,:,5) + &
	rho(2:nx+1,2:ny+1)*arx(1,1)*(ust(1:nx,2:ny+1) - ust(1:nx,1:ny)) + &
	rho(2:nx+1,2:ny+1)*ary(1,1)*(vst(2:nx+1,1:ny) - vst(1:nx,1:ny))

at(:,1,1) = at(:,1,1) + 2*tau*cp*arx(1,1)/(Pr*ny)
at(:,ny,1) = at(:,ny,1) + 2*tau*cp*arx(1,1)/(Pr*ny)

at(1,:,1) = at(1,:,1) + 2*tau*cp*arx(1,1)/(Pr*nx)
at(nx,:,1) = at(nx,:,1) + 2*tau*cp*arx(1,1)/(Pr*nx)

bt(1,:) = 2*tau*cp*arx(1,1)*T0/(Pr*nx)
bt(nx,:) = 2*tau*cp*arx(1,1)*T0/(Pr*nx)

bt(:,1) = bt(:,1) + 2*tau*cp*arx(1,1)*T1/(Pr*ny)
bt(:,ny) = bt(:,ny) + 2*tau*cp*arx(1,1)*T2/(Pr*ny)

tloop: do k = 1,itert

allocate (diag(nx), updiag(nx-1), lowdiag(nx-1), x(nx))
do i = 1,ny
	info = 0; nrhs = 1
	lowdiag(1:nx-1) = -at(2:nx,i,5)
	updiag(1:nx-1) = -at(1:nx-1,i,4)
	diag(1:nx) = at(1:nx,i,1)
	if(i==1) then
		x(1:nx) = bt(1:nx,i) + at(1:nx,i,2)*T(1:nx,i+1)
	else if(i==ny) then
		x(1:nx) = bt(1:nx,i) + at(1:nx,i,3)*T(1:nx,i-1)
	else
		x(1:nx) = bt(1:nx,i) + at(1:nx,i,2)*T(1:nx,i+1) + at(1:nx,i,3)*T(1:nx,i-1)
	end if

	call dgtsv(nx,nrhs,lowdiag,diag,updiag,x,nx,info)
	T(1:nx,i) = x(1:nx)
end do
deallocate (diag, updiag, lowdiag, x)

allocate (diag(ny), updiag(ny-1), lowdiag(ny-1), x(ny))
do i = 1,nx
	info = 0; nrhs = 1
	lowdiag(1:ny-1) = -at(i,2:ny,3)
	updiag(1:ny-1) = -at(i,1:ny-1,2)
	diag(1:ny) = at(i,1:ny,1)
	if(i==1) then
		x(1:ny) = bt(i,1:ny) + at(i,1:ny,4)*T(i+1,1:ny)
	else if(i==nx) then
		x(1:ny) = bt(i,1:ny) + at(i,1:ny,5)*T(i-1,1:ny)
	else
		x(1:ny) = bt(i,1:ny) + at(i,1:ny,4)*T(i+1,1:ny) + at(i,1:ny,5)*T(i-1,1:ny)
	end if

!	call dgtsv(ny,nrhs,lowdiag,diag,updiag,x,ny,info)
!	T(i,1:ny) = x(1:ny)
enddo
deallocate (diag, updiag, lowdiag, x)

enddo tloop
rho(2:nx+1,2:ny+1) = p_abs/(R*(273+T(1:nx,1:ny)))

end if

print*, ust(:,:)

enddo iterloop

!open(1, file = 'data-u.dat', status = 'new')  
!do i=1,nx+1
!	write(1,*) vst(i,1:ny)   
!end do 


end program cavity