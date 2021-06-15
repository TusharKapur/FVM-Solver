psi = 0.0
psiloop: do k = 1,iterpsi

allocate (diag(nx-1), updiag(nx-2), lowdiag(nx-2), x(nx-1))
do i = 1,ny-1
	info = 0; nrhs = 1
	lowdiag(1:nx-2) = -ary(1,1)/dx
	updiag(1:nx-2) = -ary(1,1)/dx
	diag(2:nx-2) = 2*ary(1,1)/dx
	diag(1) = 2*ary(1,1)/dx; diag(nx-1) = 2*ary(1,1)/dx
	
	if(i==1) then
		x(1:nx-1) = dx*dy*((ust(2:nx,i+1)-ust(1:nx-1,i+1))/dx - (vst(2:nx,i+1)-vst(2:nx,i))/dy) + &
			arx(1,1)/dy*psi(1:nx-1,i+1)
		diag(:) = diag(:) + 1*arx(1,1)/dy
	else if(i==ny-1) then
		x(1:nx-1) = dx*dy*((ust(2:nx,i+1)-ust(1:nx-1,i+1))/dx - (vst(2:nx,i+1)-vst(2:nx,i))/dy) + &
			arx(1,1)/dy*psi(1:nx-1,i-1)
		diag(:) = diag(:) + 1*arx(1,1)/dy
	else
		x(1:nx-1) = dx*dy*((ust(2:nx,i+1)-ust(1:nx-1,i+1))/dx - (vst(2:nx,i+1)-vst(2:nx,i))/dy) + &
			arx(1,1)/dy*psi(1:nx-1,i+1) + arx(1,1)/dy*psi(1:nx-1,i-1)
	end if

	call dgtsv(nx-1,nrhs,lowdiag,diag,updiag,x,nx-1,info)
	print*, info
	psi(1:nx-1,i) = x(1:nx-1)
end do
deallocate (diag, updiag, lowdiag, x)

enddo psiloop

print*, psi