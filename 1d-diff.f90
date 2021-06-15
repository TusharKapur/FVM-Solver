  !1-D temp adv-diff
program advdiff
implicit none

       double precision, parameter :: m_dot=0.0
       double precision, parameter :: rho=1.0
       integer, parameter :: nx=6, n = nx-2
     double precision, parameter :: L=1.00
     double precision, parameter :: dx=L/(nx-1)
     double precision at(nx-1,3), u(nx-1), T(nx-1), DiT(nx), F(nx), SuT(nx-1), SpT(nx-1), ar(2*nx), &
     diagt(n+1), lowdiagt(n), updiagt(n), xt(n+1), bt(n+1)

     double precision :: ur=1.0, ua=0.0, k = 235.0
     integer :: i, j, nt, time, nrhst, infot, ldat, ldbt

     double precision, dimension (:,:), allocatable :: Tmat
     allocate ( Tmat(n+1,n+1) )


    arallocationloop: do i = 2, 2*nx
                    ar(i) = 0.5 - i*(0.5-0.1)*dx/(2*L)
    end do arallocationloop
    ar(1) = 0.5
    ar(2*nx-1) = 0.1
    ar(:) = 1.0

    time = 1     

     ! U Initial conditions
     usetloop: do i = 1,nx-1

            u(i) = 0
            T(i) = 0
            
     end do usetloop

     ! BCS
     T(1) = 1.0
     T(nx-1) = 0

     F(:) = 0.0

timeloop: do nt = 1,1


     ditloop: do i = 1,nx

            !F(i) = rho*(u(i)+u(i-1))*ar(2*i-1)/2
            DiT(i) = k*nx*ar(2*i-1)/(2*L)  ! kA/(x1-x2)

     end do ditloop
     DiT(:) = 0.5

     SpT(:) = 0
     SuT(:) = 0
     SpT(1) = -(2*DiT(1) + F(1))
     SuT(1) = (2*DiT(1) + F(1))*T(1)
     SpT(nx-1) = -2*DiT(nx-1)
     SuT(nx-1) = 2*DiT(nx-1)*T(nx-1)

     atloop: do i = 2,nx-2

            at(i,1) = DiT(i) + max(F(i),0.0)
            at(i,3) = DiT(i+1) + max(-F(i+1),0.0)
            at(i,2) = at(i,1) + at(i,3) + F(i+1) - F(i) - SpT(i)

     end do atloop

     at(1,1) = 0.0
     at(1,3) = DiT(2)
     at(1,2) = at(1,1) + at(1,3) + F(2) - F(1) - SpT(1)

     at(nx-1,1) = DiT(nx-1) + F(nx-1)
     at(nx-1,3) = 0.0
     at(nx-1,2) = at(nx-1,1) + at(nx-1,3) + F(nx) - F(nx-1) - SpT(nx-1)

     ! Assembling the Tmat matrix
     Tmat(:,:) = 0.0
     Tmat(1,1) = at(1,2)
     Tmat(1,2) = -at(1,3)
     Tmat(nx-1,nx-2) = -at(nx-1,1)
     Tmat(nx-1,nx-1) = at(nx-1,2)
     
     tloop1: do i = 2,nx-2
            
            Tmat(i,i) = at(i,2)
            Tmat(i,i+1) = -at(i,3)
            Tmat(i,i-1) = -at(i,1)

     end do tloop1


     diagt(nx-1) = Tmat(nx-1,nx-1)
     do i = 1,nx-2

            diagt(i) = Tmat(i,i)
            updiagt(i) = Tmat(i,i+1)
            lowdiagt(i) = Tmat(i+1,i)
     end do

     bt(:) = 0
     bt(1) = SuT(1)
     bt(nx-1) = SuT(nx-1)

     nrhst = 1
    ldat = nx-1
    ldbt = nx-1
     xt = bt
     call dgtsv(nx-1,nrhst,lowdiagt,diagt,updiagt,xt,ldbt,infot)
     T = xt

  time = time + 1

end do timeloop

print*, xt

write (*,202)
do i = 1,n+1
write (*,201)  (Tmat(i,j),j=1,n+1)
end do
200 format (' Computing Inverse mpmattrix ',/,/, &
' Mpmattrix pmat')
201 format (6f12.6)
202 format (/,' Inverse mpmattrix pmat^{-1}')

! output data into a file 
  open(1, file = 'data1.dat', status = 'old')  
  do i=1,nx-1
  write(1,*) T(i)
  end do

deallocate (Tmat)

end program advdiff