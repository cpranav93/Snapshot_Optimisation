program reconst

implicit none

!Reconstruction parameters
integer, parameter :: nyi=130, nzi=48 !Initial condition size
integer, parameter :: nxo=60, nyo=130 !Observation size
integer, parameter :: tobs=2500 !Number of observations available
integer, parameter :: nvel=3 !Number of velocity directions (u,v,z = 3)

!Input variables
real(8), allocatable :: u_or(:,:),phi(:,:),u_rec(:,:),umean(:,:), u_ini(:,:)
real ( kind = 8 ), allocatable :: a_save(:,:), a(:,:), tco(:,:), id(:,:)
real(8) :: err

!Counters
integer, parameter :: m =nyi*3, n = tobs
integer :: tini, pos, ts, count
integer :: i,j,k

character (len=90) :: filename,path

!SVD Parameters
integer ( kind = 4 ) info
character jobu
character jobv
integer ( kind = 4 ) lda
integer ( kind = 4 ) ldu
integer ( kind = 4 ) ldv
integer ( kind = 4 ) lwork
real ( kind = 8 ) sn(m)
real ( kind = 8 ) un(m,m)
real ( kind = 8 ) v(m,n)
real ( kind = 8 ) work(5*n+m)

!open (15,file='velini.dat',form='formatted',status='unknown')
!do j=1,nyi
!do k=1,nzi
!write(15,*) ui(j,k,1),ui(j,k,2),ui(j,k,3)
!enddo
!enddo
!close(15)
allocate (u_or(nyo*3,tobs))
allocate (u_ini(nyo*3,nzi))

Print *, 'Reading Observation File'

open (15,file='velobs.dat',form='formatted',status='unknown')
do ts=1,tobs
do j=1,nyo
read(15,*) u_or(j,ts),u_or(nyo+j,ts),u_or(2*nyo+j,ts)
enddo
enddo

open (15,file='velini.dat',form='formatted',status='unknown')
do j=1,nyi
do k=1,nzi
read(15,*) u_ini(j,k),u_ini(nyi+j,k),u_ini(2*nyi+j,k)
enddo
enddo
close(15)


print*, u_ini(1:3,2)
print*, u_or(1:3,1)
!Subtracting mean from the SS
print *, 'Subtracting the mean field', maxval(u_or)
allocate (umean(nyi*3,1))
umean(:,1) = (sum(u_or,dim=2)+sum(u_ini,dim=2))/(tobs+nzi)

open (15,file='umean.dat',form='formatted',status='unknown')
do ts=1,1
do j=1,nyo
write(15,*) umean(j,ts),umean(nyo+j,ts),umean(2*nyo+j,ts)
enddo
enddo

do i = 1,tobs
u_or(:,i) = u_or(:,i) - umean(:,1)
enddo

do i = 1,nzi
u_ini(:,i) = u_ini(:,i) - umean(:,1)
enddo

allocate ( a(m,n) )
allocate ( a_save(m,n) )

Print *, 'Performing SVD', maxval(u_or)
a = u_or
!a_save = a

jobu = 'a'
jobv = 's'
lda = m
ldu = m
ldv = m
lwork = 5 * n + m

call dgesvd ( jobu, jobv, m, n, a, lda, sn, un, ldu, v, ldv, work, lwork, &
info )


if ( info == 0 ) then
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'SVD_TRUNCATED_U:'
write ( *, '(a)' ) '  DGESVD computation was successful.'
else
write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'SVD_TRUNCATED_U - Warning!'
write ( *, '(a,i8)' ) '  DGESVD returned INFO = ', info
end if

allocate (u_rec(nyi*3,tobs))

do i = 1, nyi*3
do j = 1, tobs
u_rec(i,j) = 0.0D+00
do k = 1, nyi*3
u_rec(i,j) = u_rec(i,j) + un(i,k) * sn(k) * v(k,j)
end do
end do
end do

err = maxval ( abs ( u_rec(1:m,1:n) - u_or(1:m,1:n) ) )

write ( *, '(a)' ) ' '

write ( *, '(a,g14.6)' ) '  Maximum error |A - U*S*V''| = ', err

deallocate (a)
deallocate (a_save)

path = ''

open (15,file=trim(path)//'U_USV.dat',form='formatted',status='unknown')
do i=1,m
do j=1,m
write(15,*) un(i,j)
enddo
enddo
close(15)

open (15,file=trim(path)//'S_USV.dat',form='formatted',status='unknown')
do j=1,m
write(15,*) sn(j)
enddo
close(15)

open (15,file=trim(path)//'VT_USV.dat',form='formatted',status='unknown')
do i=1,m
do j=1,n
write(15,*) v(i,j)
enddo
enddo
close(15)

allocate ( tco(m,n) )

do i = 1,m
do j = 1,n
tco(i,j) = (sn(i)) * v(i,j)
enddo
enddo

print *, 'Time Coefficient - Done'

allocate (phi(m,m))

phi = un
!do i = 1,m
!do j = 1,m
!phi(i,j) = 0.
!do k = 1,n
!phi(i,j) = phi(i,j) + u_or(i,k)*un(k,j)
!enddo
!enddo
!enddo

!do i = 1,nyi*3
!do j = 1,n
!phi(i,j) = phi(i,j)/((sn(j)**(1./2.)))
!enddo
!enddo

!allocate (u_rec(nyi*3,tobs))

do i = 1, nyi*3
do j = 1, tobs
u_rec(i,j) = 0.0D+00
do k = 1, nyi*3
u_rec(i,j) = u_rec(i,j) + phi(i,k) *tco(k,j)
end do
end do
end do

err = maxval ( abs ( u_rec(:,:) - u_or(:,:) ) )

write ( *, '(a)' ) ' '

write ( *, '(a,g14.6)' ) '  Maximum error SVD = ', err

write (filename, '( "spat_modes" )' )

print *,trim(path)//filename


OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
do i=1,nyi*3
do j=1,nyi*3
WRITE(11,REC=COUNT) phi(i,j)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)

print *, 'Writing Spatial Modes - Done'

write (filename, '( "time_coeff_obs" )' )

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
do i=1,nyi*3
do j=1,tobs
WRITE(11,REC=COUNT) tco(i,j)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)

deallocate(u_rec)
deallocate (tco)
deallocate(u_or)

allocate (tco(m,nzi))
!allocate(u_or(nyi*3,nzi))

!open (15,file='velini.dat',form='formatted',status='unknown')
!do j=1,nyi
!do k=1,nzi
!read(15,*) u_or(j,k),u_or(nyi+j,k),u_or(2*nyi+j,k)
!enddo
!enddo
!close(15)

tco = 0.

do k = 1,nzi
do i = 1,m
tco(i,k) = 0.0
do j = 1,m
tco(i,k) = tco(i,k) + phi(j,i)*u_ini(j,k)
enddo
enddo
enddo

allocate (u_rec(nyi*3,nzi))

do i = 1, nyi*3
do j = 1, nzi
u_rec(i,j) = 0.0D+00
do k = 1, nyi*3
u_rec(i,j) = u_rec(i,j) + phi(i,k)*tco(k,j)
end do
end do
end do

err = maxval ( abs ( u_rec(:,:) - u_ini(:,:) ) )

write ( *, '(a)' ) ' '

write ( *, '(a,g14.6)' ) '  Maximum error SVD = ', err

write (filename, '( "time_coeff_inlet" )' )

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
do i=1,nyi*3
do j=1,nzi
WRITE(11,REC=COUNT) tco(i,j)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)

k = 24
print *, u_ini(1:3,k)
print *, u_rec(1:3,k)

print *, u_ini(nyi+1:nyi+3,k)
print *, u_rec(nyi+1:nyi+3,k)

print *, u_ini(2*nyi+1:2*nyi+3,k)
print *, u_rec(2*nyi+1:2*nyi+3,k)

OPEN(11,FILE='ux_pod',FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
WRITE(11,REC=COUNT) u_rec(J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)
OPEN(11,FILE='uy_pod',FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
WRITE(11,REC=COUNT) u_rec(nyi+J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)
OPEN(11,FILE='uz_pod',FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
WRITE(11,REC=COUNT) u_rec(2*nyi+J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
CLOSE(11)

deallocate(u_rec)
deallocate(umean)
deallocate(u_ini)
deallocate(phi)
deallocate (tco)

print *, 'Writing Time Coefficients - Done'


end program reconst
