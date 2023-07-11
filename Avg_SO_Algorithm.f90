program reconst

implicit none

!Reconstruction parameters
integer, parameter :: nyi=163, nzi=48 !Initial condition size
integer, parameter :: nxo=1, nyo=163 !Observation size
integer, parameter :: tobs=2400 !Number of observations available
integer, parameter :: nvel=3 !Number of velocity directions (u,v,z = 3)

!Input variables
real(8), dimension(nyi,nzi,3) :: ui,uout !Initial velocity
real(8), dimension(nxo,nyo,3,tobs) :: uo !Observation velocity at each time
real(8), dimension(nxo,nyi,nzi,3) :: urec !Reconstructed velocity

!Minimization algorithm variables
integer, dimension(nzi) :: tout, tcur, dobj, cmin !Final time, current estimate and derivative of objective function
real(8), dimension(nzi) :: tact, tprev, tstart, chkerr, chkerr2 !Real time for interpolation purposes, minimum t for no convergence
real(8), dimension(nzi) :: tmin, minerr, minerr2, velerr
integer, dimension(nzi,20) :: t_corr
integer, dimension(nzi) :: num,itr
real(8) :: obj,uinter,uintero !Objective Function & interpolated velocity

!Brute Force variables
real(8), dimension(tobs) :: bfcost !Brute force cost
real(8), dimension(nyo,3,tobs) :: uxs,uxso !Observation velocity at X_0 for brute force
real(8), dimension(nyi,3) :: uiav,uiavo !Averaged velocity for BF1 method
real(8), dimension(nyi,3,tobs) :: uls !Least square estimate for BF2

!3D reconstruction variables
integer, parameter :: nx1=241, ny1=241, nz1=48 !Synthetic data size
integer, dimension(nzi) :: div
real(8),dimension(nx1,ny1,nz1,nvel) :: uc
real(8),dimension(nx1,ny1,nz1) :: u,v,w

real :: alp,dt,eps,ngbh,algerr,errmin,start1,fin1
!Weighting for time difference, dt - time step,
!algerr - GS algorithm error, ngbh - time neighbour error, errmin - min error

integer :: i,j,k,ts,dir,cnt,rnd,fnd,obs !Counters ts - time counter;
integer :: bf,tini,mm,iof,err !bf - brute force methodology; tini - initial guess; mm - minimazation algorithm, iof - inlet(1)/outlet+inlet(2)
integer :: count,cmax
character (len=150) :: filename,path

integer*4 timeArray(3)

call cpu_time(start1)
bf = 1; mm = 2; iof = 1

!Initialisation for random variable
call itime(timeArray)
i = rand ( timeArray(1)+timeArray(2)+timeArray(3) )

if (nyi .ne. nyo) then
    print *, "Error! Observation and initial domain mismatch!"
    stop
endif

!Read Observations
uo = 0.
open (15,file='velobs.dat',form='formatted',status='unknown')
do ts=1,tobs
do i=1,1
do j=1,nyo
read(15,*) uo(i,j,1,ts),uo(i,j,2,ts),uo(i,j,3,ts)
enddo
enddo
enddo
close(15)


!Read Initial Condition in a loop
do obs = 30000,30300,10
ui = 0.
write (filename, '( "velini_", I5, ".dat" )' )  obs
open (15,file=trim(filename),form='formatted',status='unknown')
do j=1,nyi
    do k=1,nzi
        read(15,*) ui(j,k,1),ui(j,k,2),ui(j,k,3)
    enddo
enddo
close(15)
print *, filename

if(iof .eq. 2) then
uout = 0.
open (15,file='velout.dat',form='formatted',status='unknown')
do j=1,nyi
do k=1,nzi
read(15,*) uout(j,k,1),uout(j,k,2),uout(j,k,3)
enddo
enddo
close(15)
endif

Print *, 'Data Read and Stored'
!call sleep(2)

tini = 0.
if (mm .eq. 0) then !BRUTE FORCE METHOD
print *,'Brute Force Employed'
do k = 1,nzi
uls = 0.
do ts = 1,tobs
uls(:,:,ts) = (uo(1,:,:,ts) - ui(:,k,:))**2
if (iof .eq. 2) then
uls(:,:,ts) = uls(:,:,ts) + (uo(nxo,:,:,ts) - uout(:,k,:))**2
endif
enddo
bfcost(:) = 0.
do ts = 1,tobs
do j = 1,nyo
do dir = 1,nvel
bfcost(ts) = bfcost(ts) + uls(j,dir,ts)
enddo
enddo
enddo
tout(k) = MINLOC(bfcost, DIM=1, MASK=(bfcost >= 0))
chkerr(k) = bfcost(tout(k))
enddo

do k = 1,nzi
print *, k,tout(k),chkerr(k)
enddo

else
!Observation velocity at X_0
uxs(:,:,:) = uo(1,:,:,:)
if (iof .eq. 2) uxso(:,:,:) = uo(nxo,:,:,:)

!Brute Force for initial guess
!Method 1 - Average along initial condition and identify best time guess
if (bf .eq. 1) then
uiav(:,:) = 0.
do k = 1,nzi
    uiav(:,:) = uiav(:,:) + ui(:,k,:)
enddo
uiav = uiav/nzi

if (iof .eq. 2) then
    uiavo(:,:) = 0.
        do k = 1,nzi
            uiavo(:,:) = uiavo(:,:) + uout(:,k,:)
        enddo
    uiavo = uiavo/nzi
endif

do ts = 1, tobs
    uxs(:,:,ts) = uxs(:,:,ts) - uiav(:,:)
    if (iof .eq. 2) uxs(:,:,ts) = uxs(:,:,ts) + uxso(:,:,ts) - uiavo(:,:)
enddo
bfcost(:) = 0.
do ts = 1, tobs
do j = 1,nyo
    do dir = 1,nvel
        bfcost(ts) = bfcost(ts) + uxs(j,dir,ts)**2
    enddo
enddo
enddo
tini = MINLOC(bfcost, DIM=1, MASK=(bfcost >= 0))

!Method 2 - Calculate minimum error for each time by calculating sum of error between obs and ini @ each z.
elseif (bf .eq. 2) then
uls(:,:,:) = 0.
do ts = 1,tobs
do k = 1,nzi
uls(:,:,ts) = uls(:,:,ts) + (uxs(:,:,ts) - ui(:,k,:))**2
if (iof .eq. 2) uls(:,:,ts) = uls(:,:,ts) + (uxso(:,:,ts) - uout(:,k,:))**2
enddo
enddo
bfcost(:) = 0.
do ts = 1,tobs
do j = 1,nyo
do dir = 1,nvel
bfcost(ts) = bfcost(ts) + uls(j,dir,ts)
enddo
enddo
enddo
tini = MINLOC(bfcost, DIM=1, MASK=(bfcost >= 0))
endif

Print*, 'Brute force Initial Guess:',tini
!call sleep(2)

!Minimization Algorithms
dt = 0.003 !Time step for each observation

tact = tini*dt !Initialising first guess with best guess from brute force
tstart = tact !Saving initial time

tout = INT(tact/dt)
cmax = 50000 !Maximum iterations

minerr = 1000000.0
minerr(24) = 0.

tout(24) = 1
tact(24) = 1.*dt
tmin = tout

!chkerr = 1000.
!fnd = 0
num = 1
t_corr = 0
do cnt = 1, cmax

alp = 0.05 !Coefficient for neighbourhood times/neighbourhood mismatch
!alp = 0.001/((((cmax/1000)-INT(cnt/1000)+1))*1000)

eps = 0.00005 !Fixing value to ensure no drastic reduction
if (mod(cnt,5000) .eq. 0) eps = 0.005 !Initial high value to enable mobility
tprev = tact !Updating current time with updated time
rnd = 0
do k = 1,nzi
    if (k .eq. 24) cycle;
!    if (chkerr(k) .lt. 300.) cycle;
    velerr(k) = 0.
    chkerr(k) = 0.
    chkerr2(k) = 0.
    !Code to caluclate DtU*(Uo-Ui)
    do j = 1,nyo
        do dir = 1,3
            if (mod(tact(k),dt) .gt. dt/50.) then
                uinter = (uo(1,j,dir,tout(k))*(dt-(tact(k)-tout(k)*dt)) + &
                        uo(1,j,dir,tout(k)+1)*(tact(k)-tout(k)*dt))/dt
                if (iof .eq. 2) then
                    uintero = (uo(nxo,j,dir,tout(k))*(dt-(tact(k)-tout(k)*dt)) + &
                        uo(nxo,j,dir,tout(k)+1)*(tact(k)-tout(k)*dt))/dt
                endif
            else
                uinter = uo(1,j,dir,tout(k))
                if (iof .eq. 2) uintero = uo(nxo,j,dir,tout(k))
            endif
            velerr(k) = velerr(k) + 2*(((uo(1,j,dir,tout(k)+1) - uo(1,j,dir,tout(k)))/dt)*&
                (uinter-ui(j,k,dir)))
            chkerr(k) = chkerr(k)+ (uo(1,j,dir,tout(k))-ui(j,k,dir))**2.
            if (iof .eq. 2) then
                velerr(k) = velerr(k) + 2*(((uo(1,j,dir,tout(k)+1) - uo(1,j,dir,tout(k)))/dt)*&
                    (uintero-uout(j,k,dir)))
                chkerr(k) = chkerr(k)+(uintero-uout(j,k,dir))**2.
            endif
        enddo
    enddo

    !Updating tout based on GS method - Not working!
    if (mm .eq. 1) then
        if (k .eq. 1) then
        tact(k) = (1/(2*alp))*(2*alp*(tprev(k+1)) - eps*2*velerr(k))
        elseif (k .gt. 1 .AND. k .lt. nzi) then
        tact(k) = (1/(2*alp))*(2*alp*(tact(k-1)+tprev(k+1)) - eps*2*velerr(k))
        else
        tact(k) = (1/(2*alp))*(2*alp*(tact(k-1)) - eps*2*velerr(k))
        endif
    !Updating based on Gradient Descent method
    elseif (mm .eq. 2) then
!        if (k .eq. 1) then
!            ngbh = 2*alp*((tprev(k)-tprev(k+1)))
!        elseif (k .gt. 1 .AND. k .lt. nzi) then
!            ngbh = 2*alp*((tprev(k)-tprev(k-1))+(tprev(k)-tprev(k+1)))/2.
!        else
!            ngbh = 2*alp*((tprev(k)-tprev(k-1)))
!        endif

!        if (k .eq. 1) then
!            ngbh = 0.
!            do j = 1, 26
!                do dir = 1,3
!                ngbh = ngbh + alp*(uo(1,j,dir,tout(k))-uo(1,j,dir,tout(nzi)))
!                enddo
!            enddo
!            do j = 96, nyo
!                do dir = 1,3
!                ngbh = ngbh + alp*(uo(1,j,dir,tout(k))-uo(1,j,dir,tout(nzi)))
!                enddo
!            enddo
!        elseif (k .gt. 1 .AND. k .le. nzi) then
!            ngbh = 0.
!            do j = 1, 26
!                do dir = 1,3
!                ngbh = ngbh + alp*(uo(1,j,dir,tout(k))-uo(1,j,dir,tout(k-1)))
!                enddo
!            enddo
!            do j = 96, nyo
!                do dir = 1,3
!                ngbh = ngbh + alp*(uo(1,j,dir,tout(k))-uo(1,j,dir,tout(k-1)))
!                enddo
!            enddo
!        endif
!            chkerr2(k) = ngbh
        tact(k) = tact(k) - (eps*(2*velerr(k)))!+ngbh
    endif
    if (tact(k) .gt. tobs*dt) then
        tact(k) = rand(0)*(tobs-1)*dt
        rnd = rnd + 1
    elseif (tact(k) .lt. dt) then
        tact(k) = rand(0)*(tobs-1)*dt
        rnd = rnd + 1
    endif
    if ((chkerr(k)+chkerr2(k)) .lt. minerr(k)) then
        cmin(k) = cnt
        minerr(k) = chkerr(k) + chkerr2(k)
        tmin(k) = tout(k)
        itr(k) = cnt
    endif
!    if (cnt .gt. 20000) then
    if (itr(k) .eq. cnt) then
        num(k) = 1
        t_corr(k,:) = 0
    endif
    if (chkerr(k) .lt. minerr(k)+0.75) then
        err = 0
        do i = 1,num(k)
        if(tout(k) .eq. t_corr(k,i)) err = 1
        enddo
        if (err .ne. 1) then
        t_corr(k,num(k)) = tout(k)
        num(k) = num(k)+1
!        print *, k,num(k), tout(k),cnt,itr(k)
        endif
    endif
!    endif
    tout(k) = INT(tact(k)/dt)
enddo

algerr = 0.
do k = 1,nzi
algerr = algerr + velerr(k) !Algorithm convergence error between previous and current time
enddo

!if (algerr .lt. errmin) then
!    errmin = algerr
!    count = cnt
!    do k = 1,nzi
!        tmin(k) = tact(k)
!        minerr(k) = chkerr(k)
!        minerr2(k) = chkerr2(k)
!    enddo
!endif

print *, cnt, eps, algerr, rnd, tout(1), (eps*(2*velerr(1))), chkerr2(1)

!if (sum(algerr) .lt. 0.000000001) then
if (fnd .eq. nyi) then
    print *, 'Exiting: ',algerr; exit !Stop when convergence is achieved
    call sleep(2)
endif
enddo
!FINISH of Loop
print *, 'eps =',eps

!Post-Proc check for convergence - Choose optimal error in case of no convergence
if (cnt .ne. cmax+1) then
do k = 1,nzi
    print *, k,tstart(k),chkerr(k),chkerr2(k), tout(k)
enddo
endif
!
if (cnt .eq. cmax+1) then
    print*, 'No Convergence - Choosing minimum error',count,sum(minerr)
    do k = 1,nzi
        tout(k) = tmin(k)!INT(tmin(k)/dt)
        print*, k,tstart(k),minerr(k),cmin(k),tout(k),num(k)
    enddo
    call sleep(2)
endif

endif

!do j = 1,nzi
!if (num(j) .eq. 1) then
!    num(j) = 2
!    t_corr(j,1) = tout(j)
!endif
!enddo
num(24) = 2
t_corr(24,1) = 1

!Reconstructing velocity from final times
urec = 0.
do k = 1,nzi
do j = 1,num(k)-1
do dir = 1,3
urec(:,:,k,dir) = urec(:,:,k,dir) + uo(:,:,dir,t_corr(k,j))
enddo
enddo
urec(:,:,k,:) = urec(:,:,k,:)/(num(k)-1)
print *, maxval(urec(:,:,k,:)),num(k)-1
enddo

!Writing reconstructed velocity to file
!open (15,file='time_recon.dat',form='formatted',status='unknown')
!do k=1,nzi
!write(15,*) tout(k)
!enddo
!close(15)
!
!open (15,file='time_avgs.dat',form='formatted',status='unknown')
!do k=1,nzi
!do j = 1,20
!write(15,*) t_corr(k,j)
!enddo
!write(15,*)
!enddo
!close(15)

!open (15,file='velrecon_2d.dat',form='formatted',status='unknown')
!do i=1,nxo
!do j=1,nyi
!do k=1,nzi
!write(15,*) urec(i,j,k,1),urec(i,j,k,2),urec(i,j,k,3)
!enddo
!enddo
!enddo
!close(15)

write (filename, '( "ux_avgrec", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
DO I=1,nxo
WRITE(11,REC=COUNT) urec(I,J,K,1)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)
write (filename, '( "uy_avgrec", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
DO I=1,nxo
WRITE(11,REC=COUNT) urec(I,J,K,2)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)
write (filename, '( "uz_avgrec", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nzi
DO J=1,nyi
DO I=1,nxo
WRITE(11,REC=COUNT) urec(I,J,K,3)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)


t_corr(:,:) = t_corr(:,:) + 599
!Prepping u*_rec for simulation i.e. initial condition *******
div = 0
uc = 0.
!do dir = 1,nvel

do ts = 1,nzi

do cnt = 1,20
if (t_corr(ts,cnt) .eq. 599) exit
print *, ts, cnt

!if (dir == 1) then

if (t_corr(ts,cnt) .lt. 10) write (filename, '( "ux000", I1 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 10 .and. t_corr(ts,cnt) .lt. 100) write (filename, '( "ux00", I2 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 100 .and. t_corr(ts,cnt) .lt. 1000) write (filename, '( "ux0", I3 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 1000 .and. t_corr(ts,cnt) .lt. 10000) write (filename, '( "ux", I4 )' )  t_corr(ts,cnt)

path = '/Volumes/LaCie/pchandra/Incompact3d Simulations/Data Assimilation/Data_Reconstruction/Data-Sets/SS_WkFl/'
print *, trim(path)//filename

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) u(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)
uc(:,:,ts,1) = uc(:,:,ts,1) + u(:,:,24)
div(ts) = div(ts) + 1
!endif

!if (dir == 2) then

if (t_corr(ts,cnt) .lt. 10) write (filename, '( "uy000", I1 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 10 .and. t_corr(ts,cnt) .lt. 100) write (filename, '( "uy00", I2 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 100 .and. t_corr(ts,cnt) .lt. 1000) write (filename, '( "uy0", I3 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 1000 .and. t_corr(ts,cnt) .lt. 10000) write (filename, '( "uy", I4 )' )  t_corr(ts,cnt)

print *, trim(path)//filename

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) v(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)
uc(:,:,ts,2) = uc(:,:,ts,2) + v(:,:,24)
!div(ts,dir) = div(ts,dir) + 1
!endif

!if (dir == 3) then

if (t_corr(ts,cnt) .lt. 10) write (filename, '( "uz000", I1 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 10 .and. t_corr(ts,cnt) .lt. 100) write (filename, '( "uz00", I2 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 100 .and. t_corr(ts,cnt) .lt. 1000) write (filename, '( "uz0", I3 )' )  t_corr(ts,cnt)
if (t_corr(ts,cnt) .ge. 1000 .and. t_corr(ts,cnt) .lt. 10000) write (filename, '( "uz", I4 )' )  t_corr(ts,cnt)

print *, trim(path)//filename

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) w(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

uc(:,:,ts,3) = uc(:,:,ts,3) + w(:,:,24)
!div(ts,dir) = div(ts,dir) + 1
!endif

enddo
enddo

do ts = 1,nzi
uc(:,:,ts,:) = uc(:,:,ts,:)/div(ts)
enddo

uc(108,40:202,:,:) = ui(:,:,:)

!Reading and adding center plane
path = '/Users/pchandra/Documents/Incompact3d_Simulations/Data_Assimilation&
/DA_LatestCode/Preproc_Dataset/6D6D_Rec_SS/4Dvar_Pre_Rec/'
write (filename, '( "ux", I5 )' )  obs

print *, trim(path)//filename

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) u(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uy", I5 )' )  obs

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) v(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uz", I5 )' )  obs

OPEN(11,FILE=trim(path)//filename,FORM='UNFORMATTED', ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
READ(11,REC=COUNT) w(I,J,K)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

uc(:,:,24,1) = u(:,:,24)
uc(:,:,24,2) = v(:,:,24)
uc(:,:,24,3) = w(:,:,24)


!ui = uc
!do j2 = 1,nyi
!do j = 1,nyi
!if (ysim(j2) .eq. yp(j)) then
!    ui(:,j2,:,:) = uc(:,j,:,:)
!elseif (ysim(j2) .gt. yp(j) .AND. ysim(j2) .lt. yp(j+1)) then
!    ui(:,j2,:,:) = ((ysim(j2)-yp(j))*(uc(:,j+1,:,:)-uc(:,j,:,:))/&
!                        (yp(j+1)-yp(j)))+uc(:,j,:,:)
!endif
!enddo
!enddo

write (filename, '( "ux_recavg_fd", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
WRITE(11,REC=COUNT) uc(I,J,K,1)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uy_recavg_fd", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
WRITE(11,REC=COUNT) uc(I,J,K,2)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uz_recavg_fd", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
DO J=1,ny1
DO I=1,nx1
WRITE(11,REC=COUNT) uc(I,J,K,3)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "ux_recavg_assim_d", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
do j=40,202
DO I=108,180
WRITE(11,REC=COUNT) uc(I,J,K,1)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uy_recavg_assim_d", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
do j=40,202
DO I=108,180
WRITE(11,REC=COUNT) uc(I,J,K,2)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

write (filename, '( "uz_recavg_assim_d", I5 )' )  obs
OPEN(11,FILE=trim(filename),FORM='UNFORMATTED',&
ACCESS='DIRECT', RECL=8)
COUNT = 1
DO K=1,nz1
do j=40,202
DO I=108,180
WRITE(11,REC=COUNT) uc(I,J,K,3)
COUNT = COUNT + 1
ENDDO
ENDDO
ENDDO
CLOSE(11)

uc = 0.
enddo

call cpu_time(fin1)
print *, 'Total time taken - ',fin1-start1



end program reconst
