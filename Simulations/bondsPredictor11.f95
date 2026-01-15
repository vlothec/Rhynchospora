subroutine extruder(N,passElements,steps,nLEFs,ncentromeres,lifetime,birth,centromeresHalt,centromeres,loopBondDictFinal)
    implicit none

    integer, intent(in) :: N,nLEFs,steps,ncentromeres,passElements
    integer, intent(in), dimension(ncentromeres) :: centromeres
    integer, intent(in), dimension(nLEFs) :: lifetime,birth,centromeresHalt
    integer, dimension(ncentromeres) :: centromeresplus,centromeresminus
    integer, intent(out), dimension(passElements) :: loopBondDictFinal
    integer, dimension(2,nLEFs,2) :: loopBondDict
    logical :: born=.false.,convergence=.false.,halt
    integer :: passSteps,b,i,j,left,right,attempt,pos,s
    integer :: death,halted,longer
    real :: uniqueSeed,randNorm,random_normal,randExp,loopSize

    print*, N,steps,passElements,nLEFs,ncentromeres
    print*, centromeres
    passSteps = (passElements/nLEFs)/2

    do i=1, ncentromeres
        centromeresplus(i) = centromeres(i) + 1
        centromeresminus(i) = centromeres(i) - 1
    end do

    ! assigning a random position to the birth of condensins
    !do i=1,nLEFs
    !    do j=1,2
    loopBondDict(:,:,:)=0
    !    end do
    !end do

    !open(10,file="loopBondDict.bin",form="unformatted",status="new")
!   do loop extrusion of condensins I for #steps
    print*, "Now performing loop extrusion of condensins.", steps
    do b=1, steps
        !print*, b,loopBondDict(b,:,:)
        do i=1, nLEFs
            if (birth(i)==b) then
                born = .false.
                do while (.not.born)
                    call init_random_seed()
                    call random_number(uniqueSeed)
                    loopBondDict(2,i,1) = int(uniqueSeed*(N-1))
                    loopBondDict(2,i,2) = loopBondDict(2,i,1) + 1
                    if (any(loopBondDict(2,1:i-1,:)==loopBondDict(2,i,1)) .or. &
                        any(loopBondDict(2,1:i-1,:)==loopBondDict(2,i,2)) .or. &
                        any(centromeres==loopBondDict(2,i,1)) .or. &
                        any(centromeres==loopBondDict(2,i,2)) .or. &
                        loopBondDict(2,i,1) <= 0 .or. &
                        loopBondDict(2,i,2) > N) then
                        born = .false.
                    else
                        born = .true.
                    end if
                end do
            end if
        end do
        do i=1, nLEFs
            if (birth(i)<b) then
                left = loopBondDict(1,i,1) - 1
                if (any(loopBondDict(1,1:nLEFs,1)==left) .or. &
                    any(loopBondDict(1,1:nLEFs,2)==left) .or. &
                    left==0 .or. &
                    any(centromeres==left)) then
                    left = loopBondDict(1,i,1)
                end if
                loopBondDict(2,i,1) = left
            end if
        end do
        do i=1, nLEFs
            if (birth(i)<b) then
                right = loopBondDict(1,i,2) + 1
                if (any(loopBondDict(1,1:nLEFs,2)==right) .or. &
                    any(loopBondDict(2,1:nLEFs,1)==right) .or. &
                    right==N+1 .or. &
                    any(centromeres==right)) then
                    right = loopBondDict(1,i,2)
                 end if
                loopBondDict(2,i,2) = right
            end if
        end do
! cheking death probability and rebirth
        death=0
        halted=0
        do i=1, nLEFs
! checking if condensin is not stalled by a centromere
            if (birth(i)<b) then
                halt=.false.
                longer=1
                !if (meta==1) then
                do j=1, ncentromeres-1
                    if (loopBondDict(2,i,2)==centromeresplus(j) .and. centromeresHalt(i)==1) halt = .true.
                    if (loopBondDict(2,i,2)==centromeresminus(j) .and. centromeresHalt(i)==1) halt = .true.
                    if (loopBondDict(2,i,1)==centromeresplus(j) .and. centromeresHalt(i)==1) halt = .true.
                    if (loopBondDict(2,i,1)==centromeresminus(j) .and. centromeresHalt(i)==1) halt = .true.
                end do
                if ((loopBondDict(2,i,2)-loopBondDict(2,i,1))>10*lifetime(i)) halt = .false.
                !if (b>30000 .and. lifetime(i)==1000) halt=.true.
                if (halt) halted = halted + 1
                !end if
                if(.not.halt) then
                    call init_random_seed()
                    call random_number(uniqueSeed)
                    if (uniqueSeed < 1./(lifetime(i)*longer)) then
                        !print*, i, 'died'
                        born = .false.
                        death = death + 1
                        do while (.not.born)
                            call init_random_seed()
                            call random_number(uniqueSeed)
                            left = int(uniqueSeed*(N-1))
                            right = left + 1
                            if (any(loopBondDict(2,1:nLEFs,1:2)==left) .or. &
                                any(loopBondDict(2,1:nLEFs,1:2)==right) .or. &
                                any(centromeres==left) .or. &
                                any(centromeres==right) .or. &
                                left <= 0 .or. &
                                right > N) then
                                born = .false.
                            else
                                born = .true.
                            end if
                        end do
                        loopBondDict(2,i,1) = left
                        loopBondDict(2,i,2) = right
                    end if
                end if
            end if
        end do
        !print*, death, "condensins I died in step.",b+1
        loopSize = 0.
        do i=1, nLEFs
            loopSize = loopSize + loopBondDict(2,i,2) - loopBondDict(2,i,1)
        end do
        loopSize = loopSize/nLEFs
        !print*, "new step",b
        do i=1, nLEFS
            !print*, i,loopBondDict(2,i,1),loopBondDict(2,i,2)
            do j=1,2
                loopBondDict(1,i,j) = loopBondDict(2,i,j)
            end do
        end do
        if (b>steps-passSteps) then
            s = b-steps+passSteps-1
            !print*, b,s
            do i=1, nLEFs
                !print*,loopBondDict(2,i,1),loopBondDict(2,i,2)
                loopBondDictFinal(s*2*nLEFs+2*i-1) = loopBondDict(2,i,1)
                loopBondDictFinal(s*2*nLEFS+2*i) = loopBondDict(2,i,2)
            end do
        end if
        !write(10) loopBondDict(1,:,:)
        !print*, death, "condensins died in step",b, " , ",halted,"were halted and mean loop size is", loopSize
    end do
    do i=1, nLEFS
        print*,loopBondDict(2,i,:)
    end do
    print*, "finished predicting bonds."
    !close(10)
end subroutine extruder

subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
        read(un) seed
        close(un)
    else
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
        call system_clock(count)
        if (count /= 0) then
            t = transfer(count, t)
        else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
                seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
        else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
    end if
    call random_seed(put=seed)
end subroutine init_random_seed

REAL FUNCTION random_normal()

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

    IMPLICIT NONE

!     Local variables
    REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
    r1 = 0.27597, r2 = 0.27846, half = 0.5, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

    DO
        call init_random_seed()
        CALL RANDOM_NUMBER(u)
        call init_random_seed()
        CALL RANDOM_NUMBER(v)
        v = 1.7156 * (v - half)

!     Evaluate the quadratic form
        x = u - s
        y = ABS(v) - t
        q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
        IF (q < r1) EXIT
!     Reject P if outside outer ellipse
        IF (q > r2) CYCLE
!     Reject P if outside acceptance region
        IF (v**2 < -4.0*LOG(u)*u**2) EXIT
    END DO

!     Return ratio of P's coordinates as the normal deviate
    random_normal = v/u
    RETURN

END FUNCTION random_normal

! Random smaple from an exponential distribution
!
FUNCTION rand_exponential(mean) RESULT(c)
real :: mean,c,temp
IF (mean <= 0.0d0) THEN
    WRITE(*,*) "mean must be positive"
ELSE
    call init_random_seed()
    CALL RANDOM_NUMBER(temp)
    c=-mean*log(temp)
END IF
END FUNCTION

