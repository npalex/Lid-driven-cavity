! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Roe-solver for the 2D incompressible Navier Stokes equations
! solve Riemann problems along one slice of data.

! waves: 2
! equations: 2

! Conserved quantities:
!       1 x_momentum
!       2 y_momentum

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


! This Riemann solver differs from the original clawpack Riemann solver
! for the interleaved indices

    implicit double precision (a-h,o-z)

    dimension   ql(meqn,           1-mbc:maxm+mbc)
    dimension   qr(meqn,           1-mbc:maxm+mbc)
    dimension    s(mwaves,         1-mbc:maxm+mbc)
    dimension wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    dimension  apdq(meqn,          1-mbc:maxm+mbc)
    dimension  amdq(meqn,          1-mbc:maxm+mbc)

!   # Re set in setprob.f90 file
    common /cparam/ Re

!   # Roe averages quantities of each interface
    parameter (maxm2 = 1800)
    double precision u(-6:maxm2+7),v(-6:maxm2+7)


!   local arrays
!   ------------
    dimension delta(2)
!    logical :: efix

!    data efix /.False./    !# Use entropy fix for transonic rarefactions

!   # Set mu to point to the component of the system that corresponds
!   # to momentum in the direction of this slice, mv to the orthogonal
!   # momentum:
    if (ixy == 1) then
        mu = 1
        mv = 2
    else
        mu = 2
        mv = 1
    endif

!   # Note that notation for u and v reflects assumption that the
!   # Riemann problems are in the x-direction with u in the normal
!   # direciton and v in the orthogonal direcion, but with the above
!   # definitions of mu and mv the routine also works with ixy=2
!   # and returns, for example, f0 as the Godunov flux g0 for the
!   # Riemann problems u_t + g(u)_y = 0 in the y-direction.

	! Compute the Roe-averaged variables needed in the Roe solver.
	! Currently, qr(1,:) = ql(1,:) = u(:) and qr(2,:) = ql(2,:) = v(:)
	! in the rpn2() subroutine call from flux2.f90, g1d are sent for both ql and qr
	! Hence, u(i) and v(i) are edge velocities at the left edge of each cell  
	do 10 i = 2-mbc, mx+mbc
		u(i) = (qr(mu,i-1)+ql(mu,i))*0.50d0
        v(i) = (qr(mv,i-1)+ql(mv,i))*0.50d0
    10 END DO 

	
!   # now split the jump in q at each interface into waves

!   # find a1 and a2, the coefficients of the 2 eigenvectors:
	do 20 i = 2-mbc, mx+mbc
        delta(1) = ql(mu,i) - qr(mu,i-1)
        delta(2) = ql(mv,i) - qr(mv,i-1)
		
		if (u(i) == 0.d0) then
			a1 = delta(2) - 2.0d0*v(i)
			a2 = 2.0d0
		else
			a1 = delta(2) - v(i)/u(i)*delta(1)
			a2 = delta(1)/u(i)
		endif
		
!      # Compute the waves.
        wave(mu,1,i) = 0.0d0
        wave(mv,1,i) = a1
        s(1,i) = u(i)

        wave(mu,2,i) = a2*u(i)
        wave(mv,2,i) = a2*v(i)
        s(2,i) = 2*u(i)

    20 END DO

!    # compute flux differences amdq and apdq.
!    ---------------------------------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do m=1,2
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            end do
        end do
    end do

    return
    end subroutine rpn2


