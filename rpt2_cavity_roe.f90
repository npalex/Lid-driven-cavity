! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
    implicit double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the incompressible 
!	  # Navier Stokes equations.
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    dimension     ql(meqn, 1-mbc:maxm+mbc)
    dimension     qr(meqn, 1-mbc:maxm+mbc)
    dimension   asdq(meqn, 1-mbc:maxm+mbc)
    dimension bmasdq(meqn, 1-mbc:maxm+mbc)
    dimension bpasdq(meqn, 1-mbc:maxm+mbc)

	dimension waveb(2,2),sb(2)

!   # Re set in setrun.py
    common /cparam/ Re

!   # Roe averages quantities of each interface
    parameter (maxm2 = 1800)
    double precision u(-6:maxm2+7),v(-6:maxm2+7)
	
		! Compute the Roe-averaged variables needed in the Roe solver.
	! Currently, qr(1,:) = ql(1,:) = u(:) and qr(2,:) = ql(2,:) = v(:)
	! in the rpn2() subroutine call from flux2.f90, g1d are sent for both ql and qr
	! Hence, u(i) and v(i) are edge velocities at the left edge of each cell  
	do 10 i = 2-mbc, mx+mbc
		u(i) = (qr(mu,i-1)+ql(mu,i))*0.50d0
        v(i) = (qr(mv,i-1)+ql(mv,i))*0.50d0
    10 END DO 
	
	if (ixy == 1) then
        mu = 1
        mv = 2
    else
        mu = 2
        mv = 1
    endif

	do 20 i = 2-mbc, mx+mbc
		if (v(i) == 0.d0) then
			a1 = asdq(mu,i) - 2.0d0*u(i)
			a2 = 2.0d0
		else
			a1 = asdq(mu,i) - u(i)/v(i)*asdq(mv,i)
			a2 = asdq(mv,i)/v(i)
		endif
		
!      # Compute the waves.
        waveb(mu,1) = a1
        waveb(mv,1) = 0.0d0
        sb(1) = v(i)

        waveb(mu,2) = a2*v(i)
        waveb(mv,2) = a2*u(i)
        sb(2) = 2*v(i)
    
!       # compute the flux differences bmasdq and bpasdq
        do m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do mw=1,mwaves
                bmasdq(m,i) = bmasdq(m,i) &
                + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(m,i) = bpasdq(m,i) &
                + dmax1(sb(mw), 0.d0) * waveb(m,mw)
            end do
        end do
    
    20 END DO
	
    return
    end subroutine rpt2
