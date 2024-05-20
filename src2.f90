subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

	!------------------------------------------------------------------
	!
	!	Written by: Nathan Alexander (May 2024)
	!
	!   The purpose of this subroutine is to peform steps 2-6 of the 
	!		numerical algorithm described in README.md
	!	
	!	Step 2: Update the cell-centered velocity distribution q 
	!			for diffusion according to 
	!
	!				u_t = 1/Re*(u_xx + u_yy)
	!
	!			and
	!
	!				v_t = 1/Re*(v_xx + v_yy),
	!
	!			using an alternating dimension implicit (ADI) method.
	!
	!	Step 3: Update the edge velocities for diffusion via 
	!			linear interpolation
	!
	!	Step 4: Update the pressure distribution by solving the Poisson
	!			equation with a fictitious source:
	!
	!				p_xx + p_yy + c_0*p =  1/dt*(u_x + v_y)
	!
	!	Step 5: Update the edge velocities for pressure-driven flow
	!
	!	Step 6: Update the cell-centered velocities for pressure-driven flow
	!	
	!	Inputs:
	!		meqn   - number of equations
	!		mbc    - number of ghost cells per side of spatial domain
	!		mx     - number of grid cells in x direction
	!		my     - number of grid cells in y direction
	!		xlower - lower bound of spatia domain in x-direction
	!		ylower - lower bound of spatia domain in y-direction
	!		dx     - grid spacing in x-direction
	!		dy     - grid spacing in y-direction
	!		q      - solution vector q = (u, v, p, du/dx + dv/dy)
	!		maux   - not used
	!		aux    - not used
	!		t      - current time
	!		dt     - time interval
	!
	!	Output:
	!			q - updated solution vector q = (u, v, p, du/dx + dv/dy)
	!
	!------------------------------------------------------------------
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
	real(kind=8) :: u_edge(1-mbc:mx,1-mbc:my+mbc), v_edge(1-mbc:my+mbc,1-mbc:my)

	real(kind=8) :: Re	
	common /cparam/ Re 
	
	integer :: i, j, k, info
	
	!-- ADI method memory allocation
	real(kind=8) :: a1(mx*my,mx*my), a2(mx*my,mx*my), c1(mx*my,mx*my), c2(mx*my,mx*my)
	real(kind=8) :: a1_main(mx*my), a1_ld(mx*my-1), a1_ud(mx*my-1), a1_main_v(mx*my), a1_ld_v(mx*my-1), a1_ud_v(mx*my-1)
	real(kind=8) :: error, b(2, mx*my), b2(2, mx*my)
	real(kind=8) :: ab(mx+1,mx*my), q_col(2, mx*my), a2_bs(mx+1,mx*my), a2_bs_v(mx+1,mx*my)
	
	!-- Poisson solver memory allocation
	real(kind=8) :: a_bs(mx+1,mx*my), bu(mx*my,mx*my), bv(mx*my,mx*my)
	real(kind=8) :: bp(mx*mx), bp2(mx*mx-1)
	real(kind=8) :: beta
	integer :: ipiv(mx*my-1)

	!-- initialize arrays
	a1 = 0.d0
	a1_ld = 0.d0
	a1_main = 0.d0
	a1_ud = 0.d0
	a1_ld_v = 0.d0
	a1_main_v = 0.d0
	a1_ud_v = 0.d0	

	a2 = 0.d0
	a2_bs = 0.d0
	a2_bs_v = 0.d0

	c1 = 0.d0
	c2 = 0.d0
	
	b = 0.d0
	b2 = 0.d0

	ab = 0.d0
	q_col = 0.d0
	
	bu = 0.d0
	bv = 0.d0
	a_bs = 0.d0
	
	beta = dx/dt/2.d0
	
	!----------------------------------
	! calculate ADI method matrices
	!----------------------------------
	call adi_matrices(a1, c1, a2, a2_bs, c2, b, b2 ,a1_ld ,a1_main ,a1_ud ,mx, my ,dx ,dy ,dt)
	!-- Note, a second set of a1 diagonals and a second a2 matrix in banded storage is required
	!	because the originals are overwritten by each call to DGTSV() and DPBSV(), respectively
	a1_ld_v = a1_ld
	a1_main_v = a1_main
	a1_ud_v = a1_ud
	a2_bs_v = a2_bs
	
	!--------------------------------------
	! calculate Poisson solver matrices
	!--------------------------------------
	call poisson_matrices(bu, bv, a_bs, mx, my, dx, dy, dt)
	
	!---------------------------------------------------------
    !	Step 2: Update the cell-centered velocity distribution 
	!			for diffusion (here, q = (u,v) pressure not 
	!			updated during this step)
	!		a) solve a1*q_intm = c1*q_n + b
	!		b) solve a2*q_star = c2*q_intm + b
    !---------------------------------------------------------
	
	!-- recast q in form required by dgemv()
	do k = 1, 2
		do j = 1, my
			q_col(k,((j-1)*mx+1):(j*mx)) = q(k,1:mx,j)
		enddo
	enddo	

	!-- calcualte b = c1*q_n + b
	do k = 1, 2
		call dgemv('n', mx*my, mx*my, 1.d0, c1, mx*my, q_col(k,:), 1, 1.d0, b(k,:), 1)
	enddo
	
	!-- a) solve a1*q_intm = c1*q_n + b, store result in b
	call dgtsv(mx*my, 1, a1_ld, a1_main, a1_ud, b(1,:), mx*my, info)
	call dgtsv(mx*my, 1, a1_ld_v, a1_main_v, a1_ud_v, b(2,:), mx*my, info)
	
	!-- calculate c2*q_intm + b2 and store result in b2
	do k = 1, 2
		call dgemv('n', mx*my, mx*my, 1.d0, c2, mx*my, b(k,:), 1, 1.d0, b2(k,:), 1)
	enddo

	!-- solve a2*q_star = c2*q_intm + b and store result in b2
	call dpbsv('u', mx*my, mx, 1, a2_bs, mx + 1, b2(1,:), mx*my, info)
	call dpbsv('u', mx*my, mx, 1, a2_bs_v, mx + 1, b2(2,:), mx*my, info)
	
	!-- store velocities updated for diffusion in q
	do k = 1, 2
		do j = 1, my
			q(k, 1:mx, j) = b2(k,((j-1)*mx+1):(j*mx))
		enddo
	enddo	

	!---------------------------------------------------------
    !	Step 4: Update pressure distribution by solving:
	!
	!			A*p = beta*bu*u_star + beta*bv*v_star
    !---------------------------------------------------------

	!-- calculate beta*bu*u_star and store result in bp
	call dgemv('n', mx*my, mx*my, beta, bu, mx*my, b2(1,:), 1, 0.d0, bp, 1)
	
	!-- calculate beta*bu*u_star + beta*bv*v_starand store the result in bp
	call dgemv('n', mx*my, mx*my, beta, bv, mx*my, b2(2,:), 1, 1.d0, bp, 1)

	!-- solve a*p = bp = beta*bu*u_star + beta*bv*v_star
	call dpbsv('u', mx*my, mx, 1, a_bs, mx + 1, bp, mx*my, info)

	!-- store the pressure distribution in q(3,:,:)
	do j = 1, my
		q(3,1:mx, j) = bp(((j-1)*mx+1):(j*mx))	
	enddo
	
    !-- Update ghost cells (required to calculate staggered edge velocities)
	call bc2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,[0,0,0,0])
	
	!---------------------------------------------------------
    !	Step 3: Update edge velocities for diffusion
    !---------------------------------------------------------
	!-- calculate x-direction edge velocities
    do j = 1-mbc,my+1
        do i = 1-mbc,mx
			u_edge(i,j) = (q(1,i+1,j) + q(1,i,j))/2.d0
		enddo
	enddo
	
	!-- calculate y-direction edge velocities
	do j = 1-mbc,my
        do i = 1-mbc,mx+1
			v_edge(i,j) = (q(2,i,j+1) + q(2,i,j))/2.d0
        enddo
	enddo
	
	!---------------------------------------------------------
    !	Step 5: Update edge velocities for pressure-driven flow
    !---------------------------------------------------------
	!-- x-direction edge velocities
	do j = 1-mbc, mx+mbc
		do i = 1-mbc,mx
			u_edge(i,j) = u_edge(i,j) - dt/dx*(q(3,i+1,j) - q(3,i,j))
		enddo
	enddo
	
	!-- y-direction edge velocities
	do j = 1-mbc, mx
		do i = 1-mbc,mx+mbc
			v_edge(i,j) = v_edge(i,j) - dt/dx*(q(3,i,j+1) - q(3,i,j))
		enddo
	enddo
	
	!-- Calculate the divergence in each cell of physical domain (excluding ghost cells)
	do j = 1,mx 
		do i = 1,mx
			q(4,i,j) = (u_edge(i,j) - u_edge(i-1,j) + v_edge(i,j) - v_edge(i,j-1))/dx
		enddo
	enddo

	!---------------------------------------------------------
    !	Step 6: Update the cell-centered velocities for 
	!			pressure-driven flow
    !---------------------------------------------------------
	do j = 1,mx
		do i = 1, mx
			q(1,i,j) =  q(1,i,j) - dt/dx/2.d0*(q(3,i+1,j) - q(3,i-1,j))
			q(2,i,j) = q(2,i,j) - dt/dx/2.d0*(q(3,i,j+1) - q(3,i,j-1))
		enddo
	enddo
	
end subroutine src2