subroutine poisson_matrices(bu, bv, a_bs, mx, my, dx, dy, dt)

	!------------------------------------------------------------------
	!
	!	Written by: Nathan Alexander (May 2024)
	!
	!   The purpose of this subroutine is to define the matrices used
	!   	in the calling program to update the pressure distribution.
	!		The distribution is updated by the calling program by solving
	!
	!			p_xx + p_yy + c_0*p =  1/dt*(u_x + v_y)
	!
	!		which is discretized using central differences and 
	!		recast in the following matrix form:
	!
	!			a*p = beta*bu*u + beta*bu*v,
	!
	!		where beta = dx/dt/2.
	!	
	!	input/output:
	!			bu - u derivative matrix on r.h.s. of Poisson equation 
	!			bv - v derivative matrix on r.h.s. of Poisson equation
	!			a_bs - the matrix a in banded storage
	!
	!------------------------------------------------------------------
	
	!-- allocate memory
	real(kind=8) :: bu(mx*my,mx*my), bv(mx*my,mx*my)
	real(kind=8) :: bu_ud(mx*my-1), bu_ld(mx*my-1), bv_ud(mx*my-mx), bv_ld(mx*my-mx), bu_main(mx*my), bv_main(mx*my)
	
	real(kind=8) :: a_bs(mx+1,mx*my)
	real(kind=8) :: a_uud(mx*(mx-1)), a_ud(mx*my-1), a_main(mx*my), a_ld(mx*my-1), a_lld(mx*(mx-1))
	
	real(kind=8), intent(in) :: dx,dy,dt
	real(kind=8) :: c_0


	!-- define parameters
	c_0 = 0.000000001d0					! proportionality constant for fictiticous pressure source  
	
	!-- initialize arrays
	a_bs = 0.d0							! matrix A in banded storage	
	a_main = 3.d0 + c_0*dx**2.d0		! main diagonal for poisson matrix A
	a_uud = -1.d0						! mx elements above main diagonal
	a_ud = -1.d0						! 1 elements above main diagonal
	a_ld = -1.d0						! 1 element below main diagonal
	a_lld = -1.d0						! mx elements below main diagonal
	
	bu = 0.d0										
	bv = 0.d0										
	
	bu_ud = -1.d0						! upper diag. of bu matrix
	bu_ld = 1.d0						! lower diag. of bu matrix
	bv_ud = -1.d0						! upper diag. of bv matrix			
	bv_ld = 1.d0						! lower diag. of bv matrix
	bu_main = 0.d0						! main diag. of bu matrix
	bv_main = 0.d0						! main diag. of bv matrix
	
	!----------------------
	!-- Construct bu matrix
	!----------------------
	!-- bu outer diagonals
	do i = 1, mx*my-1
		if (mod(i,mx) == 0) then
			bu_ud(i) = 0.d0
			bu_ld(i) = 0.d0
		endif
	enddo
	
	!-- bu main diagonal
	do i = 1, mx*my
		if (mod(i,mx) == 0) then
			bu_main(i) = 1.d0
		else if (mod(i-1,mx) == 0) then
			bu_main(i) = -1.d0
		endif
		bu(i,i) = bu_main(i)
	enddo
	
	!-- bu outer diagonals
	do j = 1, mx*my-1
		bu(j,j+1) = bu_ud(j)
		bu(j+1,j) = bu_ld(j)
	enddo
	
	!----------------------
	!-- Construct bv matrix
	!----------------------
	!-- bv main diagonal
	bv_main(1:mx) = -1.d0
	bv_main(mx*my-mx+1:mx*my) = 1.d0
	
	!-- populate main diagonal of bv
	do i = 1,mx*my
		bv(i,i) = bv_main(i)
	enddo

	!-- bv outer diagonals
	do j = 1, mx*my-mx
		bv(j,j+mx) = bv_ud(j)
		bv(j+mx,j) = bv_ld(j)
	enddo
	
	!---------------------
	!-- Construct matrix A
	!---------------------
	
	!-- determina diagonals of matrix A
	!-- corners
	a_main(1) = 2.d0 + c_0*dx**2.d0
	a_main(mx) = 2.d0 + c_0*dx**2.d0
	a_main(mx*my-(mx-1)) = 2.d0 + c_0*dx**2.d0
	a_main(mx*my) = 2.d0 + c_0*dx**2.d0
	
	!-- center points
	do i = 1, mx-2
		a_main((mx+2+(i-1)*mx):(2*mx + (i-1)*mx -1)) = 4.d0 + c_0*dx**2.d0
	enddo
	
	do i = 1, mx-1
		if (i < mx) then
			a_ud(i*mx) = 0.d0
			a_ld(i*mx) = 0.d0
		endif
	enddo

	!-- construct banded storage array for the matrix A	
	a_bs(1,(mx+1):) = a_uud(:)
	a_bs(mx,2:) = a_ud(:)
	a_bs(mx+1,:) = a_main(:)
	
end subroutine poisson_matrices