subroutine adi_matrices(a1, c1, a2, a2_bs, c2, b, b2 ,a1_ld ,a1_main ,a1_ud ,mx,my ,dx ,dy ,dt)

	!------------------------------------------------------------------
	!
	!	Written by: Nathan Alexander (May 2024)
	!
	!   The purpose of this subroutine is to define the matrices used
	!   	in the calling program to update the velocity distribution.
	!		for diffusion. The dimensionally split diffusion equation 
	!		is discretized using central differences and recast in the 
	!		following matrix form:
	!
	!			a1*q_intm = c1*q_n + b
	!
	!		and
	!		
	!			a2*q_star = c2*q_intm + b
	!	
	!	input/output:
	!			a1		- see above
	!			c1 		- see above
	!			a2_bs	- the matrix a2 in banded storage
	!			c2		- see above
	!			b		- see above
	!			b2		- temporarily stores b
	!			a1_ld	- lower diagonal of matrix a1
	!			a1_main	- main diagonal of matrix a1
	!			a1_ud	- upper diagonal of matrix a1
	!			
	!------------------------------------------------------------------

	!-- allocate memory
	real(kind=8) :: a1(mx*my,mx*my), a2(mx*my,mx*my), c1(mx*my,mx*my), c2(mx*my,mx*my)
	real(kind=8) :: a1_ld(mx*my-1), a1_main(mx*my), a1_ud(mx*my-1), c1_main(mx*my)
	real(kind=8) :: a2_main(mx*my), c2_main(mx*my), c2_od(mx*my-1), a2_bs(mx+1,mx*my)
	real(kind=8) :: alpha, b(2, mx*my), b2(2, mx*my)
	real(kind=8), intent(in) :: dx,dy,dt
	
	integer :: i
	
	real(kind=8) :: Re	
	common /cparam/ Re 
		
	!-- parameters
	alpha = dt/(dx**2)/Re
	
	!-- initialize arrays
	a1_ld = 0.d0
	a1_ud = 0.d0
	c2_od = 0.d0
	a1 = 0.d0 
	a2 = 0.d0
	c1 = 0.d0
	c2 = 0.d0
	a1_main = 0.d0
	c2_main = 0.d0
	b = 0.d0
	b2 = 0.d0
	
	!---------------------
	!-- Construct b vector, (imposes no-slip and no-flow bcs)
	!---------------------
	do i = 1, mx*my
		if (i .gt. mx*my - mx) then
			b(1, i) = alpha              
			b2(1, i) = alpha       
		else
			b(1, i) = 0.d0
			b2(1, i) = 0.d0
		endif
	enddo
	
	!------------------------------------------------
	!-- define outer-diagonals for matrices a1 and c2
	!------------------------------------------------
	do i = 1, mx*my-1                    ! mod(a,p) division of a by p
		if (mod(i,mx) .eq. 0.d0) then
			a1_ld(i) = 0.d0
			a1_ud(i) = 0.d0
			c2_od(i) = 0.d0
		else
			a1_ld(i) = -alpha/2.d0
			a1_ud(i) = -alpha/2.d0
			c2_od(i) = alpha/2.d0
		endif
	enddo
	
	!----------------------
	!-- Construct a1 matrix
	!----------------------
	! edges for a1_main
	a1_main = 1.d0 + 3.d0/2.d0*alpha
	
	! center points for a1_main
	do i = 0, mx-1
		a1_main((mx+2+(i-1)*mx):(2*mx + (i-1)*mx -1)) = 1.d0 + alpha
	enddo
	
	! upper diagonal
	a1(1,2) = a1_ud(1)
	
	! main diagonal
	a1(1,1) = a1_main(1)
	a1(mx*my, mx*my) = a1_main(mx*my)
	
	! lower diagonal
	a1(2,1) = a1_ld(1)
	
	do i = 2, mx*my-1
		a1(i, i+1) = a1_ud(i)   ! upper diagonal
		a1(i, i) = a1_main(i)   ! main diagonal
		a1(i+1, i) = a1_ld(i)	! lower diagonal
	enddo
	
	!----------------------
	!-- Construct c1 matrix
	!----------------------
	!-- center points for c1
	c1_main = 1.d0 - alpha
	
	!-- edges for c1
	c1_main(1:mx) = 1.d0 - 3.d0/2.d0*alpha
	c1_main((mx*my-mx+1):mx*my) = 1.d0 - 3.d0/2.d0*alpha
	
	! define main diagonal for c1
	do i = 1, mx*my
		c1(i, i) = c1_main(i)	! main diagonal
	enddo
	
	! efine outer diagonals for c1
	do i = 1, mx*my-mx
		c1(i, i+mx) = alpha/2.d0	! upper diagonal
		c1(i+mx,i) = alpha/2.d0 	! lower diagonal
	enddo

	!----------------------
	!-- Construct a2 matrix
	!----------------------
	!-- center points for a2
	a2_main = 1.d0 + alpha
	
	!-- edges of physical domain for a2
	a2_main(1:mx) = 1.d0 + 3.d0/2.d0*alpha
	a2_main((mx*my-mx+1):mx*my) = 1.d0 + 3.d0/2.d0*alpha
	
	! main diagonal
	do i = 1, mx*my
		a2(i, i) = a2_main(i)
	enddo
	
	! outer diagonals
	do i = 1, mx*my-mx
		a2(i, i+mx) = -alpha/2.d0
		a2(i+mx,i) = -alpha/2.d0 
	enddo
	
	!----------------------
	!-- Construct c2 matrix
	!----------------------
	!-- edges for c2
	c2_main = 1.d0 - 3.d0/2.d0*alpha
	
	!-- center points for c2
	do i = 0, mx-1
		c2_main((mx+2+(i-1)*mx):(2*mx + (i-1)*mx -1)) = 1.d0 - alpha
	enddo
	
	! upper diagonal
	c2(1,2) = c2_od(1)
	
	! main diagonal
	c2(1,1) = c2_main(1)
	c2(mx*my, mx*my) = c2_main(mx*my)
	
	! lower diagonal
	c2(2,1) = c2_od(1)
	
	do i = 2, mx*my-1
		c2(i, i+1) = c2_od(i)   ! upper diagonal
		c2(i, i) = c2_main(i)   ! main diagonal
		c2(i+1, i) = c2_od(i)	! lower diagonal
	enddo
	
	do i = 1, mx*my
		do j = 1, mx*my
			if (i .ge. max(1,j-mx) .and. i .le. j) then
				a2_bs(mx+1+i-j,j) = a2(i,j)
			endif
		enddo
	enddo
	
end subroutine adi_matrices