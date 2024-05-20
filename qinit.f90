!
!
!
!=====================================================
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
!=====================================================
!
!     # Set initial conditions for q.
	implicit none
	integer :: i, j
	real (kind=8) :: xcell, ycell
	
	integer, intent(in) :: meqn,mbc,mx,my,maux
	real(kind=8), intent(in) :: xlower,dx,ylower,dy
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
	
	!-- fluid at rest
	q = 0.
	
return
end
