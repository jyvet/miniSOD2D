module elem_qua

use mod_constants
use mod_maths

	implicit none
	!integer(4), parameter :: quad_order_edges(4,2) = transpose(reshape([1,2,2,3,3,4,4,1],(/2,4/)))
	integer(4), allocatable :: quad_order_edges(:,:)

	contains

		subroutine init_basic_qua()
			implicit none
			allocate(quad_order_edges(4,2))
			quad_order_edges = transpose(reshape([1,2,2,3,3,4,4,1],(/2,4/)))
		end subroutine init_basic_qua

		subroutine quad_highorder(mporder,mnpbou,xi,eta,atoIJ,N,dN) ! QUA16 element
			implicit none
						integer(4),intent(in) :: mporder,mnpbou
			real(rp),intent(in)   :: xi,eta
			integer(4),intent(in) :: atoIJ(mnpbou)
			real(rp),intent(out)  :: N(mnpbou), dN(2,mnpbou)
			real(rp)              :: xi_grid(mporder+1)

			call getGaussLobattoLegendre_roots(mporder,xi_grid)
			call DoubleTensorProduct(mporder,mnpbou,xi_grid,xi,eta,atoIJ,N,dN)
		end subroutine quad_highorder

		subroutine quad_edges(ielem,nelem,npoin,connec,coord,ncorner,nedge,dist)

			implicit none

			integer(4), intent(in)            :: ielem, nelem, npoin
			integer(4), intent(in)            :: connec(nelem,nnode)
			real(rp),    intent(in)            :: coord(npoin,ndime)
			integer(4), intent(out)           :: ncorner, nedge
			real(rp),    intent(out)           :: dist(4,ndime)
			integer(4)                        :: ind(nnode)
			real(rp)                           :: xp(4,ndime)

			ind = connec(ielem,:)
			ncorner = 4
			nedge = 4

			xp(1:4,1:ndime) = coord(ind(1:4),1:ndime) ! Corner coordinates
			dist(1,:) = xp(2,:)-xp(1,:)
			dist(2,:) = xp(3,:)-xp(2,:)
			dist(3,:) = xp(4,:)-xp(3,:)
			dist(4,:) = xp(1,:)-xp(4,:)

		end subroutine quad_edges

end module
