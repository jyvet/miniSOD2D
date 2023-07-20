module mod_mesh

    use mod_constants
    use mod_nvtx
	use mod_maths

    implicit none

	contains

		subroutine create_mesh(nelem,npoin,xyzBase,connec,xyz)

			implicit none
			integer(4), intent(in)  :: nelem, npoin
			real(rp),   intent(in)  :: xyzBase(nnode,ndime)
			integer(4), intent(out) :: connec(nelem,nnode)
			real(rp),   intent(out) :: xyz(npoin,ndime)

			call nvtxStartRange("create_mesh")
			call gen_connectivity(nelem,connec)
			call gen_coordinates(nelem,npoin,connec,xyzBase,xyz)
			call nvtxEndRange

		end subroutine create_mesh

		subroutine gen_connectivity(nelem,connec)

			implicit none
			integer(4), intent(in)  :: nelem
			integer(4), intent(out) :: connec(nelem,nnode)
			integer(4)              :: ielem, inode

			call nvtxStartRange("gen_connectivity")
			do ielem = 1, nelem
				do inode = 1,nnode
					connec(ielem,inode) = (ielem-1)*nnode + inode
				end do
			end do
			call nvtxEndRange

		end subroutine gen_connectivity

		subroutine gen_coordinates(nelem,npoin,connec,xyzBase,xyz)

			implicit none
			integer(4), intent(in)    :: nelem, npoin, connec(nelem,nnode)
			real(rp),   intent(in)    :: xyzBase(nnode,ndime)
			real(rp),   intent(out)   :: xyz(npoin,ndime)
			integer(4)                :: ielem, inode, idime, ipoin
			real(rp)                  :: xyz0

			xyz0 = 0.0_rp

			call nvtxStartRange("gen_coordinates")
			do ielem = 1, nelem
				do inode = 1,nnode
					xyz(connec(ielem,inode),1) = xyzBase(inode,1)+xyz0
					xyz(connec(ielem,inode),2) = xyzBase(inode,2)+xyz0
					xyz(connec(ielem,inode),3) = xyzBase(inode,3)+xyz0
				end do
				xyz0 = xyz0+3.0_rp
			end do
			call nvtxEndRange

		end subroutine gen_coordinates

end module mod_mesh