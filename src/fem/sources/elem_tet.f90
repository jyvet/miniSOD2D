module elem_tet

use mod_constants

implicit none

contains

    subroutine tet_04(xi,eta,zeta,N,dN)
    
        implicit none
        real(rp), intent(in)  :: xi,eta,zeta
        real(rp), intent(out) :: N(4),dN(3,4)

        N(1) = 1.0_rp-xi-eta-zeta
        N(2) = xi
        N(3) = eta
        N(4) = zeta

        dN(1,1) = -1.0_rp
        dN(2,1) = -1.0_rp
        dN(3,1) = -1.0_rp

        dN(1,2) = 1.0_rp
        dN(2,2) = 0.0_rp
        dN(3,2) = 0.0_rp

        dN(1,3) = 0.0_rp
        dN(2,3) = 1.0_rp
        dN(3,3) = 0.0_rp

        dN(1,4) = 0.0_rp
        dN(2,4) = 0.0_rp
        dN(3,4) = 1.0_rp

    end subroutine tet_04

end module elem_tet