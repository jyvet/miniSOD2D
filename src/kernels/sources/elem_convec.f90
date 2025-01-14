module elem_convec

   use mod_nvtx
   use mod_constants

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Computes convective term for Euler/NS equation system, as well             !
      ! as for any generic scalar transport that might occur. Based                !
      ! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
      ! This module can be passed to CUDA in order to do fine-grained parallelism. !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains
         subroutine full_convec_ijk_momentum(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,Rmom,iteration)
            implicit none

            integer(4), intent(in)  :: nelem, npoin, iteration
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime),rho(npoin),pr(npoin)!!,qr,ur !! Reordered matrices for regular memory accesses
            real(rp),    intent(out) :: Rmom(npoin,ndime)
            integer(4)              :: ielem,igaus,idime,jdime,kdime,inode,isoI,isoJ,isoK,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime),gradIsoU(ndime,ndime),gradIsoF(ndime,ndime,ndime),gradIsoQ(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradP(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ,gradQ(ndime,ndime),gradIsoFuu(ndime,ndime,ndime),gradRE(ndime),divFk
            real(rp)                 :: ul(nnode,ndime),ql(nnode,ndime),rhol(nnode),prl(nnode),fl(nnode,ndime,ndime),fuul(nnode,ndime,ndime),divFuu(ndime)
            real(rp), dimension(porder+1) :: dlxi_ip,dleta_ip,dlzeta_ip

            !$acc kernels
            Rmom(:,:) = 0.0_rp
            !$acc end kernels
            if (iteration>1) then
               !$acc parallel loop gang private(ipoin,idime,jdime,kdime,Re_mom,ul,ql,rhol,prl,fl,fuul) !!vector_length(vecLength)
               do ielem = 1,nelem
                  !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoRho,gradIsoP,gradIsoU,gradIsoF,gradIsoFuu,gradIsoQ,gradRho,gradP,gradRE,gradU,divF,divU,divQ,gradQ,divFuu,isoI,isoJ,isoK,ii)
                  do igaus = 1,ngaus
                     !$acc loop seq
                     do ii=1,porder+1
                        dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                        dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                        dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                     end do
                     isoI = gmshAtoI(igaus)
                     isoJ = gmshAtoJ(igaus)
                     isoK = gmshAtoK(igaus)

                     gradIsoRho(:) = 0.0_rp
                     gradIsoP(:) = 0.0_rp
                     gradIsoU(:,:) = 0.0_rp
                     gradIsoF(:,:,:) = 0.0_rp
                     gradIsoFuu(:,:,:) = 0.0_rp
                     gradIsoQ(:,:) = 0.0_rp

                     !$acc loop seq
                     do ii=1,porder+1
                        gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rho(connec(ielem,invAtoIJK(ii,isoJ,isoK)))
                        gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rho(connec(ielem,invAtoIJK(isoI,ii,isoK)))
                        gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rho(connec(ielem,invAtoIJK(isoI,isoJ,ii)))

                        gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*pr(connec(ielem,invAtoIJK(ii,isoJ,isoK)))
                        gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*pr(connec(ielem,invAtoIJK(isoI,ii,isoK)))
                        gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*pr(connec(ielem,invAtoIJK(isoI,isoJ,ii)))

                        !$acc loop seq
                        do idime=1,ndime
                           gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)
                           gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)
                           gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)

                           gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*q(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)
                           gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)
                           gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)

                           do jdime=1, ndime
                              gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*q(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),jdime)
                              gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),jdime)
                              gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),jdime)
                              gradIsoFuu(idime,jdime,1) = gradIsoFuu(idime,jdime,1) + dlxi_ip(ii)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),jdime)
                              gradIsoFuu(idime,jdime,2) = gradIsoFuu(idime,jdime,2) + dleta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),jdime)
                              gradIsoFuu(idime,jdime,3) = gradIsoFuu(idime,jdime,3) + dlzeta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),jdime)
                           end do
                        end do
                     end do

                     gradRho(:) = 0.0_rp
                     gradP(:) = 0.0_rp
                     gradU(:,:) = 0.0_rp
                     gradQ(:,:) = 0.0_rp
                     divF(:) = 0.0_rp
                     divFuu(:) = 0.0_rp
                     divQ = 0.0_rp

                     !$acc loop seq
                     do idime=1, ndime
                        !$acc loop seq
                        do jdime=1, ndime
                           gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                           gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                           !$acc loop seq
                           do kdime=1,ndime
                              gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                              gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                              divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                              divFuu(idime) = divFuu(idime) + He(jdime,kdime,igaus,ielem)*gradIsoFuu(idime,jdime,kdime)
                           end do
                        end do
                     end do
                     divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                     divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3)

                     !$acc loop seq
                     do idime=1, ndime
                        Re_mom(igaus,idime) = 0.25_rp*(divF(idime)+u(connec(ielem,igaus),idime)*divQ+q(connec(ielem,igaus),idime)*divU+rho(connec(ielem,igaus))*divFuu(idime)) + gradP(idime)
                        do jdime=1, ndime
                           Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.25_rp*(u(connec(ielem,igaus),idime)*u(connec(ielem,igaus),jdime)*gradRho(jdime)  &
                                                                     + q(connec(ielem,igaus),idime)*gradU(idime,jdime)+u(connec(ielem,igaus),jdime)*gradQ(idime,jdime))
                        end do
                        Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                     end do
                  end do
                  !
                  ! Final assembly
                  !
                  !$acc loop vector collapse(2)
                  do idime = 1,ndime
                     do inode = 1,nnode
                        !$acc atomic update
                        Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                        !$acc end atomic
                     end do
                  end do
               end do
               !$acc end parallel loop 
            else 
               !$acc parallel loop gang private(ipoin,idime,jdime,kdime,Re_mom,ul,ql,rhol,prl,fl,fuul) !!vector_length(vecLength)
               do ielem = 1,nelem
                  !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoRho,gradIsoP,gradIsoU,gradIsoF,gradIsoFuu,gradIsoQ,gradRho,gradP,gradRE,gradU,divF,divU,divQ,gradQ,divFuu,isoI,isoJ,isoK,ii)
                  do igaus = 1,ngaus
                     !$acc loop seq
                     do ii=1,porder+1
                        dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                        dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                        dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                     end do
                     isoI = gmshAtoI(igaus)
                     isoJ = gmshAtoJ(igaus)
                     isoK = gmshAtoK(igaus)

                     gradIsoRho(:) = 0.0_rp
                     gradIsoP(:) = 0.0_rp
                     gradIsoU(:,:) = 0.0_rp
                     gradIsoF(:,:,:) = 0.0_rp
                     gradIsoFuu(:,:,:) = 0.0_rp
                     gradIsoQ(:,:) = 0.0_rp

                     !$acc loop seq
                     do ii=1,porder+1
                        gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rho(connec(ielem,invAtoIJK(ii,isoJ,isoK)))
                        gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rho(connec(ielem,invAtoIJK(isoI,ii,isoK)))
                        gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rho(connec(ielem,invAtoIJK(isoI,isoJ,ii)))

                        gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*pr(connec(ielem,invAtoIJK(ii,isoJ,isoK)))
                        gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*pr(connec(ielem,invAtoIJK(isoI,ii,isoK)))
                        gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*pr(connec(ielem,invAtoIJK(isoI,isoJ,ii)))

                        !$acc loop seq
                        do idime=1,ndime
                           gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)
                           gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)
                           gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)

                           gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*q(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)
                           gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)
                           gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)

                           do jdime=1, ndime
                              gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*q(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),jdime)
                              gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),jdime)
                              gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*q(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),jdime)
                              gradIsoFuu(idime,jdime,1) = gradIsoFuu(idime,jdime,1) + dlxi_ip(ii)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),idime)*u(connec(ielem,invAtoIJK(ii,isoJ,isoK)),jdime)
                              gradIsoFuu(idime,jdime,2) = gradIsoFuu(idime,jdime,2) + dleta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),idime)*u(connec(ielem,invAtoIJK(isoI,ii,isoK)),jdime)
                              gradIsoFuu(idime,jdime,3) = gradIsoFuu(idime,jdime,3) + dlzeta_ip(ii)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),idime)*u(connec(ielem,invAtoIJK(isoI,isoJ,ii)),jdime)
                           end do
                        end do
                     end do

                     gradRho(:) = 0.0_rp
                     gradP(:) = 0.0_rp
                     gradU(:,:) = 0.0_rp
                     gradQ(:,:) = 0.0_rp
                     divF(:) = 0.0_rp
                     divFuu(:) = 0.0_rp
                     divQ = 0.0_rp

                     !$acc loop seq
                     do idime=1, ndime
                        !$acc loop seq
                        do jdime=1, ndime
                           gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                           gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
                           !$acc loop seq
                           do kdime=1,ndime
                              gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                              gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                              divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
                              divFuu(idime) = divFuu(idime) + He(jdime,kdime,igaus,ielem)*gradIsoFuu(idime,jdime,kdime)
                           end do
                        end do
                     end do
                     divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                     divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3)

                     !$acc loop seq
                     do idime=1, ndime
                        Re_mom(igaus,idime) = 0.25_rp*(divF(idime)+u(connec(ielem,igaus),idime)*divQ+q(connec(ielem,igaus),idime)*divU+rho(connec(ielem,igaus))*divFuu(idime)) + gradP(idime)
                        do jdime=1, ndime
                           Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.25_rp*(u(connec(ielem,igaus),idime)*u(connec(ielem,igaus),jdime)*gradRho(jdime)  &
                                                                     + q(connec(ielem,igaus),idime)*gradU(idime,jdime)+u(connec(ielem,igaus),jdime)*gradQ(idime,jdime))
                        end do
                        Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
                     end do
                  end do
                  !
                  ! Final assembly
                  !
                  !$acc loop vector collapse(2)
                  do idime = 1,ndime
                     do inode = 1,nnode
                        !$acc atomic update
                        Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
                        !$acc end atomic
                     end do
                  end do
               end do
               !$acc end parallel loop           
            end if
         end subroutine full_convec_ijk_momentum

         subroutine full_convec_ijk_energy(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,E,Rener)
            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime),u(npoin,ndime),rho(npoin),E(npoin)
            real(rp),    intent(out) :: Rener(npoin)
            integer(4)              :: ielem,igaus,idime,jdime,kdime,inode,isoI,isoJ,isoK,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_ener(nnode),aux
            real(rp)                 :: gradIsoRho(ndime),gradIsoE(ndime),gradIsoU(ndime,ndime),gradIsoF(ndime,ndime,ndime),gradIsoQ(ndime,ndime),gradIsoFe(ndime,ndime),gradIsoRE(ndime),gradIsoFue(ndime,ndime),gradIsok(ndime),gradIsoRk(ndime),gradIsoFuk(ndime,ndime),gradIsoFk(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradE(ndime),gradU(ndime,ndime),divU,divFe,divQ,gradQ(ndime,ndime),gradIsoFuu(ndime,ndime,ndime),gradRE(ndime),divFue,gradk(ndime),divFk,gradRk(ndime),divFuk
            real(rp)                 :: ul(nnode,ndime),ql(nnode,ndime),rhol(nnode),El(nnode),fel(nnode,ndime),REl(nnode),fuel(nnode,ndime),kl(nnode),fkl(nnode,ndime),fukl(nnode,ndime),Rkl(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip,dleta_ip,dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rener(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,idime,jdime,kdime,Re_ener,ul,ql,rhol,El,REl,fel,fuel,kl,Rkl,fkl,fukl) !!vector_length(vecLength)
            do ielem = 1,nelem
               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
               end do
               !$acc loop vector
               do inode = 1,nnode
                  kl(inode) = 0.0_rp
                  !$acc loop seq
                  do idime = 1,ndime
                     kl(inode) = kl(inode) + u(ipoin(inode),idime)*u(ipoin(inode),idime)*0.5_rp
                  end do
               end do
               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime) = u(ipoin(inode),idime)
                     ql(inode,idime) = q(ipoin(inode),idime)
                     fel(inode,idime) = rho(ipoin(inode))*E(ipoin(inode))*u(ipoin(inode),idime)
                     fuel(inode,idime) = E(ipoin(inode))*u(ipoin(inode),idime)
                     fkl(inode,idime) = rho(ipoin(inode))*kl(inode)*u(ipoin(inode),idime)
                     fukl(inode,idime) = kl(inode)*u(ipoin(inode),idime)
                  end do
               end do
               !$acc loop vector
               do inode = 1,nnode
                  rhol(inode) = rho(ipoin(inode))
                  El(inode) = E(ipoin(inode))
                  REl(inode) = rho(ipoin(inode))*E(ipoin(inode))
                  Rkl(inode) = rho(ipoin(inode))*kl(inode)
               end do
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoRho,gradIsoRE,gradIsoRk,gradIsoE,gradIsok,gradIsoU,gradIsoF,gradIsoFuu,gradIsoQ,gradIsoFe,gradIsoFue,gradIsoFk,gradIsoFuk,gradRho,gradE,gradRE,gradk,gradRk,gradU,divU,divQ,divFe,divFue,divFk,divFuk,gradQ,isoI,isoJ,isoK,ii)
               do igaus = 1,ngaus
                  !$acc loop seq
                  do ii=1,porder+1
                     dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                     dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                     dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                  end do
                  isoI = gmshAtoI(igaus)
                  isoJ = gmshAtoJ(igaus)
                  isoK = gmshAtoK(igaus)

                  gradIsoRho(:) = 0.0_rp
                  gradIsoE(:) = 0.0_rp
                  gradIsoRE(:) = 0.0_rp
                  gradIsok(:) = 0.0_rp
                  gradIsoRk(:) = 0.0_rp
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0.0_rp
                  gradIsoFe(:,:) = 0.0_rp
                  gradIsoFue(:,:) = 0.0_rp
                  gradIsoFk(:,:) = 0.0_rp
                  gradIsoFuk(:,:) = 0.0_rp
                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                     gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                     gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                     gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(El(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(El(invAtoIJK(isoI,ii,isoK)))
                     gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(El(invAtoIJK(isoI,isoJ,ii)))

                     gradIsoRE(1) = gradIsoRE(1) + dlxi_ip(ii)*(REl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoRE(2) = gradIsoRE(2) + dleta_ip(ii)*(REl(invAtoIJK(isoI,ii,isoK)))
                     gradIsoRE(3) = gradIsoRE(3) + dlzeta_ip(ii)*(REl(invAtoIJK(isoI,isoJ,ii)))

                     gradIsok(1) = gradIsok(1) + dlxi_ip(ii)*(kl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsok(2) = gradIsok(2) + dleta_ip(ii)*(kl(invAtoIJK(isoI,ii,isoK)))
                     gradIsok(3) = gradIsok(3) + dlzeta_ip(ii)*(kl(invAtoIJK(isoI,isoJ,ii)))

                     gradIsoRk(1) = gradIsoRk(1) + dlxi_ip(ii)*(Rkl(invAtoIJK(ii,isoJ,isoK)))
                     gradIsoRk(2) = gradIsoRk(2) + dleta_ip(ii)*(Rkl(invAtoIJK(isoI,ii,isoK)))
                     gradIsoRk(3) = gradIsoRk(3) + dlzeta_ip(ii)*(Rkl(invAtoIJK(isoI,isoJ,ii)))

                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFue(idime,1) = gradIsoFue(idime,1) + dlxi_ip(ii)*fuel(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFue(idime,2) = gradIsoFue(idime,2) + dleta_ip(ii)*fuel(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFue(idime,3) = gradIsoFue(idime,3) + dlzeta_ip(ii)*fuel(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFk(idime,1) = gradIsoFk(idime,1) + dlxi_ip(ii)*fkl(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFk(idime,2) = gradIsoFk(idime,2) + dleta_ip(ii)*fkl(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFk(idime,3) = gradIsoFk(idime,3) + dlzeta_ip(ii)*fkl(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoFuk(idime,1) = gradIsoFuk(idime,1) + dlxi_ip(ii)*fukl(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoFuk(idime,2) = gradIsoFuk(idime,2) + dleta_ip(ii)*fukl(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoFuk(idime,3) = gradIsoFuk(idime,3) + dlzeta_ip(ii)*fukl(invAtoIJK(isoI,isoJ,ii),idime)
                     end do
                  end do

                  gradRho(:) = 0.0_rp
                  gradE(:) = 0.0_rp
                  gradRE(:) = 0.0_rp
                  gradk(:) = 0.0_rp
                  gradRk(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  gradQ(:,:) = 0.0_rp
                  divQ = 0.0_rp
                  divFe = 0.0_rp
                  divFue = 0.0_rp
                  divFk = 0.0_rp
                  divFuk = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                         gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
                         gradRE(idime)  = gradRE(idime) + He(idime,jdime,igaus,ielem) * gradIsoRE(jdime)
                         divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
                         divFue = divFue + He(idime,jdime,igaus,ielem) * gradIsoFue(idime,jdime)
                         gradk(idime)   = gradk(idime) + He(idime,jdime,igaus,ielem) * gradIsok(jdime)
                         gradRk(idime)  = gradRk(idime) + He(idime,jdime,igaus,ielem) * gradIsoRk(jdime)
                         divFk = divFk + He(idime,jdime,igaus,ielem) * gradIsoFk(idime,jdime)
                         divFuk = divFuk + He(idime,jdime,igaus,ielem) * gradIsoFuk(idime,jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                         end do
                     end do
                  end do
                  divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                  divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3)
                  Re_ener(igaus) = 0.25_rp*(divFe+El(igaus)*divQ+REl(igaus)*divU+rhol(igaus)*divFue) + &
                                   0.25_rp*(divFk+kl(igaus)*divQ+Rkl(igaus)*divU+rhol(igaus)*divFuk)
                  !$acc loop seq
                  do idime=1, ndime
                     Re_ener(igaus) = Re_ener(igaus) + 0.25_rp*(ql(igaus,idime)*gradE(idime)+ul(igaus,idime)*gradRE(idime)+&
                                                                fuel(igaus,idime)*gradRho(idime)) +&
                                                     0.25_rp*(ql(igaus,idime)*gradk(idime)+ul(igaus,idime)*gradRk(idime)+&
                                                                fukl(igaus,idime)*gradRho(idime))
                  end do
                  Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
               end do
               !
               ! Final assembly
               !
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rener(ipoin(inode)) = Rener(ipoin(inode))+Re_ener(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
         end subroutine full_convec_ijk_energy

         subroutine full_convec_ijk_mass(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,Rmass)
            implicit none

            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode),gmshAtoJ(nnode),gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime),u(npoin,ndime),rho(npoin)
            real(rp),    intent(out) :: Rmass(npoin)
            integer(4)              :: ielem,igaus,idime,jdime,kdime,inode,isoI,isoJ,isoK,ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mass(nnode)
            real(rp)                 :: gradIsoRho(ndime),gradIsoU(ndime,ndime),gradIsoQ(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradU(ndime,ndime),divU,divQ,gradQ(ndime,ndime)
            real(rp)                 :: ul(nnode,ndime),ql(nnode,ndime),rhol(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip,dleta_ip,dlzeta_ip

            call nvtxStartRange("Full convection")
            !$acc kernels
            Rmass(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(ipoin,idime,jdime,kdime,Re_mass,ul,ql,rhol) !!vector_length(vecLength)
            do ielem = 1,nelem

               !$acc loop vector
               do inode = 1,nnode
                  ipoin(inode) = connec(ielem,inode)
                  rhol(inode) = rho(ipoin(inode))                 
               end do

               !$acc loop vector collapse(2)
               do idime = 1,ndime
                  do inode = 1,nnode
                     ul(inode,idime) = u(ipoin(inode),idime)
                     ql(inode,idime) = q(ipoin(inode),idime)
                  end do
               end do
 
               !$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip,gradIsoRho,gradIsoU,gradIsoQ,gradRho,gradU,divU,divQ,gradQ,isoI,isoJ,isoK,ii)
               do igaus = 1,ngaus
                  !$acc loop seq
                  do ii=1,porder+1
                     dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
                     dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
                     dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
                  end do
                  isoI = gmshAtoI(igaus)
                  isoJ = gmshAtoJ(igaus)
                  isoK = gmshAtoK(igaus)

                  gradIsoRho(:) = 0.0_rp
                  gradIsoU(:,:) = 0.0_rp
                  gradIsoQ(:,:) = 0.0_rp

                  !$acc loop seq
                  do ii=1,porder+1
                     gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                     gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
                     gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                     !$acc loop seq
                     do idime=1,ndime
                        gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

                        gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
                        gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
                        gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)
                     end do
                  end do

                  gradRho(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  gradQ(:,:) = 0.0_rp
                  divQ = 0.0_rp
                  !$acc loop seq
                  do idime=1, ndime
                     !$acc loop seq
                     do jdime=1, ndime
                         gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                         !$acc loop seq
                         do kdime=1,ndime
                            gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                            gradQ(idime,jdime) = gradQ(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoQ(idime,kdime)
                         end do
                     end do
                  end do
                  divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3)
                  divQ  = gradQ(1,1)  + gradQ(2,2)  + gradQ(3,3)
                  Re_mass(igaus) = 0.5_rp*(divQ+rhol(igaus)*divU)
                  !$acc loop seq
                  do idime=1, ndime
                     Re_mass(igaus) = Re_mass(igaus) + 0.5_rp*(ul(igaus,idime)*gradRho(idime))
                  end do
                  Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
               end do
               !
               ! Final assembly
               !
               !$acc loop vector
               do inode = 1,nnode
                  !$acc atomic update
                  Rmass(ipoin(inode)) = Rmass(ipoin(inode))+Re_mass(inode)
                  !$acc end atomic
               end do
            end do
            !$acc end parallel loop
         end subroutine full_convec_ijk_mass

         subroutine full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener,iteration)
            implicit none

            integer(4), intent(in)  :: nelem, npoin,iteration
            integer(4), intent(in)  :: connec(nelem,nnode)
            real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
            real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
            real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
            integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
            real(rp),    intent(out) :: Rmass(npoin)
            real(rp),    intent(out) :: Rmom(npoin,ndime)
            real(rp),    intent(out) :: Rener(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, kdime, inode, isoI, isoJ, isoK, ii
            integer(4)              :: ipoin(nnode)
            real(rp)                 :: Re_mom(nnode,ndime)
            real(rp)                 :: Re_mass(nnode), Re_ener(nnode),aux
            real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime),gradIsoRE(ndime),gradIsoFue(ndime,ndime), gradIsok(ndime),gradIsoRk(ndime),gradIsoFuk(ndime,ndime),gradIsoFk(ndime,ndime)
            real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ,gradQ(ndime,ndime),gradIsoFuu(ndime,ndime,ndime),gradRE(ndime),divFue,gradk(ndime),divFk,gradRk(ndime),divFuk
            real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime),fuul(nnode,ndime,ndime),divFuu(ndime),REl(nnode),fuel(nnode,ndime),kl(nnode),fkl(nnode,ndime),fukl(nnode,ndime),Rkl(nnode)
            real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

            call nvtxStartRange("Full convection")
            call full_convec_ijk_momentum(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,Rmom,iteration)
            call full_convec_ijk_energy(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,E,Rener)
            call full_convec_ijk_mass(nelem,npoin,connec,He,gpvol,dlxigp_ip,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,Rmass)
            call nvtxEndRange
         end subroutine full_convec_ijk

         subroutine fem_convec(nelem,npoin,connec,Ngp,dNgp,He,gpvol,u,q,rho,pr,E,Rmass,Rmom,Rener)

            implicit none
            integer(4), intent(in)  :: nelem, npoin
            integer(4), intent(in)  :: connec(nelem,4)
            real(rp),    intent(in)  :: Ngp(4,4), dNgp(ndime,4,4)
            real(rp),    intent(in)  :: He(ndime,ndime,4,nelem)
            real(rp),    intent(in)  :: gpvol(1,4,nelem)
            real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
            real(rp),    intent(out) :: Rmass(npoin)
            real(rp),    intent(out) :: Rmom(npoin,ndime)
            real(rp),    intent(out) :: Rener(npoin)
            integer(4)              :: ielem, igaus, idime, jdime, kdime, inode, jnode
            integer(4)              :: ipoin(4)
            real(rp)                 :: Re_mom(4,ndime)
            real(rp)                 :: Re_mass(4), Re_ener(4)
            real(rp)                 :: ul(4,ndime), ql(4,ndime), rhol(4), prl(4),El(4),fel(4,ndime),fl(4,ndime,ndime), fpl(4)
            real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ
            real(rp)                 :: gpcar(4,ndime), aux1, aux2, aux3, aux4

            call nvtxStartRange("FEM convection")
            !$acc kernels
            Rmom(:,:) = 0.0_rp
            Rmass(:) = 0.0_rp
            Rener(:) = 0.0_rp
            !$acc end kernels

            !$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel,ipoin) present(connec,u,q,rho,pr,E,Rmass,Rmom,Rener)
            do ielem = 1,nelem
               Re_mass(:) = 0.0_rp
               Re_ener(:) = 0.0_rp
               Re_mom(:,:) = 0.0_rp

               !$acc loop seq private(divQ, divFe, gradRho, gradP, gradE, gradU, divF, divU, gpcar, aux1, aux2, aux3, aux4)
               do igaus = 1,4
                  divQ = 0.0_rp
                  divFe = 0.0_rp
                  gradRho(:) = 0.0_rp
                  gradE(:) = 0.0_rp
                  gradP(:) = 0.0_rp
                  gradU(:,:) = 0.0_rp
                  divF(:) = 0.0_rp

                  !$acc loop seq
                  do idime = 1,ndime
                     !$acc loop seq
                     do inode = 1,4
                        gpcar(inode,idime) = DOT_PRODUCT(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                     end do
                  end do

                  !$acc loop seq
                  do idime = 1,ndime
                     !$acc loop seq
                     do jnode = 1,4
                        divQ = divQ + gpcar(jnode,idime)*q(connec(ielem,jnode),idime)
                        divFe = divFe + gpcar(jnode,idime)*(E(connec(ielem,jnode))+pr(connec(ielem,jnode)))*u(connec(ielem,jnode),idime)
                        gradRho(idime) = gradRho(idime) + gpcar(jnode,idime)*rho(connec(ielem,jnode))
                        gradE(idime) = gradE(idime) + gpcar(jnode,idime)*(E(connec(ielem,jnode))+pr(connec(ielem,jnode)))
                        gradP(idime) = gradP(idime) + gpcar(jnode,idime)*pr(connec(ielem,jnode))

                        !$acc loop seq
                        do jdime = 1,ndime
                           gradU(idime,jdime) = gradU(idime,jdime) + gpcar(jnode,jdime)*u(connec(ielem,jnode),idime)
                           divF(idime) = divF(idime) + gpcar(jnode,jdime)*q(connec(ielem,jnode),idime)*u(connec(ielem,jnode),jdime)
                        end do
                     end do
                  end do

                  divU = gradU(1,1) + gradU(2,2) + gradU(3,3)

                  !$acc loop seq
                  do inode = 1,4
                     aux1 = 0.0_rp
                     aux2 = 0.0_rp

                     !$acc loop seq
                     do idime = 1,ndime
                        aux1 = aux1 + gradRho(idime)*u(connec(ielem,inode),idime)
                        aux2 = aux2 + gradE(idime)*u(connec(ielem,inode),idime)
                     end do

                     Re_mass(inode) = Re_mass(inode) + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*0.5_rp*(divQ + rho(connec(ielem,inode))*divU + aux1)
                     Re_ener(inode) = Re_ener(inode) + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*0.5_rp*(divFe + (E(connec(ielem,inode))+pr(connec(ielem,inode)))*divU + aux2)

                     !$acc loop seq
                     do idime = 1,ndime
                        aux3 = 0.0_rp
                        aux4 = 0.0_rp

                        !$acc loop seq
                        do jdime = 1,ndime
                           aux3 = aux3 + gradU(idime,jdime)*q(connec(ielem,inode),jdime)
                           aux4 = aux4 + u(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime)*gradRho(jdime)
                        end do

                        Re_mom(inode,idime) = Re_mom(inode,idime) + gpvol(1,igaus,ielem)*Ngp(igaus,inode)*0.5_rp*(divF(idime) + q(connec(ielem,inode),idime)*divU + aux3 + aux4 + gradP(idime))
                     end do
                  end do
               end do

               !
               ! Assembly
               !
               do inode = 1,4
                  !$acc atomic update
                  Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode)) + Re_mass(inode)
                  !$acc end atomic
                  !$acc atomic update
                  Rener(connec(ielem,inode)) = Rener(connec(ielem,inode)) + Re_ener(inode)
                  !$acc end atomic
                  do idime = 1,ndime
                     !$acc atomic update
                     Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime) + Re_mom(inode,idime)
                     !$acc end atomic
                  end do
               end do
            end do
            !$acc end parallel loop

            call nvtxEndRange

         end subroutine fem_convec

end module elem_convec
