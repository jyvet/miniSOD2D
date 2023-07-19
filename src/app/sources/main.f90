program main

	use mod_constants
	use mod_nvtx
	use mod_gmsh_indices
	use quadrature_rules
	use elem_hex
	use mod_mesh
	use jacobian_oper
	use elem_convec
	use elem_diffu
	use mpi

    implicit none

	!
	! Mesh vars not in mod_constants
	!
	integer(4), parameter   :: nelem = 10000
	integer(4), parameter   :: npoin = nelem * nnode
	integer(4), allocatable :: connec(:,:)
	real(rp)  , allocatable :: coord(:,:), He(:,:,:,:), gpvol(:,:,:)

	!
	! Element characteristics
	!
	integer(4), allocatable :: gmshIJK(:,:), invAtoIJK(:,:,:), AtoIJK(:), AtoI(:), AtoJ(:), AtoK(:), gmsh2sodIJK(:,:)
	real(rp)  , allocatable :: Ngp(:,:), Ngp_l(:,:), dNgp(:,:,:), dNgp_l(:,:,:), dlxigp_ip(:,:,:)
	real(rp)  , allocatable :: wgp(:), xgp(:,:)

	!
	! Loop variables
	!
	integer(4)              :: ielem, inode, igaus, ipoin, iorder, jnode, i, j, k

	!
	! Case variables and residuals
	!
	real(rp),   allocatable :: u(:,:), q(:,:), rho(:), pr(:), E(:), Tem(:)
	real(rp),   allocatable :: Rmass(:), Rmom(:,:), Rener(:)
	real(rp),   allocatable :: Dmass(:), Dmom(:,:), Dener(:)

	!
	! Fluid properties
	!
	real(rp)                :: Cp, Pra
	real(rp),   allocatable :: mu_fluid(:), mu_e(:,:), mu_sgs(:,:)

	!
	! Timing
	!
	real(8) :: tstart, tend, tconvec, tdiffu, tavg_convec, tavg_diffu

	!
	! MPI vars
	!
	integer(4) :: ierr, myrank, nprocs

	!
	! Initialize the MPI environment
	!
	call MPI_Init(ierr)

	!
	! Get number of processes and process rank
	!
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
	if (nprocs > 1) then
		write(*,*) 'This program only works with 1 MPI process'
		stop
	end if

	!
	! Print run config to screen
	!
	write(*,*) 'Number of elements = ', nelem
	write(*,*) 'Element order = ', porder
	write(*,*) 'Number of nodes per element = ', nnode
	write(*,*) 'Nodes on mesh = ', npoin
	write(*,*)

	!
	! Generate isopar. element
	!

		! Allocate memory for element characteristics vars.
		allocate(gmshIJK(nnode,3), invAtoIJK(porder+1,porder+1,porder+1))
		allocate(AtoIJK(nnode), AtoI(nnode), AtoJ(nnode), AtoK(nnode), gmsh2sodIJK(porder+1,2))
		allocate(Ngp(ngaus,nnode), Ngp_l(ngaus,nnode), dNgp(ndime,nnode,ngaus), dNgp_l(ndime,nnode,ngaus))
		allocate(wgp(ngaus), xgp(ngaus,ndime), dlxigp_ip(ngaus,ndime,porder+1))

		! Generate element gmshIJK
		call nvtxStartRange("Generate gmshIJK Gmsh")
		call genHighOrderHex(porder, gmshIJK)
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, nnode
			write(*,*) 'gmshIJK(',inode,') = ', gmshIJK(inode,:)
		end do
		write(*,*)
#endif

		! Generate gmsh2sodIJK conversion
		call nvtxStartRange("Generate gmsh2sodIJK")
		gmsh2sodIJK(1,:) = [0,1]
		gmsh2sodIJK(2,:) = [porder,2]
		do inode = 3,porder+1
			gmsh2sodIJK(inode,:) = [inode-2,inode]
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, porder+1
			write(*,*) 'gmsh2sodIJK(',inode,') = ', gmsh2sodIJK(inode,:)
		end do
		write(*,*)
#endif

		! Convert gmshIJK to SOD_IJK
		call nvtxStartRange("Convert gmshIJK to SOD_IJK")
		do inode = 1,nnode
			do i = 1,3
				do iorder = 1,porder+1
					if ( gmshIJK(inode,i) .eq. gmsh2sodIJK(iorder,1) ) then
						gmshIJK(inode,i) = gmsh2sodIJK(iorder,2)
						exit
					end if
				end do
			end do
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, nnode
			write(*,*) 'gmshIJK(',inode,') = ', gmshIJK(inode,:)
		end do
		write(*,*)
#endif

		! Generate AtoIJK and invAtoIJK
		call nvtxStartRange("Generate invAtoIJK")
		jnode = 1
		do k = 1,porder+1
			do i = 1,porder+1
				do j = 1,porder+1
					do inode = 1,nnode
						if ( gmshIJK(inode,1) .eq. i .and. gmshIJK(inode,2) .eq. j .and. gmshIJK(inode,3) .eq. k ) then
							invAtoIJK(i,j,k) = inode
							exit
						end if
					end do
#ifdef __DEBUG__
					write(*,*) i, j, k, invAtoIJK(i,j,k)
#endif
				AtoIJK(jnode) = invAtoIJK(i,j,k)
				jnode = jnode + 1
				end do
			end do
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		write(*,*)
		do inode = 1, nnode
			write(*,*) 'AtoIJK(',inode,') = ', AtoIJK(inode)
		end do
		write(*,*)
#endif

		! Get separate I, J, K arrays per node
		call nvtxStartRange("Get separate I, J, K arrays per node")
		do inode = 1,nnode
			AtoI(inode) = gmshIJK(inode,1)
			AtoJ(inode) = gmshIJK(inode,2)
			AtoK(inode) = gmshIJK(inode,3)
		end do
		call nvtxEndRange

		! Generate isopar. coordinates and GLL quadrature info
		call nvtxStartRange("Generate isopar. coordinates and GLL quadrature info")
		call GaussLobattoLegendre_hex(porder,ngaus,AtoIJK,xgp,wgp)
		call nvtxEndRange
#ifdef __DEBUG__
		write(*,*)
		do igaus = 1,ngaus
			write(*,*) 'xgp(',igaus,') = ', xgp(igaus,:)
		end do
		write(*,*)
		do igaus = 1,ngaus
			write(*,*) 'wgp(',igaus,') = ', wgp(igaus)
		end do
		write(*,*)
#endif

		! Generate shape functions and derivatives
		call nvtxStartRange("Generate shape functions and derivatives")
		do igaus = 1,ngaus
			call hex_highorder(porder,nnode,xgp(igaus,1),xgp(igaus,2),xgp(igaus,3), &
								AtoIJK,Ngp(igaus,:),dNgp(:,:,igaus), &
								Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do igaus = 1,ngaus
			do inode = 1,nnode
				write(*,*) 'Ngp(',igaus,',',inode,') = ', Ngp(igaus,inode)
			end do
		end do
		write(*,*)
#endif

	!
	! Generate mesh
	!
	allocate(connec(nelem,nnode))
	allocate(coord(npoin,ndime))
	call nvtxStartRange("Generate mesh")
	call create_mesh(nelem, npoin, xgp, connec, coord)
	call nvtxEndRange
#ifdef __DEBUG__
	do ielem = 1, nelem
		do inode = 1, nnode
			write(*,*) 'connec(',ielem,',',inode,') = ', connec(ielem,inode)
		end do
	end do
	write(*,*)
	do ipoin = 1, npoin
		write(*,*) 'coord(',ipoin,') = ', coord(ipoin,:)
	end do
	write(*,*)
#endif

	!
	! Compute Jcobian info
	!
	allocate(He(ndime,ndime,ngaus,nelem))
	allocate(gpvol(1,ngaus,nelem))
	call nvtxStartRange("Compute Jcobian info")
	call elem_jacobian(nelem,npoin,connec,coord,dNgp,wgp,gpvol,He)
	call nvtxEndRange
#ifdef __DEBUG__
	do ielem = 1,nelem
		do igaus = 1,ngaus
			do i = 1,ndime
				do j = 1,ndime
					write(*,*) 'He(',i,',',j,',',igaus,',',ielem,') = ', He(i,j,igaus,ielem)
				end do
			end do
		end do
	end do
	write(*,*)
	do ielem = 1,nelem
		do igaus = 1,ngaus
			write(*,*) 'gpvol(',igaus,',',ielem,') = ', gpvol(1,igaus,ielem)
		end do
	end do
	write(*,*)
#endif

	!
	! Generate initial conditions
	!
	allocate(u(npoin,ndime), q(npoin,ndime), rho(npoin), pr(npoin), E(npoin), Tem(npoin))
	allocate(mu_fluid(npoin), mu_e(nelem,nnode), mu_sgs(nelem,nnode))
	call nvtxStartRange("Generate initial conditions")
	!$acc kernels
	u(:,:) = 1.0_rp
	q(:,:) = 1.0_rp
	rho(:) = 1.0_rp
	pr(:)  = 1.0_rp
	E(:)   = 1.0_rp
	Tem(:)   = 1.0_rp
	mu_fluid(:) = 1.0_rp
	mu_e(:,:) = 1.0_rp
	mu_sgs(:,:) = 1.0_rp
	!$acc end kernels
	call nvtxEndRange

	!
	! Fluid properties
	!
	Cp = 1.0_rp
	Pra = 1.0_rp

	!
	! Call the convective term multiple times
	!
	allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
	allocate(Dmass(npoin), Dmom(npoin,ndime), Dener(npoin))
	open(unit=1,file='timers.dat',status='replace')
	tavg_convec = 0.0d0
	tavg_diffu = 0.0d0
	call nvtxStartRange("Loop kernels")
	do i = 1,100
		call nvtxStartRange("Call convective term")
		tstart = MPI_Wtime()
		call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)
		tend = MPI_Wtime()
		tconvec = tend-tstart
		tavg_convec = tavg_convec + tconvec
		call nvtxEndRange
		call nvtxStartRange("Call diffusive term")
		tstart = MPI_Wtime()
		call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,Cp,Pra,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Dmass,Dmom,Dener)
		tend = MPI_Wtime()
		tdiffu = tend-tstart
		tavg_diffu = tavg_diffu + tdiffu
		write(1,*) i, tconvec, tdiffu
		call nvtxEndRange
	end do
	call nvtxEndRange
	close(1)

	!
	! Write avg times to screen
	!
	tavg_convec = tavg_convec / 100.0d0
	tavg_diffu = tavg_diffu / 100.0d0
	write(*,*) 'Avg. convective time = ', tavg_convec
	write(*,*) 'Avg. diffusive time = ', tavg_diffu

	!
	! Write results to file
	!
#ifdef __DEBUG__
	call nvtxStartRange("Write results to file")
	open(unit=10,file='results_convec.dat',status='replace')
	do ipoin = 1,npoin
		write(10,*) ipoin, Rmass(ipoin), Rmom(ipoin,1), Rmom(ipoin,2), Rmom(ipoin,3), Rener(ipoin)
	end do
	close(10)
	open(unit=11,file='results_diffu.dat',status='replace')
	do ipoin = 1,npoin
		write(10,*) ipoin, Dmass(ipoin), Dmom(ipoin,1), Dmom(ipoin,2), Dmom(ipoin,3), Dener(ipoin)
	end do
	close(11)
	call nvtxEndRange
#endif

	!
	! Finalize MPI environment
	!
	call MPI_Finalize(ierr)

end program main