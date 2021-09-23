PROGRAM MAIN
#IFDEF _OPENACC
    USE CUBLAS
#ENDIF
    IMPLICIT NONE
	Type :: matrix 
		Real*8, Allocatable :: vals(:,:)
	End Type matrix

    Type(matrix), allocatable, Target :: matrices(:)
    Type(matrix), allocatable, Target :: lhs(:), rhs(:)
    Integer :: nt, ntheta, nmatrices, iter,i, nvec, ntests, ndim2, rc
    Real*8 :: alpha, beta
    Character*1 :: transa='N', transb='N'
    Real*8, Pointer :: mmatrix(:,:), lh(:,:), rh(:,:)
    Real*8, Allocatable :: mmatrix2(:,:), lh2(:,:), rh2(:,:)
    !Real*8, Allocatable :: square(:,:)



    alpha = 1.0d0
    beta  = 0.0d0

    ntests = 32
    ntheta = 128
    nvec   = 32
    nt = ntheta
    nmatrices = 2 !ntheta/2
    
    Allocate(mmatrix2(1:ntheta,1:ntheta))
    Allocate(lh2(1:ntheta,1:nvec))
    Allocate(rh2(1:ntheta,1:nvec))

    ! Allocate on the host first
    Allocate(matrices(nmatrices))

    Allocate(     lhs(nmatrices))

    Allocate(     rhs(nmatrices))
    ! NOW, we copy in
    !!$ACC enter data copyin(matrices, lhs, rhs)

    Do i = 1, nmatrices
        ndim2 = ntheta ! -i+1
        Allocate(matrices(i)%vals(1:ntheta,ndim2))
        Allocate(lhs(i)%vals(1:ntheta,1:nvec))
        Allocate(rhs(i)%vals(ndim2,1:nvec))
        matrices(i)%vals(:,:) = 1.0d0
        lhs(i)%vals(:,:) = 1.0d0
        rhs(i)%vals(:,:) = 1.0d0
    !!$ACC enter data copyin(matrices(i)%vals)
    !!$ACC enter data copyin(lhs(i)%vals)
    !!$ACC enter data copyin(rhs(i)%vals)
    Enddo


    
    Do iter = 1, ntests
        write(6,*)iter
        ! do something with lhs and rhs
        Do i = 1, nmatrices
        lhs(i)%vals = lhs(i)%vals+1d0
        Enddo

        !///////////////////////////////////////
        ! This is legendre transform routine...

        Do i = 1, nmatrices
                !This tells the compiler that matrices, lhs, and rhs are on the device
                mmatrix => matrices(i)%vals
                lh => lhs(i)%vals
                rh => rhs(i)%vals
                lh2 = lhs(i)%vals
                rh2 = rhs(i)%vals
                mmatrix2 = matrices(i)%vals
                ndim2 = ntheta
             !$ACC data copyin(mmatrix,lh), copyout(rh)
                !!$ACC data present(matrices,lhs,rhs)
                ! this is a contunitation of the last line !$ACC    copyin(alpha,beta)
                !!$ACC update device(lhs(i)%vals)
                !!$ACC host_data use_device(matrices(i)%vals,lhs(i)%vals,rhs(i)%vals)
				CALL DGEMM(transa,transb,ndim2,nvec,ntheta, alpha, mmatrix, &
                    &  ntheta,lh , ntheta, beta,rh,ndim2)

!				CALL cublasDGEMMcu_dpm(transa,transb,ndim2,nvec,ntheta, alpha, mmatrix, &
!                    &  ntheta,lh , ntheta, beta,rh,ndim2)
                !!$ACC end host_datal
                !!$ACC update host(rhs(i)%vals)
                !$ACC end data
                Write(6,*) "RH=", rh(1,1)
        Enddo
        
        !//////////////////////
    Enddo

    ! FOLLOW THIS ORDER WHEN DEALLOCATING
    Do i = 1, nmatrices
        !!$ACC exit data delete(matrices(i)%vals)
        DeAllocate(matrices(i)%vals)
        !!$ACC exit data delete(RHS(i)%vals)
        DeAllocate(rhs(i)%vals)
        !!$ACC exit data delete(lhs(i)%vals)
        DeAllocate(lhs(i)%vals)
    Enddo
    !!$ACC exit data delete(matrices,rhs,lhs)
    DeAllocate(matrices)
END PROGRAM MAIN
