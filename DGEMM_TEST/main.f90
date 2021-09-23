PROGRAM MAIN
!IFDEF _OPENACC
    USE CUBLAS
!ENDIF
    IMPLICIT NONE
	Type :: matrix 
		Real*8, Allocatable :: vals(:,:)
	End Type matrix

    Type(matrix), allocatable :: matrices(:)
    Type(matrix), allocatable :: lhs(:), rhs(:)
    Integer :: ntheta, nmatrices, iter,i, nvec, ntests, ndim2, rc
    Real*8 :: alpha, beta
    Character*1 :: transa='N', transb='N'
    alpha = 1.0d0
    beta  = 1.0d0

    ntests = 32
    ntheta = 128
    nvec   = 32
    nmatrices = 2 !ntheta/2
    

    ! Allocate on the host first
    Allocate(matrices(nmatrices))

    Allocate(     lhs(nmatrices))

    Allocate(     rhs(nmatrices))
    ! NOW, we copy in
    !$ACC enter data copyin(matrices, lhs, rhs)

    Do i = 1, nmatrices
        ndim2 = ntheta-i+1
        Allocate(matrices(i)%vals(1:ntheta,ndim2))
        Allocate(lhs(i)%vals(1:ntheta,1:nvec))
        Allocate(rhs(i)%vals(ndim2,1:nvec))
        matrices(i)%vals(:,:) = 1.0d0
        lhs(i)%vals(:,:) = 1.0
        rhs(i)%vals(:,:) = 1.0
    !$ACC enter data copyin(matrices(i)%vals)
    !$ACC enter data copyin(lhs(i)%vals)
    !$ACC enter data copyin(rhs(i)%vals)
    Enddo


    
    Do iter = 1, ntests
        write(6,*)iter
        ! do something with lhs and rhs
        Do i = 1, nmatrices
        lhs(i)%vals = lhs(i)%vals+1
        Enddo

        !///////////////////////////////////////
        ! This is legendre transform routine...

        Do i = 1, nmatrices
                !This tells the compiler that matrices, lhs, and rhs are on the device
                !!$ACC data copy(matrices,lhs,rhs,matrices(i)%vals,lhs(i)%vals,rhs(i)%vals), copyin(alpha,beta)
                !!$ACC data present(matrices,lhs,rhs)
                ! this is a contunitation of the last line !$ACC    copyin(alpha,beta)
                !!$ACC update device(lhs(i)%vals)
                !$ACC host_data use_device(matrices(i)%vals,lhs(i)%vals,rhs(i)%vals)
                ndim2 = ntheta-i+1
!				CALL DGEMM(transa,transb,ndim2,nvec,ntheta, alpha, matrices(i)%vals, &
!                    &  ntheta,lhs(i)%vals , ntheta, beta,rhs(i)%vals,ndim2)

				CALL cublasDGEMMcu_dpm('N','N',ndim2,nvec,ntheta, alpha, matrices(i)%vals, &
                    &  ntheta,lhs(i)%vals , ntheta, beta,rhs(i)%vals,ndim2)
                !$ACC end host_data
                !!$ACC update host(rhs(i)%vals)
                !!$ACC end data
        Enddo
        
        !//////////////////////
    Enddo

    ! FOLLOW THIS ORDER WHEN DEALLOCATING
    Do i = 1, nmatrices
        !$ACC exit data delete(matrices(i)%vals)
        DeAllocate(matrices(i)%vals)
        !$ACC exit data delete(RHS(i)%vals)
        DeAllocate(rhs(i)%vals)
        !$ACC exit data delete(lhs(i)%vals)
        DeAllocate(lhs(i)%vals)
    Enddo
    !$ACC exit data delete(matrices,rhs,lhs)
    DeAllocate(matrices)
END PROGRAM MAIN
