PROGRAM MAIN
    USE CUBLAS
    USE OPENACC
    IMPLICIT NONE
    Type :: matrix 
	Real*8, Allocatable :: vals(:,:)
    End Type matrix
    Type(matrix), allocatable, Target :: matrices(:)
    Type(matrix), allocatable, Target :: lhs(:), rhs(:)
    Integer(4):: rc, nt, ntheta, nmatrices, iter,i, nvec, ntests
    Integer(4):: ndim2,ndim3
    integer :: numHandles, queue
    Real*8 :: alpha, beta

!    Character*1 :: transa='T', transb='N'
    integer(4) :: transa=CUBLAS_OP_T, transb=CUBLAS_OP_N
    Real*8, Pointer :: mmatrix(:,:), lh(:,:), rh(:,:)
    Real*8, Pointer :: lh2(:,:), rh2(:,:), mmatrix2(:,:)
    type(cublasHandle), allocatable, dimension(:) :: cuHandles
    integer(acc_handle_kind) :: stream
    !Real*8, Allocatable :: square(:,:)

    alpha = 1.0d0
    beta  = 0.0d0
    ntests = 20 
    ntheta = 512
    nvec   = 32
    nt = ntheta
    nmatrices = ntheta !ntheta/2
    numHandles = 3
   
    allocate(cuHandles(numHandles)) 
    ! Allocate on the host first
    Allocate(matrices(nmatrices))
    Allocate(     lhs(nmatrices))
    Allocate(     rhs(nmatrices))
    ! NOW, we copy in
    !$ACC enter data copyin(matrices, lhs, rhs)

    Do i = 1, nmatrices
        if (MOD(i,2) .eq. 0) Then
            ndim2 = ntheta -i+1
        Else
            ndim2 = i
        Endif
        Allocate(matrices(i)%vals(1:ntheta,ndim2))
        Allocate(lhs(i)%vals(1:ntheta,1:nvec))
        Allocate(rhs(i)%vals(ndim2,1:nvec))
        matrices(i)%vals(:,:) = 1.0d0
        lhs(i)%vals(:,:) = 1.0d0
        rhs(i)%vals(:,:) = 1.0d0
    !$ACC enter data copyin(matrices(i)%vals)
    !$ACC enter data copyin(lhs(i)%vals)
    !$ACC enter data copyin(rhs(i)%vals)
    Enddo


    
    Do iter = 1, ntests
        write(6,*)iter
        ! Do something with lhs
        ! This mimics the rest of the progam
        Do i = 1, nmatrices
        lhs(i)%vals = lhs(i)%vals+1d0
        Enddo

        !///////////////////////////////////////////////////////////
        ! This portion of the loop mimics the legendre transform routine...



        If (iter.eq.1) Then
            ! create the cublas handles once
            Do i = 1, numHandles
                rc = cublasCreate(cuHandles(i))
                stream = acc_get_cuda_stream(i)
                rc = cublasSetStream(cuHandles(i), stream) 
            Enddo
        Endif

        Do i = 1, nmatrices,2

            !This tells the compiler that matrices, lhs, and rhs are on the device
            mmatrix => matrices(i)%vals
            lh => lhs(i)%vals
            rh => rhs(i)%vals
            mmatrix2 => matrices(i+1)%vals
            lh2 => lhs(i+1)%vals
            rh2 => rhs(i+1)%vals
            !ndim2 = ntheta-i+1
            queue = mod(i/2,numHandles) + 1
            !$ACC update device(lhs(i)%vals,lhs(i+1)%vals) async(queue)
            !if (MOD(i,2) .eq. 0) Then
            !    ndim2 = ntheta -i+1
            !Else
            !    ndim2 = i
            !Endif
            ndim2 = i
            ndim3 = ntheta-i+2

            !$ACC data present(mmatrix,lh,rh,mmatrix2,lh2,rh2)
            rc = cublasDGEMM_v2(cuHandles(queue),transa,transb,ndim2,nvec,ntheta, alpha, mmatrix, &
                &  ntheta,lh , ntheta, beta,rh,ndim2)
            rc = cublasDGEMM_v2(cuHandles(queue),transa,transb,ndim3,nvec,ntheta, alpha, mmatrix2, &
                &  ntheta,lh2 , ntheta, beta,rh2,ndim3)
            !$ACC END DATA
            !Write(6,*) "RH=", rh(1,1)

            !$ACC update host(rhs(i)%vals,rhs(i+1)%vals) async(queue)
        Enddo
        !$ACC WAIT 
        !//////////////////////
    Enddo
    do i = 1, numHandles
         rc = cublasDestroy(cuHandles(i))
    enddo

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
