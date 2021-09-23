PROGRAM MAIN

    IMPLICIT NONE
	Type :: matrix 
		Real*8, Allocatable :: vals(:,:)
	End Type matrix

    Type(matrix), allocatable, Target :: matrices(:)
    Type(matrix), allocatable, Target :: lhs(:), rhs(:)
    Integer :: nt, ntheta, nmatrices, iter,i, nvec, ntests, ndim2
    Real*8 :: alpha, beta
    Character*1 :: transa='T', transb='N'

    !Real*8, Allocatable :: square(:,:)



    alpha = 1.0d0
    beta  = 0.0d0

    ntests = 100
    ntheta = 128
    ntheta = 512
    nvec   = 32
    nt = ntheta
    nmatrices = ntheta !ntheta/2
    


    ! Allocate on the host first
    Allocate(matrices(nmatrices))

    Allocate(     lhs(nmatrices))

    Allocate(     rhs(nmatrices))
    ! NOW, we copy in


    Do i = 1, nmatrices
        if (MOD(i,2) .eq. 0) Then
            ndim2 = ntheta -i+1
        Else
            ndim2 = i
        Endif
        !ndim2 = ntheta
        Allocate(matrices(i)%vals(1:ntheta,ndim2))
        Allocate(lhs(i)%vals(1:ntheta,1:nvec))
        Allocate(rhs(i)%vals(ndim2,1:nvec))
        matrices(i)%vals(:,:) = 1.0d0
        lhs(i)%vals(:,:) = 1.0d0
        rhs(i)%vals(:,:) = 1.0d0

    Enddo


    
    Do iter = 1, ntests
        write(6,*)iter
        ! do something with lhs and rhs
        Do i = 1, nmatrices
        lhs(i)%vals = lhs(i)%vals+1d0
        Enddo


        Do i = 1, nmatrices
                !This tells the compiler that matrices, lhs, and rhs are on the device
     

                if (MOD(i,2) .eq. 0) Then
                    ndim2 = ntheta -i+1
                Else
                    ndim2 = i
                Endif
                !ndim2 = ntheta
				CALL DGEMM(transa,transb,ndim2,nvec,ntheta, alpha, matrices(i)%vals, &
                    &  ntheta,lhs(i)%vals , ntheta, beta,rhs(i)%vals,ndim2)

        Enddo

    Enddo


    Do i = 1, nmatrices

        DeAllocate(matrices(i)%vals)

        DeAllocate(rhs(i)%vals)

        DeAllocate(lhs(i)%vals)
    Enddo

    DeAllocate(matrices)
END PROGRAM MAIN
