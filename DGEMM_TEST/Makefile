################################################################################
#
# Copyright (c) 2016, NVIDIA Corporation.  All rights reserved.
#
# Please refer to the NVIDIA end user license agreement (EULA) associated
# with this source code for terms and conditions that govern your use of
# this software. Any use, reproduction, disclosure, or distribution of
# this software and related documentation outside the terms of the EULA
# is strictly prohibited.
#
################################################################################

TEST 	 = test_dgemm
TEST2 	 = test_dgemm2
TEST3    = test_dgemm3
TEST4    = test_dgemm4
TEST5    = async
TEST6    = async2
AS3      = async3
AS4      = async4

FC       = pgfortran
EXE	 = exe


FCFLAGS  = -fast -ta=tesla:cc50,cuda8.0 -Mpreprocess
PCFLAGS  = -fast -ta=tesla:cc50,cuda8.0,pinned -Mpreprocess
LIBR     = -Mcuda -Mcudalib=cublas ##-lblas
SFLAGS = -O0 -Mpreprocess -Mbounds
SFLAGS = -fast
LIBS = -lblas

all: test serial

test: main.f90
	$(FC) $(FCFLAGS) -o $(TEST).$(EXE) $< $(LIBR)
test2: main2.f90
	$(FC) $(FCFLAGS) -o $(TEST2).$(EXE) $< $(LIBR)
test3: main3.f90
	$(FC) $(FCFLAGS) -o $(TEST3).$(EXE) $< $(LIBR)
test4: main4.f90
	$(FC) $(PCFLAGS) -o $(TEST4).$(EXE) $< $(LIBR)
serial: main_serial.f90
	$(FC) $(SFLAGS) -o serial.$(EXE) $< $(LIBS)
async: cublasAsync.F90
	$(FC) $(PCFLAGS) -o $(TEST5).$(EXE) $< $(LIBR)
async2: cublasAsync2.F90
	$(FC) $(PCFLAGS) -o $(TEST6).$(EXE) $< $(LIBR)
async3: cublasAsync3.F90
	$(FC) $(PCFLAGS) -o $(AS3).$(EXE) $< $(LIBR)
async4: cublasAsync4.F90
	$(FC) $(PCFLAGS) -o $(AS4).$(EXE) $< $(LIBR)

clean:
	@echo 'Cleaning up...'
	@rm -rf *.exe
