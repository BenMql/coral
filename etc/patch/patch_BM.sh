#!/bin/bash
sed -i '45i \ \ integer, parameter, public :: IN_CORE_CHEBYSHEV = 4 !BM_patch' src/fft_fftw3_f03.f90
sed -i '148s/.*/\ \ \ else if \(\(format==PHYSICAL_IN_Z\).or.\(format==IN_CORE_CHEBYSHEV\)\) then !BM patch/' src/fft_fftw3_f03.f90
sed -i '165s/.*/\ \ \ else if \(\(format==PHYSICAL_IN_Z\).or.\(format==IN_CORE_CHEBYSHEV\)\) then !BM patch/' src/fft_fftw3_f03.f90
sed -i '215s/.*/\ \ \ else if \(\(format==PHYSICAL_IN_Z\).or.\(format==IN_CORE_CHEBYSHEV\)\) then !BM patch/' src/fft_fftw3_f03.f90
sed -i '523i \ \ \ \ else if \(format == IN_CORE_CHEBYSHEV\) then           ! BM patch' src/fft_fftw3_f03.f90
sed -i '524i \ \ \ \ \ \ \ ! For R2C/C2R tranforms                           ! BM patch' src/fft_fftw3_f03.f90
sed -i '525i \ \ \ \ \ \ \ call r2c_1m_z_plan\(plan\(0,3\), ph, sp\)             ! BM patch' src/fft_fftw3_f03.f90
sed -i '526i \ \ \ \ \ \ \ call c2c_1m_y_plan\(plan\(0,2\), sp, FFTW_FORWARD \)  ! BM patch' src/fft_fftw3_f03.f90
sed -i '527i \ \ \ \ \ \ \ call c2c_1m_y_plan\(plan\(2,2\), sp, FFTW_BACKWARD\)  ! BM patch' src/fft_fftw3_f03.f90
sed -i '528i \ \ \ \ \ \ \ call c2r_1m_z_plan\(plan\(2,3\), sp, ph\)             ! BM patch' src/fft_fftw3_f03.f90
sed -i '852i \ \ \ \ \ \ \    ! ############################################### ! BM patch' src/fft_fftw3_f03.f90
sed -i '853i \ \ \ \ \ \ \    else if \(format==IN_CORE_CHEBYSHEV\) then          ! BM patch' src/fft_fftw3_f03.f90
sed -i '854i \ \ \ \ \ \ \                                                      ! BM patch' src/fft_fftw3_f03.f90
sed -i '855i \ \ \ \ \ \ \       ! ===== 1D FFTs in Z =====                     ! BM patch' src/fft_fftw3_f03.f90
sed -i '856i \ \ \ \ \ \ \       call r2c_1m_z\(in_r,wk13\)                       ! BM patch' src/fft_fftw3_f03.f90
sed -i '857i \ \ \ \ \ \ \                                                      ! BM patch' src/fft_fftw3_f03.f90
sed -i '858i \ \ \ \ \ \ \       ! ===== Swap Z --> Y; 1D FFTs in Y =====       ! BM patch' src/fft_fftw3_f03.f90
sed -i '859i \ \ \ \ \ \ \       if \(dims\(1\)>1\) then                            ! BM patch' src/fft_fftw3_f03.f90
sed -i '860i \ \ \ \ \ \ \          call transpose_z_to_y\(wk13,wk2_r2c,sp\)      ! BM patch' src/fft_fftw3_f03.f90
sed -i '861i \ \ \ \ \ \ \          call c2c_1m_y\(wk2_r2c,-1,plan\(0,2\)\)         ! BM patch' src/fft_fftw3_f03.f90
sed -i '862i \ \ \ \ \ \ \       else  ! out_c==wk2_r2c if 1D decomposition     ! BM patch' src/fft_fftw3_f03.f90
sed -i '863i \ \ \ \ \ \ \          call transpose_z_to_y\(wk13,out_c,sp\)        ! BM patch' src/fft_fftw3_f03.f90
sed -i '864i \ \ \ \ \ \ \          call c2c_1m_y\(out_c,-1,plan\(0,2\)\)           ! BM patch' src/fft_fftw3_f03.f90
sed -i '865i \ \ \ \ \ \ \       end if                                         ! BM patch' src/fft_fftw3_f03.f90
sed -i '866i \ \ \ \ \ \ \                                                      ! BM patch' src/fft_fftw3_f03.f90
sed -i '867i \ \ \ \ \ \ \       ! ===== Swap Y --> X; 1D FFTs in X =====       ! BM patch' src/fft_fftw3_f03.f90
sed -i '868i \ \ \ \ \ \ \       if \(dims\(1\)>1\) then                            ! BM patch' src/fft_fftw3_f03.f90
sed -i '869i \ \ \ \ \ \ \          call transpose_y_to_x\(wk2_r2c,out_c,sp\)     ! BM patch' src/fft_fftw3_f03.f90
sed -i '870i \ \ \ \ \ \ \       end if                                         ! BM patch' src/fft_fftw3_f03.f90
sed -i '871i \ \ \ \ \ \ \       ! the dct transform takes place in the         ! BM patch' src/fft_fftw3_f03.f90
sed -i '872i \ \ \ \ \ \ \       ! calling routine.                             ! BM patch' src/fft_fftw3_f03.f90
sed -i '873i \ \ \ \ \ \ \    ! ############################################### ! BM patch' src/fft_fftw3_f03.f90
sed -i '945i \ \ \ \ \ \ \    ! ################################################### ! BM patch' src/fft_fftw3_f03.f90
sed -i '946i \ \ \ \ \ \ \    else if \(format==IN_CORE_CHEBYSHEV\) then              ! BM patch' src/fft_fftw3_f03.f90
sed -i '947i \ \ \ \ \ \ \                                                          ! BM patch' src/fft_fftw3_f03.f90
sed -i '948i \ \ \ \ \ \ \       ! ===== 1D FFTs in X =====                         ! BM patch' src/fft_fftw3_f03.f90
sed -i '949i #ifdef OVERWRITE                                             ' src/fft_fftw3_f03.f90
sed -i '950i \ \ \ \ \ \ \       ! the dct transform takes place in the             ! BM patch' src/fft_fftw3_f03.f90
sed -i '951i \ \ \ \ \ \ \       ! calling routine.                                 ! BM patch' src/fft_fftw3_f03.f90
sed -i '952i #else                                                       ' src/fft_fftw3_f03.f90
sed -i '953i \ \ \ \ \ \ \       sz = sp%xsz(1)*sp%xsz(2)*sp%xsz(3)                 ! BM patch' src/fft_fftw3_f03.f90
sed -i '954i \ \ \ \ \ \ \       wk1_p = fftw_alloc_complex(sz)                     ! BM patch' src/fft_fftw3_f03.f90
sed -i '955i \ \ \ \ \ \ \       call c_f_pointer(wk1_p, wk1, &                     ! BM patch' src/fft_fftw3_f03.f90
sed -i '956i \ \ \ \ \ \ \                   [sp%xsz(1),sp%xsz(2),sp%xsz(3)])       ! BM patch' src/fft_fftw3_f03.f90
sed -i '957i \ \ \ \ \ \ \       wk1 = in_c                                         ! BM patch' src/fft_fftw3_f03.f90
sed -i '958i \ \ \ \ \ \ \       ! the dct transform takes place in the             ! BM patch' src/fft_fftw3_f03.f90
sed -i '959i \ \ \ \ \ \ \       ! calling routine.                                 ! BM patch' src/fft_fftw3_f03.f90
sed -i '960i #endif                                                      ' src/fft_fftw3_f03.f90
sed -i '961i \ \ \ \ \ \ \                                                          ! BM patch' src/fft_fftw3_f03.f90
sed -i '962i \ \ \ \ \ \ \       ! ===== Swap X --> Y; 1D FFTs in Y =====           ! BM patch' src/fft_fftw3_f03.f90
sed -i '963i \ \ \ \ \ \ \       if (dims(1)>1) then                                ! BM patch' src/fft_fftw3_f03.f90
sed -i '964i #ifdef OVERWRITE                                            ' src/fft_fftw3_f03.f90
sed -i '965i \ \ \ \ \ \ \          call transpose_x_to_y(in_c,wk2_r2c,sp)          ! BM patch' src/fft_fftw3_f03.f90
sed -i '966i #else                                                       ' src/fft_fftw3_f03.f90
sed -i '967i \ \ \ \ \ \ \          call transpose_x_to_y(wk1,wk2_r2c,sp)           ! BM patch' src/fft_fftw3_f03.f90
sed -i '968i #endif                                                      ' src/fft_fftw3_f03.f90
sed -i '969i \ \ \ \ \ \ \          call c2c_1m_y(wk2_r2c,1,plan(2,2))              ! BM patch' src/fft_fftw3_f03.f90
sed -i '970i \ \ \ \ \ \ \       else  ! in_c==wk2_r2c if 1D decomposition          ! BM patch' src/fft_fftw3_f03.f90
sed -i '971i #ifdef OVERWRITE                                            ' src/fft_fftw3_f03.f90
sed -i '972i \ \ \ \ \ \ \          call c2c_1m_y(in_c,1,plan(2,2))                 ! BM patch' src/fft_fftw3_f03.f90
sed -i '973i #else                                                       ' src/fft_fftw3_f03.f90
sed -i '974i \ \ \ \ \ \ \          call c2c_1m_y(wk1,1,plan(2,2))                  ! BM patch' src/fft_fftw3_f03.f90
sed -i '975i #endif                                                      ' src/fft_fftw3_f03.f90
sed -i '976i \ \ \ \ \ \ \       end if                                             ! BM patch' src/fft_fftw3_f03.f90
sed -i '977i \ \ \ \ \ \ \                                                          ! BM patch' src/fft_fftw3_f03.f90
sed -i '978i \ \ \ \ \ \ \       ! ===== Swap Y --> Z; 1D FFTs in Z =====           ! BM patch' src/fft_fftw3_f03.f90
sed -i '979i \ \ \ \ \ \ \       if (dims(1)>1) then                                ! BM patch' src/fft_fftw3_f03.f90
sed -i '980i \ \ \ \ \ \ \          call transpose_y_to_z(wk2_r2c,wk13,sp)          ! BM patch' src/fft_fftw3_f03.f90
sed -i '981i \ \ \ \ \ \ \       else                                               ! BM patch' src/fft_fftw3_f03.f90
sed -i '982i #ifdef OVERWRITE                                            ' src/fft_fftw3_f03.f90
sed -i '983i \ \ \ \ \ \ \          call transpose_y_to_z(in_c,wk13,sp)             ! BM patch' src/fft_fftw3_f03.f90
sed -i '984i #else                                                       ' src/fft_fftw3_f03.f90
sed -i '985i \ \ \ \ \ \ \          call transpose_y_to_z(wk1,wk13,sp)              ! BM patch' src/fft_fftw3_f03.f90
sed -i '986i #endif                                                      ' src/fft_fftw3_f03.f90
sed -i '987i \ \ \ \ \ \ \       end if                                             ! BM patch' src/fft_fftw3_f03.f90
sed -i '988i \ \ \ \ \ \ \       call c2r_1m_z(wk13,out_r)                          ! BM patch' src/fft_fftw3_f03.f90
sed -i '989i \ \ \ \ \ \ \    ! ################################################### ! BM patch' src/fft_fftw3_f03.f90




