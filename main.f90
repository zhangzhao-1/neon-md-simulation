program final
    implicit none
    ! 参数
    integer, parameter:: N = 100                 ! 原子数
    real*8, parameter:: length = 1.0e-7          ! 边长
    real*8, parameter:: mass = 6.668e-27         ! He原子质量
    real*8, parameter:: kb = 1.380649D-23        ! 玻尔兹曼常数 
    real*8, parameter:: pi = 4.0d0*atan(1.0d0)   
    real*8, parameter:: sigma = 1.48e-10         ! σ
    real*8, parameter:: eps = 1.48e-22            ! ε
    real*8, parameter :: dt = 1.0e-16            ! 时间步长
    integer, parameter :: nsteps = 500000        ! 总模拟步数
    integer, parameter :: n_eq = 20000           ! 平衡步数
    integer, parameter :: diagstep = 50000       ! 输出间隔
    integer, parameter :: r2_diagstep = 10000    ! 粒子间距诊断间隔

    ! 温度数组
    real*8, dimension(5) :: temp_array = (/200.0d0, 250.0d0, 300.0d0, 350.0d0, 400.0d0/)
    real*8 :: TT_current  ! 当前温度
    integer :: temp_idx    

    ! 核心变量
    integer :: i, j, k, seedsize, seedd, step, atomi, idx_max, iter
    real*8 :: xx(N), yy(N), zz(N)                ! 当前位置
    real*8 :: xxpre(N), yypre(N), zzpre(N)       ! 前一时刻位置
    real*8 :: xxnew(N), yynew(N), zznew(N)       ! 下一时刻位置
    real*8 :: vx(N), vy(N), vz(N)                ! 速度
    real*8 :: fx(N), fy(N), fz(N)                ! 力
    real*8 :: accelx(N), accely(N), accelz(N)    ! 加速度

    ! 能量与统计变量
    real*8 :: E_kin, E_pot, E_tot, T_inst
    real*8 :: meanE, M2, delta, delta2, varE, cv_fluc
    integer :: n_sample                          ! 采样步数

    real*8 :: min_r2, max_v, tmpv, sum_vx, sum_vy, sum_vz
    real*8 :: xstep, ystep, zstep, randd, lambda, ll
    real*8 :: smal, dx, dy, dz  
    integer, allocatable :: seed(:)

    ! 输出文件
    integer :: iunit
    character(len=100) :: filename

    ! 获取随机种子
    print*, "Please enter an integer as the random seed :"
    read(*, *) seedd

    ! 打开输出文件
    filename = 'heat_capacity_results.txt'
    open(newunit=iunit, file=filename, status='replace')
    write(iunit, '(A)') 'Temperature(K)    Heat_Capacity(J/K)    Heat_Capacity_per_atom(J/K/atom)'
    close(iunit)

    ! 温度循环
    do temp_idx = 1, size(temp_array)
        TT_current = temp_array(temp_idx)
        
        print "(A,F8.4)", 'temperature: ', TT_current

        ! 重置统计变量
        meanE = 0.0d0    
        M2 = 0.0d0       
        n_sample = 0     

        ! 距离初始化，考虑缩小原子初始间距
        smal = 2.0d0 * 2.0d0*sigma 
        xstep=smal/5.0d0  
        ystep=smal/5.0d0
        zstep=smal/4.0d0
        atomi = 1
        do i = 1, 5
            do j = 1, 5
                do k = 1, 4
                    if (atomi <= N) then
                        call random_number(randd)     ! 位置随机排列，提高代码可重复性
                        xx(atomi)=(length/2.0d0-smal/2.0d0)+(i-1)*xstep+randd*0.1d0*xstep
                        call random_number(randd)
                        yy(atomi)=(length/2.0d0-smal/2.0d0)+(j-1)*ystep+randd*0.1d0*ystep
                        call random_number(randd)
                        zz(atomi) = (length/2.0d0-smal/2.0d0)+(k-1)*zstep+randd*0.1d0*zstep
                        atomi = atomi + 1
                    end if
                end do
            end do
        end do

        ! 初始化速度
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        seed = seedd  ! 所有温度使用相同的种子
        call random_seed(put=seed)
        deallocate(seed)

        do atomi = 1, N
            ll=box()
            vx(atomi)=ll*sqrt(kb*TT_current/mass)
            ll=box()
            vy(atomi)=ll*sqrt(kb*TT_current/mass)
            ll=box()
            vz(atomi)=ll*sqrt(kb*TT_current/mass)
        end do

        ! 速度缩放
        E_kin=0.5d0*mass*sum(vx**2+vy**2+vz**2) ! 理论动能
        lambda=sqrt((1.5d0*N*TT_current*kb)/E_kin) ! 放缩因子
        vx=vx*lambda
        vy=vy*lambda
        vz=vz*lambda

        ! 初始动能
        E_kin = 0.5d0 * mass * sum(vx**2 + vy**2 + vz**2)
        print'(A, ES15.8)', "Initial E_kin:", E_kin
        print'(A, ES15.8)', "Theoretical E_kin:", 1.5d0*N*TT_current*kb
        if (abs(E_kin-1.5d0*N*TT_current*kb)/abs(1.5d0*N*TT_current*kb) > 0.01d0) then
            print *, "WARNING: Initial kinetic energy mismatch! Difference = ", abs(E_kin-1.5d0*N*TT_current*kb)
        end if

        ! 初始力与加速度
        call force(xx, yy, zz, fx, fy, fz, E_pot, N, sigma, eps, length, min_r2, 0)
        accelx=fx/mass
        accely=fy/mass
        accelz=fz/mass

        do atomi = 1, N
            xxpre(atomi)=xx(atomi)-vx(atomi)*dt+ 0.5d0*accelx(atomi)*dt**2
            yypre(atomi)=yy(atomi)-vy(atomi)*dt+ 0.5d0*accely(atomi)*dt**2
            zzpre(atomi)=zz(atomi)-vz(atomi)*dt+ 0.5d0*accelz(atomi)*dt**2
            xxpre(atomi)=modulo(xxpre(atomi)+length, length)
            yypre(atomi)=modulo(yypre(atomi)+length, length)
            zzpre(atomi)=modulo(zzpre(atomi)+length, length)
        end do

        ! 上面都是初始化，下面进入正式计算
        print *, " "
        print  "(A,I6,A,I6)", "Start simulation: Total steps ", nsteps, ", Equilibrium steps", n_eq
        print *, " "
        do step = 1, nsteps
            ! 更新位置
            do atomi = 1, N
                xxnew(atomi) = 2*xx(atomi) - xxpre(atomi) + accelx(atomi)*dt**2
                yynew(atomi) = 2*yy(atomi) - yypre(atomi) + accely(atomi)*dt**2
                zznew(atomi) = 2*zz(atomi) - zzpre(atomi) + accelz(atomi)*dt**2
                xxnew(atomi) = modulo(xxnew(atomi) + length, length)
                yynew(atomi) = modulo(yynew(atomi) + length, length)
                zznew(atomi) = modulo(zznew(atomi) + length, length)
            end do

            ! 计算力，势能
            call force(xxnew, yynew, zznew, fx, fy, fz, E_pot, N, sigma, eps, length, min_r2, step)

            ! 更新速度，加速度
            do atomi = 1, N
                dx=xxnew(atomi)-xxpre(atomi)
                dx=dx-length*nint(dx/length)
                vx(atomi)=dx/(2.0d0*dt)
                dy=yynew(atomi)-yypre(atomi)
                dy= dy-length*nint(dy/length)
                vy(atomi)=dy/(2.0d0*dt)
                dz=zznew(atomi)-zzpre(atomi)
                dz=dz-length*nint(dz/length)
                vz(atomi)=dz/(2.0d0*dt)
                accelx(atomi)=fx(atomi)/mass
                accely(atomi)=fy(atomi)/mass
                accelz(atomi)=fz(atomi)/mass
            end do

            ! 更新位置
            xxpre = xx
            yypre = yy
            zzpre = zz
            xx = xxnew
            yy = yynew
            zz = zznew

            ! 计算能量
            E_kin=0.5d0*mass*sum(vx**2+vy**2+vz**2) ! 动能
            E_tot=E_kin + E_pot ! 总能量
            T_inst=(2.0d0*E_kin)/(3.0d0*N*kb)

            ! 整个过程中要保持温度在目标温度
            if (step <= n_eq) then ! 平衡阶段
                do iter = 1, 5
                    if (abs(T_inst-TT_current)>5.0d0) then
                        lambda=sqrt(TT_current/T_inst)
                        vx=vx*lambda
                        vy=vy*lambda
                        vz=vz*lambda
                        ! 重新计算动能和温度
                        E_kin = 0.5d0*mass*sum(vx**2+vy**2+vz**2)
                        T_inst=(2.0d0*E_kin)/(3.0d0*N*kb)
                    else 
                        exit
                    end if
                end do
            else
                if (abs(T_inst-TT_current)>10.0d0) then
                    lambda=sqrt(TT_current/T_inst)
                    vx=vx*lambda
                    vy=vy*lambda
                    vz=vz*lambda
                    ! 重新计算动能和温度
                    E_kin=0.5d0*mass*sum(vx**2+vy**2+vz**2)
                    T_inst=(2.0d0*E_kin)/(3.0d0*N*kb)
                end if
            end if

            ! 消除漂移速度
            sum_vx = sum(vx)
            sum_vy = sum(vy)
            sum_vz = sum(vz)
            vx = vx - sum_vx / real(N, 8)
            vy = vy - sum_vy / real(N, 8)
            vz = vz - sum_vz / real(N, 8)

            ! 最终动能
            E_kin=0.5d0*mass*sum(vx**2+vy**2+vz**2)
            T_inst=(2.0d0*E_kin)/(3.0d0*N*kb)
            E_tot=E_kin + E_pot

            ! Welford算法
            if (step > n_eq) then
                n_sample=n_sample+1
                delta=E_tot-meanE
                meanE=meanE+delta/dble(n_sample)
                delta2=E_tot-meanE
                M2=M2+delta*delta2
            end if

            ! 设置最大速度，防止数值爆炸
            max_v = 0.0d0
            idx_max = -1
            do i = 1, N
                tmpv = sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
                if (tmpv > max_v) then
                    max_v = tmpv
                    idx_max = i
                end if
            end do
            if (max_v > 1.0d5) then
                print *, "ERROR: Too large max_v at step", step, "max_v = ", max_v
                stop
            end if

            ! 300K的时候每diagstep输出一次中间结果
            if (TT_current == 300.0d0 .and. mod(step, diagstep) == 0) then
                write(*, '(I8, 3X, A, F8.2, 3X, A, ES12.5, 3X, A, ES12.5, 3X, A, ES12.5)') &
                    step, "T_inst", T_inst, "E_kin", E_kin, "E_pot", E_pot, "E_tot", E_tot
            end if
        end do  

        ! 计算热容 - 修正公式
        if (n_sample > 1) then
            varE = M2/dble(n_sample-1)  ! 使用n-1进行无偏估计
            cv_fluc = varE/(kb*TT_current**2)
            
            ! 检查热容量级是否合理
            if (cv_fluc < 0) then
                print *, "WARNING: Negative heat capacity detected! Using absolute value."
                cv_fluc = abs(cv_fluc)
            end if

            ! 输出最终结果 
            print *, " "
            print *, "Simulation Results T", TT_current, "K"
            print '(A, I8)', "Total steps", n_sample
            print '(A, ES15.8)', "<E>", meanE
            print '(A, ES17.8)', "<E*E>-<E>*<E>", varE
            print '(A, ES16.8)', "Cv", cv_fluc
            print '(A, ES14.8)', "Cv/N ", cv_fluc / N
            print '(A, ES16.8)', "Theoretical Cv", 1.5d0*N*kb
            ! print '(A, ES15.8)', "error", abs(1.5d0*N*kb-cv_fluc)
            print *, " "

            ! 检查热容量级,如果和标准值差十倍以上会警告,实验证明警告挺常见的
            if (abs(cv_fluc - 1.5d0*N*kb) / (1.5d0*N*kb) > 10.0d0) then
                print *, "Warning:Expected", 1.5d0*N*kb, "but got", cv_fluc
            end if

            ! 写入文件
            open(newunit=iunit, file=filename, status='old', position='append')
            write(iunit, '(F8.2, 5X, ES16.8, 5X, ES16.8)') TT_current, cv_fluc, cv_fluc/N
            close(iunit)
        else
            print *, "ERROR: Not enough valid sampling for T = ", TT_current
        end if

    end do  ! 结束温度循环

    contains
    ! Box-Muller变换
    function box() result(ran)
        implicit none
        real*8 :: ran, u1, u2
        call random_number(u1)
        call random_number(u2)
        if (u1 < 1.0d-10) u1 = 1.0d-10
        ran = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * pi * u2)
    end function box

    subroutine force(xxx, yyy, zzz, ffx, ffy, ffz, E_pot_out, NN, sigma1, eps1, length1, min_r2_out, sstep)
        implicit none
        integer, intent(in) :: NN, sstep
        real*8, intent(in)  :: xxx(NN), yyy(NN), zzz(NN),sigma1, eps1, length1
        real*8, intent(out) :: ffx(NN), ffy(NN), ffz(NN),E_pot_out, min_r2_out
        integer :: ii, jj
        real*8 :: ddx, ddy, ddz, r2, r_inv2, r_inv6, r_inv12,V_lj, f_lj_over_r, f_pair,r_cut2, r_safe2, V_cut 

        ! 截断以及参数初始化
        r_cut2=(3.0d0*sigma1)**2
        r_safe2=(0.8d0*sigma1)**2
        V_cut=4.0d0*eps1*((sigma1**12/r_cut2**6)-(sigma1**6/r_cut2**3))
        ffx = 0.0d0
        ffy = 0.0d0
        ffz = 0.0d0
        E_pot_out = 0.0d0
        min_r2_out = 1.0d99

        ! 计算相互作用
        do ii = 1, NN-1
            do jj = ii+1, NN
                ! 镜像法计算粒子间距
                ddx=xxx(jj)-xxx(ii)
                ddy=yyy(jj)-yyy(ii)
                ddz=zzz(jj)-zzz(ii)
                ddx=ddx-length1*nint(ddx/length1)
                ddy=ddy-length1*nint(ddy/length1)
                ddz=ddz-length1*nint(ddz/length1)

                r2=ddx*ddx+ddy*ddy+ddz*ddz
                min_r2_out = min(min_r2_out, r2)

                ! 仅处理截断距离内的相互作用
                if (r2 > r_cut2) cycle
                if (r2 < r_safe2) r2 = r_safe2  ! 避免数值发散

                ! 计算势能
                r_inv2=1.0d0/r2
                r_inv6=r_inv2**3
                r_inv12=r_inv6**2
                V_lj=4.0d0*eps1*(sigma1**12*r_inv12-sigma1**6*r_inv6)-V_cut
                E_pot_out=E_pot_out+V_lj 

                ! 计算力
                f_lj_over_r= 48.0d0*eps1*(sigma1**12*r_inv12-0.5d0*sigma1**6*r_inv6)*r_inv2
                f_pair = f_lj_over_r

                ! 分配力到两个粒子
                ffx(ii)=ffx(ii)+f_pair*ddx
                ffy(ii)=ffy(ii)+f_pair*ddy
                ffz(ii)=ffz(ii)+f_pair*ddz
                ffx(jj)=ffx(jj)-f_pair*ddx
                ffy(jj)=ffy(jj)-f_pair*ddy
                ffz(jj)=ffz(jj)-f_pair*ddz
            end do
        end do

        ! 用来输出最小粒子间距,如果报错可以打开这一部分找出问题在哪一步
        ! if (mod(sstep, r2_diagstep)==0) then
            ! print*, "Step ", sstep, " Min r2: ", min_r2_out, "  3sigma*sigma : ", r_cut2
        ! end if

        ! 检查势能是否稳定
        if (E_pot_out /= E_pot_out) then
            print *, "ERROR: NaN potential energy at step", sstep
            stop
        end if
    end subroutine force

end program final

