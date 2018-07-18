program FF_HI_semimp
use omp_lib
implicit none
real*8, parameter :: PI = 4*atan(1.0D0)
integer*8 :: Nblocks, Nsteps, steps, block, time(1:8), seed, i, Ntraj, VR_opt
integer :: NTimeSteps, timestep, trajectories, Ndtwidths
real*8 :: sr, b, h, a, Q0, Nrelax_times, dt, Ql, Ql2, F(3), dW(3), Qtemp(3), Bq, Bs, delX(3)
real*8 :: a2, a4, time1, time2, B_eta, Bpsi, Bpsi2
real*8, dimension(3,3) :: k, delT, tau
!large arrays must be declared allocatable so they go into heap, otherwise
!OpenMP threads run out of memory
real*8, dimension(:, :), allocatable :: Q, Q_eq_VR
real*8, dimension(:), allocatable :: timestepwidths, Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2
real*8, dimension(:), allocatable :: Qavg, Vqavg, S, Serr

call date_and_time(values=time)
seed = time(8)*100 + time(7)*10

open(unit=31, file='inputparameters.inp')
open(unit=30, file='timestepdata.inp')
open(unit=32, file='options.inp')

read (32, *) VR_opt
read (31, *) sr, b, h, Q0
read (30, *) Ntraj, Ndtwidths, Nrelax_times
allocate(timestepwidths(Ndtwidths))
do i=1,Ndtwidths
    read(30, *) timestepwidths(i)
end do

allocate(Q(3,1:Ntraj))
if (VR_opt.eq.1) then
    allocate(Q_eq_VR(3,1:Ntraj))
end if
allocate(Aeta(Ndtwidths), Veta(Ndtwidths), Apsi(Ndtwidths), Vpsi(Ndtwidths), Apsi2(Ndtwidths), Vpsi2(Ndtwidths))
allocate(Qavg(Ndtwidths), Vqavg(Ndtwidths), S(Ndtwidths), Serr(Ndtwidths))

Aeta = 0.D0
Apsi = 0.D0
Apsi2 = 0.D0
Veta = 0.D0
Vpsi = 0.D0
Vpsi2 = 0.D0
Qavg = 0.D0
Vqavg = 0.D0
S = 0.D0
Serr = 0.D0

a = h*sqrt(PI)
a2 = a**2
a4 = a2**2

k(:,:) = 0.D0
k(1,2) = 1.D0
k = sr*k

delT(:,:) = 0.D0
delT(1,1) = 1.D0
delT(2,2) = 1.D0
delT(3,3) = 1.D0

delX = (/1.D0, 0.D0, 0.D0/)

!open(unit=21, file='Qresults_Q0=10.dat')

looptimesteps: do timestep=1,Ndtwidths

    dt = timestepwidths(timestep)
    NtimeSteps = int(Nrelax_times/dt)

    !DOES NOT WORK WHEN Q0-sqrt(b) < 0
    Q(:,:) = generate_Q_FF(Q0, b, Ntraj, seed, 10000)

    !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, VR_opt)
    !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

    !equilibration for 10 relaxation times or 100 time steps, whichever is larger
    do steps=1,max(int(10/dt),100)
        k(1,2) = 0

        !$OMP DO schedule(dynamic, 100)
        do i=1,Ntraj
            dW = Wiener_step(seed, dt)
            Q(:,i) =  step(Q(:,i), k, dt, Q0, delT, b, a, dW)
        end do
        !$OMP END DO

    end do

    if (VR_opt.eq.0) then
        do steps=1,Ntimesteps
            k(1,2) = sr
            !$OMP DO schedule(dynamic, 100)
            do i=1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, Q0, delT, b, a, dW)
            end do
            !$OMP END DO
        end do

    elseif (VR_opt.eq.1) then
        !call step_with_VR()
        !$OMP single
        Q_eq_VR(:,:) = Q(:,:)
        !$OMP end single

        do steps=1,Ntimesteps

            !$OMP DO schedule(dynamic, 100)
            do i=1,Ntraj
                dW = Wiener_step(seed, dt)
                k(1,2) = sr
                Q(:,i) =  step(Q(:,i), k, dt, Q0, delT, b, a, dW)
                k(1,2) = 0.D0
                Q_eq_VR(:,i) = step(Q_eq_VR(:,i), k, dt, Q0, delT, b, a, dW)
            end do
            !$OMP END DO
        end do

    else
        print *, "Variance Reduction option not set, place in options.inp file"
    end if

    !$OMP END PARALLEL

    if (VR_opt.eq.0) then
        call measure_no_VR()
    elseif (VR_opt.eq.1) then
        call measure_with_VR()
    else
        print *, "Variance Reduction option not set, place in options.inp file"
    end if

end do looptimesteps

!close(unit=21)

deallocate(Q)
if (VR_opt.eq.1) then
    deallocate(Q_eq_VR)
end if

Aeta = Aeta/(Ntraj*sr)
Veta = Veta/(Ntraj*sr**2)
Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

Apsi = Apsi/(Ntraj*sr**2)
Vpsi = Vpsi/(Ntraj*sr**4)
Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

Apsi2 = Apsi2/(Ntraj*sr**2)
Vpsi2 = Vpsi2/(Ntraj*sr**4)
Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

Qavg = sqrt(Qavg/Ntraj)
Vqavg = Vqavg/Ntraj
Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

S = S/Ntraj
Serr = Serr/Ntraj
Serr = sqrt((Serr - S**2)/(Ntraj-1))

call Write_data()

contains

subroutine measure_with_VR()
    implicit none
        !Measurements
    tau = 0

    do i=1,Ntraj
        !Add from shear-flow dumbbell
        Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
        Ql = sqrt(Ql2)
        Bs = dot_product(delX, Q(:,i))**2/Ql2
        F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:,i)/Ql
        tau(:,:) = dyadic_prod(Q(:,i), F)

        Qavg(timestep) = Qavg(timestep) + Ql2
        Vqavg(timestep) = Vqavg(timestep) + Ql
        S(timestep) = S(timestep) + 0.5*(3*Bs - 1)
        Serr(timestep) = Serr(timestep) + 0.25*(9*Bs**2 - 6*Bs + 1)

        B_eta = tau(1,2)
        Bpsi = (tau(1,1) - tau(2,2))
        Bpsi2 = (tau(2,2) - tau(3,3))
!        Aeta(timestep) = Aeta(timestep) + B_eta
!        Apsi(timestep) = Apsi(timestep) + Bpsi
!        Apsi2(timestep) = Apsi2(timestep) + Bpsi2
!        Veta(timestep) = Veta(timestep) + B_eta**2
!        Vpsi(timestep) = Vpsi(timestep) + Bpsi**2
!        Vpsi2(timestep) = Vpsi2(timestep) + Bpsi2**2

        !subtract from equilibrium Dumbbell
        Ql2 = Q_eq_VR(1,i)**2 + Q_eq_VR(2,i)**2 + Q_eq_VR(3,i)**2
        Ql = sqrt(Ql2)
        Bs = dot_product(delX, Q_eq_VR(:,i))**2/Ql2
        F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q_eq_VR(:,i)/Ql
        tau(:,:) = dyadic_prod(Q_eq_VR(:,i), F)

        !Qavg(timestep) = Qavg(timestep) - Ql2
        !Vqavg(timestep) = Vqavg(timestep) - Ql
        S(timestep) = S(timestep) - 0.5*(3*Bs - 1)
        Serr(timestep) = Serr(timestep) - 0.25*(9*Bs**2 - 6*Bs + 1)

        B_eta = B_eta - tau(1,2)
        Bpsi = Bpsi - (tau(1,1) - tau(2,2))
        Bpsi2 = Bpsi2 - (tau(2,2) - tau(3,3))
        Aeta(timestep) = Aeta(timestep) + B_eta
        Apsi(timestep) = Apsi(timestep) + Bpsi
        Apsi2(timestep) = Apsi2(timestep) + Bpsi2
        Veta(timestep) = Veta(timestep) + B_eta**2
        Vpsi(timestep) = Vpsi(timestep) + Bpsi**2
        Vpsi2(timestep) = Vpsi2(timestep) + Bpsi2**2
    end do

end subroutine

subroutine measure_no_VR()
    implicit none
        !Measurements
    tau = 0

    do i=1,Ntraj
        Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
        Ql = sqrt(Ql2)
        Bs = dot_product(delX, Q(:,i))**2/Ql2
        F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:,i)/Ql
        tau(:,:) = dyadic_prod(Q(:,i), F)

        Qavg(timestep) = Qavg(timestep) + Ql2
        Vqavg(timestep) = Vqavg(timestep) + Ql
        S(timestep) = S(timestep) + 0.5*(3*Bs - 1)
        Serr(timestep) = Serr(timestep) + 0.25*(9*Bs**2 - 6*Bs + 1)

        B_eta = tau(1,2)
        Bpsi = (tau(1,1) - tau(2,2))
        Bpsi2 = (tau(2,2) - tau(3,3))
        Aeta(timestep) = Aeta(timestep) + B_eta
        Apsi(timestep) = Apsi(timestep) + Bpsi
        Apsi2(timestep) = Apsi2(timestep) + Bpsi2
        Veta(timestep) = Veta(timestep) + B_eta**2
        Vpsi(timestep) = Vpsi(timestep) + Bpsi**2
        Vpsi2(timestep) = Vpsi2(timestep) + Bpsi2**2
    end do

end subroutine

subroutine Write_data()
    implicit none
    open(unit=20, file='Ql.dat')
    open(unit=21, file='S.dat')
    open(unit=22, file='eta.dat')
    open(unit=23, file='psi.dat')
    open(unit=24, file='psi2.dat')

    !10 format(F4.2,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5)
    10 format(F6.3,4X, 2(E12.5, 2X, E12.5, 4X))
    write(*,*) 'dtw       Q             err             S             err'
    do i = 1,Ndtwidths
        write(*,10) timestepwidths(i), Qavg(i), Vqavg(i), S(i), Serr(i)
    end do

    12 format(F6.3,4X, 3(E12.5, 2X, E12.5, 4X))
    write(*,*) 'dtw       eta           err             psi           err       psi2        err'
    do i = 1,Ndtwidths
        write(*,12) timestepwidths(i), Aeta(i), Veta(i), Apsi(i), Vpsi(i), Apsi2(i), Vpsi2(i)
    end do


    11 format(F6.3,4X, E15.8, 2X, E15.8, 4X)
    do i=20,24
        write(i,*) 'timestepwidth    avg    err'
    end do
    do i = 1,Ndtwidths
        write(20,11) timestepwidths(i), Qavg(i), Vqavg(i)
        write(21,11) timestepwidths(i), S(i), Serr(i)
        write(22,11) timestepwidths(i), Aeta(i), Veta(i)
        write(23,11) timestepwidths(i), Apsi(i), Vpsi(i)
        write(24,11) timestepwidths(i), Apsi2(i), Vpsi2(i)
    end do

    close(unit=20)
    close(unit=21)
    close(unit=22)
    close(unit=23)
    close(unit=24)

end subroutine

function construct_B(Q, a, delT)
    implicit none
    real*8, intent(in) :: Q(3), a
    real*8, intent(in) :: delT(3,3)
    real*8, dimension(3,3) :: construct_B
    real*8 :: Ql, aux, g, g_til, Ql2, Ql4, Ql6
    real*8, save :: C43 = 4.D0/3.D0
    real*8, save :: C83 = 8.D0/3.D0
    real*8, save :: C143 = 14.D0/3.D0

    Ql2 = Q(1)**2 + Q(2)**2 + Q(3)**2
    Ql = sqrt(Ql2)
    Ql4 = Ql2**2
    Ql6 = Ql2**3

    aux = a/(C43*Ql*(Ql2 + C43*a2)**3)
    g = 1.D0-aux*(Ql6 + C143*a2*Ql4 + 8.D0*a4*Ql2)
    g_til = -aux*(Ql6 + 2.D0*a2*Ql4 - C83*a4*Ql2)
    construct_B = sqrt(g)*delT + (sqrt(g+g_til)-sqrt(g))*dyadic_prod(Q,Q)/Ql2
    return

end function construct_B

function shift_xor(val,shift)
    integer*8 :: shift_xor
    integer*8, intent(in) :: val, shift
    shift_xor = ieor(val,ishft(val,shift))
end function

function rand_floats(seed, N)
    implicit none
    integer*8, intent(in) :: N
    integer*8, intent(inout) :: seed
    real*8, dimension(N) :: rand_floats
    integer :: i

    do i=1,N
        !Generates a random number between 0 and 1
        !Using xorshift and one round of 64-bit MCG
        seed = shift_xor(shift_xor(shift_xor(seed, 13_8),-17_8),43_8)
        rand_floats(i) = seed * 2685821657736338717_8 * 5.42101086242752217D-20 + 0.5D0
    end do

end function rand_floats

function Wiener_step(seed, dt)
    implicit none
    integer*8, intent(inout) :: seed
    real*8, intent(in) :: dt
    real*8, dimension(3) :: Wiener_step
    real*8, dimension(3) :: dW
    integer :: i

    do i=1,3
        !Generates a random number between -0.5 and 0.5
        !Using xorshift and one round of 64-bit MCG
        seed = shift_xor(shift_xor(shift_xor(seed, 13_8),-17_8),43_8)
        dW(i) = seed * 2685821657736338717_8 * 5.42101086242752217D-20
    end do

    !Generates an approximately gaussian-distributed number dW
    Wiener_step = dW*sqrt(dt)*(14.14855378D0*dW*dW + 1.21569221D0)
    return

end function Wiener_step

function step(Q, k, dt, Q0, delT, b, a, dW)
    implicit none
    real*8, intent(in) :: Q(3), dt, Q0, b, a, dW(3)
    real*8, intent(in) :: delT(3,3), k(3,3)
    real*8 :: Ql, L, RdotB2_L, Qlength, temp_R1, temp_R2
    real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
    real*8, dimension(3,3) :: B1, B2, B1B1, B2B2
    real*8, dimension(3) :: step

    Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
    F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q/Ql

    B1 = construct_B(Q, a, delT)

    B1B1 = matmul(B1,B1)
    Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
              + ten_vec_dot(B1,dW)

    B2 = construct_B(Qpred, a, delT)
    B2B2 = matmul(B2,B2)
    RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
            + ten_vec_dot(B1,dW)

    L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
    u = RHS/L
    RdotB2 = ten_vec_dot(b*B2B2*dt/4.D0,u)
    RdotB2_L = sqrt(RdotB2(1)**2 + RdotB2(2)**2 + RdotB2(3)**2)

    Qlength = find_roots( -(2.D0*Q0+L), -b+Q0**2-RdotB2_L+2.d0*L*Q0, &
                            RdotB2_L*Q0+L*(b-Q0**2), Q0, sqrt(b) )

    step = RHS*Qlength/L

end function step

pure function find_roots(a, b, c, Q0, lim)
    implicit none
    real*8, intent(in) :: a, b, c, Q0, lim
    real*8 :: find_roots
    real*8 :: Q, R, theta, x
    integer :: i

    Q = (a**2 - 3.D0*b)/9.D0
    R = (2.D0*a**3 - 9.D0*a*b + 27.D0*c)/54.D0

    theta = acos(R/sqrt(Q**3))

    do i=-1,1
        x = -2.D0*sqrt(Q)*cos((theta + real(i)*PI*2.D0)/3.D0)-a/3.D0
        if ((x.ge.(Q0-lim)).and.(x.le.(Q0+lim))) then
            find_roots = x
            EXIT
        end if
   end do

end function find_roots

pure function dyadic_prod(vec1,vec2)
    implicit none
    real*8, dimension(3), intent(in) :: vec1, vec2
    real*8, dimension(3,3) :: dyadic_prod
    integer :: i, j

    do i=1,3
        do j=1,3
            dyadic_prod(i,j) = vec1(i)*vec2(j)
        end do
    end do
end function dyadic_prod

pure function ten_vec_dot(tensor, vector)
    implicit none
    real*8, dimension(3,3), intent(in) :: tensor
    real*8, dimension(3), intent(in) :: vector
    real*8, dimension(3) :: ten_vec_dot
    integer :: i

    do i=1,3
        ten_vec_dot(i) = dot_product(tensor(i,:), vector)
    end do
end function ten_vec_dot

pure function beta(x,y)
    implicit none
    real*8, intent(in) :: x, y
    real*8 :: beta

    beta = log_gamma(x) + log_gamma(y) - log_gamma(x+y)
    beta = exp(beta)
end function beta

function psiQ_FF(Q, b, Q0)
    implicit none
    real*8, intent(in) :: Q, b, Q0
    real*8 :: Jeq
    real*8 :: psiQ_FF

    Jeq = (1.D0/(b+3.D0)+Q0**2/b)*beta(0.5D0,(b+2.D0)/2.D0)*b**(1.5D0)

    psiQ_FF = Q**2*(1.D0-(Q-Q0)**2.D0/b)**(b/2.D0)/Jeq
end function psiQ_FF

function integral_psiQ_FF(Q, b, Q0)
    !Simple cumulative trapezoidal integral of psiQ_FF at
    !points specified in Q
    implicit none
    real*8, dimension(:), intent(in) :: Q
    real*8, intent(in) :: b, Q0
    real*8, dimension(size(Q)) :: integral_psiQ_FF
    integer :: k

    integral_psiQ_FF(1) = 0.D0
    do k=2,size(Q)
        integral_psiQ_FF(k) = integral_psiQ_FF(k-1) + &
                              (psiQ_FF(Q(k-1),b,Q0) + psiQ_FF(Q(k),b,Q0))*(Q(k)-Q(k-1))/2.D0
    end do
    !Trapezoidal rule is far from perfect, but we must have int from 0 to 1
    integral_psiQ_FF = integral_psiQ_FF/integral_psiQ_FF(size(Q))

end function integral_psiQ_FF

function generate_Ql_eq_FF(N, b, Q0, seed, Nsteps)
    implicit none
    real*8, intent(in) :: b, Q0
    integer*8, intent(in) :: N, Nsteps
    integer*8, intent(inout) :: seed
    integer :: k
    real*8, dimension(N) :: generate_Ql_eq_FF, rands
    real*8, dimension(Nsteps) :: Q, intpsiQ
    real*8 :: width

    !Generate a Q vector between Q0-sqrt(b) and Q0+sqrt(b) with Nsteps steps
    !Round off a little at the end to prevent singularities
    width = 2.D0*sqrt(b)/(Nsteps-1)
    Q = 0.D0
    Q(1) = Q0-sqrt(b)
    do k=2,Nsteps
        Q(k) = Q(k-1) + width
    end do
    Q(1) = Q(1) + 0.0000001D0
    Q(Nsteps) = Q(Nsteps) - 0.0000001D0

    intpsiQ = integral_psiQ_FF(Q,b,Q0)

    rands(:) = rand_floats(seed, N)

    !Inverse of integral returns distribution of psiQ_FF
    do k=1,N
        generate_Ql_eq_FF(k) = inverse_lin_interp(Q,intpsiQ,rands(k))
    end do

end function generate_Ql_eq_FF

function inverse_lin_interp(x, fx, fxval)
    implicit none
    real*8, dimension(:), intent(in) :: x, fx
    real*8, intent(in) :: fxval
    real*8 :: inverse_lin_interp
    integer :: i, j

    do i=1,size(x)
        if (fx(i) > fxval) then
            inverse_lin_interp = (fxval-fx(i-1))/(fx(i)-fx(i-1))*(x(i)-x(i-1)) + x(i-1)
            EXIT
        end if
   end do

end function inverse_lin_interp

function generate_Q_FF(Q0,b, N, seed, Nsteps)
    implicit none
    real*8, intent(in) :: Q0, b
    integer*8, intent(inout) :: seed
    integer*8, intent(in) :: N, Nsteps
    real*8, dimension(3,N) :: generate_Q_FF
    real*8 :: Ql(N)

    generate_Q_FF(:,:) = spherical_unit_vectors(N, seed)

    Ql(:) = generate_Ql_eq_FF(N, b, Q0, seed, Nsteps)

    generate_Q_FF(1,:) = generate_Q_FF(1,:)*Ql
    generate_Q_FF(2,:) = generate_Q_FF(2,:)*Ql
    generate_Q_FF(3,:) = generate_Q_FF(3,:)*Ql

end function generate_Q_FF

function spherical_unit_vectors(N, seed)
    implicit none
    integer*8, intent(inout) :: seed
    integer*8, intent(in) :: N
    real*8 :: x1, x2, r1(1), r2(1)
    real*8, dimension(3,N) :: spherical_unit_vectors

    do i=1,N
        r1 = rand_floats(seed,1)
        r2 = rand_floats(seed,1)
        x1 = r1(1)*2.D0 - 1.D0
        x2 = r2(1)*2.D0 - 1.D0
        do while ((x1**2 + x2**2).ge.1.D0)
            r1 = rand_floats(seed,1)
            r2 = rand_floats(seed,1)
            x1 = r1(1)*2.D0 - 1.D0
            x2 = r2(1)*2.D0 - 1.D0
        end do
        spherical_unit_vectors(1,i) = 2.D0*x1*sqrt(1.D0-x1**2-x2**2)
        spherical_unit_vectors(2,i) = 2.D0*x2*sqrt(1.D0-x1**2-x2**2)
        spherical_unit_vectors(3,i) = 1.D0-2.D0*(x1**2 + x2**2)

    end do

end function spherical_unit_vectors

end program FF_HI_semimp
