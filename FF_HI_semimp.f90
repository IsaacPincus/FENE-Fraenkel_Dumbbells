program FF_HI_semimp
implicit none
real*8, parameter :: PI = 4*atan(1.0D0)
integer :: Nblocks, Nsteps, steps, block, time(1:8), kseed, i
real*8 :: dt, sr, b, h, a, Q0, Ql
real*8 :: Aeta, Veta, Apsi, Vpsi, Beta, Bpsi, Apsi2, Vpsi2, Bpsi2
real*8, dimension(3,3) :: k, delT, tau
real*8, dimension(3) :: Q, F
integer, dimension(:), allocatable :: seed

call date_and_time(values=time)
call random_seed(size=kseed)
allocate(seed(1:kseed))
forall (i=1:kseed) seed(i) =time(8)*100000+time(7)*1000+time(6)*10+time(5)+12341293*real(i)
call random_seed(put=seed)

Nblocks = 5000
Nsteps = 1000000
dt = 0.1D0
sr = 1.D0
b = 10000000D0
h = 0.15D0
a = h*sqrt(PI)
Q0 = 0

k(:,:) = 0.D0
k(1,2) = 1.D0
k = sr*k

tau(:,:) = 0.D0

delT(:,:) = 0.D0
delT(1,1) = 1.D0
delT(2,2) = 1.D0
delT(3,3) = 1.D0

Aeta = 0.D0
Apsi = 0.D0
Apsi2 = 0.D0
Veta = 0.D0
Vpsi = 0.D0
Vpsi2 = 0.D0

Q(:) = (/1,1,1/)

mainloop: do block=1,Nblocks

    tau(:,:) = 0.D0

    equilibration: do steps=1,100
        k(:,:) = 0.D0
        call step(Q, k, dt, Q0, delT, b, a)
    end do equilibration

    do steps=1,Nsteps
        k(1,2) = sr
        call step(Q, k, dt, Q0, delT, b, a)
        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:)/Ql
        tau(:,:) = tau(:,:) + dyadic_prod(Q,F)

    end do

    Beta = tau(1,2)/(Nsteps*sr)
    Bpsi = (tau(1,1) - tau(2,2))/(Nsteps*sr**2)
    Bpsi2 = (tau(2,2) - tau(3,3))/(Nsteps*sr**2)
    Aeta = Aeta + Beta
    Apsi = Apsi + Bpsi
    Apsi2 = Apsi2 + Bpsi2
    Veta = Veta + Beta**2
    Vpsi = Vpsi + Bpsi**2
    Vpsi2 = Vpsi2 + Bpsi2**2

    if (mod(block,100).eq.0) then
        write(*,*) "block = "
        write(*,*) block
    endif

end do mainloop

Aeta = Aeta/Nblocks
Veta = Veta/Nblocks
Veta = sqrt((Veta - Aeta**2)/(Nblocks-1))

Apsi = Apsi/Nblocks
Vpsi = Vpsi/Nblocks
Vpsi = sqrt((Vpsi - Apsi**2)/(Nblocks-1))

Apsi2 = Apsi2/Nblocks
Vpsi2 = Vpsi2/Nblocks
Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Nblocks-1))

write(*,*) Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2

contains

function construct_B(Q, a, delT)
    implicit none
    real*8, intent(in) :: Q(3), a
    real*8, intent(in) :: delT(3,3)
    real*8, dimension(3,3) :: construct_B
    real*8 :: Ql, aux, g, g_til, a2, a4, Ql2, Ql4, Ql6, C43

    Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
    Ql2 = Ql**2
    Ql4 = Ql**4
    Ql6 = Ql**6
    a2 = a**2
    a4 = a**4
    C43 = 4.D0/3.D0

    aux = a/(C43*Ql*(Ql2 + C43*a2)**3)
    g = 1.D0-aux*(Ql6 + 14.D0/3.D0*a2*Ql4 + 8.D0*a4*Ql2)
    g_til = -aux*(Ql6 + 2.D0*a2*Ql4 - 8.D0/3.D0*a4*Ql2)
    construct_B = sqrt(g)*delT + (sqrt(g+g_til)-sqrt(g))*dyadic_prod(Q,Q)/Ql2
    return

end function construct_B

function Wiener_step(dt)
    implicit none
    real*8, intent(in) :: dt
    real*8, dimension(3) :: Wiener_step
    real*8, dimension(3) :: dW

    call random_number(dW)
    dW = dW - 0.5D0

    Wiener_step = dW*(14.14855378D0*sqrt(dt)*dW*dW + 1.21569221D0*sqrt(dt))
    return

end function Wiener_step

subroutine step(Q, k, dt, Q0, delT, b, a)
    implicit none
    real*8, intent(in) :: dt, Q0, b, a
    real*8, intent(in) :: delT(3,3), k(3,3)
    real*8 :: Ql, L, RdotB2_L, Qlength, temp_R1, temp_R2
    real*8, dimension(3) :: F, dW, u, RHS, Qpred, RdotB2
    real*8, dimension(3,3) :: B1, B2, B1B1, B2B2
    real*8 :: Q(3)

    Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
    F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q/Ql

    B1 = construct_B(Q, a, delT)
    dW = Wiener_step(dt)

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

    Q = RHS*Qlength/L

end subroutine step

function find_roots(a, b, c, Q0, lim)
    implicit none
    real*8, parameter :: PI = 4*atan(1.0D0)
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
        endif
   end do

end function find_roots

function dyadic_prod(vec1,vec2)
    implicit none
    real*8, dimension(3), intent(in) :: vec1, vec2
    real*8, dimension(3,3) :: dyadic_prod
    integer :: i, j

    do j=1,3
        do i=1,3
            dyadic_prod(i,j) = vec1(i)*vec2(j)
        end do
    end do
end function dyadic_prod

function ten_vec_dot(tensor, vector)
    implicit none
    real*8, dimension(3,3), intent(in) :: tensor
    real*8, dimension(3), intent(in) :: vector
    real*8, dimension(3) :: ten_vec_dot
    integer :: i

    do i=1,3
        ten_vec_dot(i) = dot_product(tensor(i,:), vector)
    end do
end function ten_vec_dot

end program FF_HI_semimp
