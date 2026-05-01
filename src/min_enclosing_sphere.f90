! Implementation of the smallest enclosing sphere Welzl algorithm (equivalent to Python miniball library)
subroutine minimum_enclosing_sphere(P, N, C, R)
    USE iso_fortran_env
    Implicit NONE 
    
    integer, intent(in) :: N
    real(8), intent(inout) :: P(N,3)
    real(8), intent(out) :: C(3), R

    call random_shuffle(P, N)
    call welzl(P, N, 0, C, R)
end subroutine

recursive subroutine welzl(P, n, m, C, R)
    real(8), intent(in) :: P(n,3)
    integer, intent(in) :: n, m
    real(8), intent(out) :: C(3), R
    real(8) :: Ctmp(3), Rtmp

    if (n == 0 .or. m == 4) then
        call sphere_from_points(P(1:m,:), m, C, R)
        return
    end if

    call welzl(P, n-1, m, C, R)

    if (distance(P(n,:), C) > R + EPS) then
        call welzl([P(1:n-1,1:3), P(n,1:3)], n-1, m+1, C, R)
    end if
end subroutine

subroutine sphere_from_points(P, m, C, R)
    real(8), intent(in) :: P(:,:)
    integer, intent(in) :: m
    real(8), intent(out) :: C(3), R

    select case (m)
    case (0)
        C = 0.0d0
        R = 0.0d0
    case (1)
        C = P(1,:)
        R = 0.0d0
    case (2)
        C = 0.5d0 * (P(1,:) + P(2,:))
        R = distance(P(1,:), C)
    case (3)
        call circumcenter3(P(1,:), P(2,:), P(3,:), C)
        R = distance(P(1,:), C)
    case (4)
        call circumsphere4(P, C)
        R = distance(P(1,:), C)
    end select
end subroutine

subroutine circumcenter3(A, B, Cc, C)
    real(8), intent(in) :: A(3), B(3), Cc(3)
    real(8), intent(out) :: C(3)
    real(8) :: aloc(3), bloc(3), axb(3), denom

    aloc = B - A
    bloc = Cc - A
    axb = cross(aloc, bloc)
    denom = 2.0d0 * dot(axb, axb)

    C = A + (cross(axb, a) * dot(b, b) + &
             cross(b, axb) * dot(a, a)) / denom
end subroutine

subroutine circumsphere4(P, C)
    real(8), intent(in) :: P(4,3)
    real(8), intent(out) :: C(3)
    real(8) :: A(4,4), Dx, Dy, Dz, D

    A(:,1:3) = P
    A(:,4) = 1.0d0

    Dx = det4(A(:,[2,3,4,1]))
    Dy = -det4(A(:,[1,3,4,2]))
    Dz = det4(A(:,[1,2,4,3]))
    D  = -det4(A(:,1:4))

    C = (/ Dx, Dy, Dz /) / (2.0d0 * D)
end subroutine

function distance(a, b) result(d)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: d
    d = sqrt(sum((a - b)**2))
end function

function dot(a, b) result(d)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: d
    d = sum(a*b)
end function

function cross(a, b) result(c)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: c(3)
    c = (/ a(2)*b(3)-a(3)*b(2), &
           a(3)*b(1)-a(1)*b(3), &
           a(1)*b(2)-a(2)*b(1) /)
end function

subroutine random_shuffle(P, N)
    real(8), intent(inout) :: P(N,3)
    integer :: i, j
    real(8) :: tmp(3), r

    do i = N, 2, -1
        call random_number(r)
        j = 1 + int(r*i)
        tmp = P(i,:)
        P(i,:) = P(j,:)
        P(j,:) = tmp
    end do
end subroutine

function det4(M) result(d)
    real(8), intent(in) :: M(4,4)
    real(8) :: d

    d = M(1,1)*det3(M(2:4,2:4)) - M(1,2)*det3(M(2:4,[1,3,4])) + &
        M(1,3)*det3(M(2:4,[1,2,4])) - M(1,4)*det3(M(2:4,1:3))
end function

function det3(M) result(d)
    real(8), intent(in) :: M(3,3)
    real(8) :: d

    d = M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2)) - &
        M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1)) + &
        M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
end function

