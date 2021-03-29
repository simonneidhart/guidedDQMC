==============================================================================!
!rmsdmin
!  This subroutine calculates the least square rmsd of two coordinate
!  sets coord1(3,n) and coord2(3,n) using a method based on quaternion.
! adapted from ls_rmsd module of Chaok Seok and Evangelos Coutsias 2004
!==============================================================================!
subroutine rmsdmin(nat, coord1, coord2, rxyzrot, rxyz2, rmsd)
!   use metric_rmsd_ls_rmsd
   implicit none
   !parameters
   real*8 :: coord1(3,nat),coord2(3,nat)
   real*8 :: rmsd,xa,ya,za,xb,yb,zb
   !internal
   integer :: iat, nat
   real*8 :: rxyz1(3,nat), rxyz2(3,nat), rxyzrot(3,nat)
   real*8 :: Rmatrix(3,3), S(4,4)
   real*8 :: q(4)
!   integer,parameter :: lwork=104 !104 is optimal
!   integer :: info
!   real*8 :: work(lwork)
   real*8 :: U(3,3)
   real*8 :: lambda
!   real*8 ::x_norm,y_norm

   ! make copies of the original coordinates with center at origin
        xa=0.d0 ; ya=0.d0 ; za=0.d0
        xb=0.d0 ; yb=0.d0 ; zb=0.d0
        do iat=1,nat
        rxyz1(1,iat) = coord1(1,iat)
        rxyz1(2,iat) = coord1(2,iat)
        rxyz1(3,iat) = coord1(3,iat)
        rxyz2(1,iat) = coord2(1,iat)
        rxyz2(2,iat) = coord2(2,iat)
        rxyz2(3,iat) = coord2(3,iat)
        xa=xa+rxyz1(1,iat)
        ya=ya+rxyz1(2,iat)
        za=za+rxyz1(3,iat)
        xb=xb+rxyz2(1,iat)
        yb=yb+rxyz2(2,iat)
        zb=zb+rxyz2(3,iat)
        enddo

        xa=xa/nat ; ya=ya/nat ; za=za/nat
        xb=xb/nat ; yb=yb/nat ; zb=zb/nat
        do iat=1,nat
        rxyz1(1,iat) = rxyz1(1,iat)-xa
        rxyz1(2,iat) = rxyz1(2,iat)-ya
        rxyz1(3,iat) = rxyz1(3,iat)-za
        rxyz2(1,iat) = rxyz2(1,iat)-xb
        rxyz2(2,iat) = rxyz2(2,iat)-yb
        rxyz2(3,iat) = rxyz2(3,iat)-zb
        enddo

!!  ! calculate the R matrix
!!!  x_norm = 0.0d0
!!!  y_norm = 0.0d0
!!   do i = 1, 3
!!!     x_norm = x_norm + dot_product(rxyz1(i,:), rxyz1(i,:))
!!!     y_norm = y_norm + dot_product(rxyz2(i,:), rxyz2(i,:))
!!      do j = 1, 3
!!         Rmatrix(i,j) = dot_product(rxyz1(i,:),rxyz2(j,:))
!!      end do
!!   end do

   do iat=1,nat
   Rmatrix(1,1)=Rmatrix(1,1)+rxyz1(1,iat)*rxyz2(1,iat)
   Rmatrix(2,1)=Rmatrix(2,1)+rxyz1(2,iat)*rxyz2(1,iat)
   Rmatrix(3,1)=Rmatrix(3,1)+rxyz1(3,iat)*rxyz2(1,iat)
   Rmatrix(1,2)=Rmatrix(1,2)+rxyz1(1,iat)*rxyz2(2,iat)
   Rmatrix(2,2)=Rmatrix(2,2)+rxyz1(2,iat)*rxyz2(2,iat)
   Rmatrix(3,2)=Rmatrix(3,2)+rxyz1(3,iat)*rxyz2(2,iat)
   Rmatrix(1,3)=Rmatrix(1,3)+rxyz1(1,iat)*rxyz2(3,iat)
   Rmatrix(2,3)=Rmatrix(2,3)+rxyz1(2,iat)*rxyz2(3,iat)
   Rmatrix(3,3)=Rmatrix(3,3)+rxyz1(3,iat)*rxyz2(3,iat)
   enddo


   ! S matrix
   S(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
   S(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2)
   S(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3)
   S(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1)

   S(1, 2) = S(2, 1)
   S(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
   S(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1)
   S(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1)

   S(1, 3) = S(3, 1)
   S(2, 3) = S(3, 2)
   S(3, 3) =-Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
   S(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2)

   S(1, 4) = S(4, 1)
   S(2, 4) = S(4, 2)
   S(3, 4) = S(4, 3)
   S(4, 4) =-Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3)

   ! Calculate eigenvalues and eigenvectors, and
   ! take the maximum eigenvalue lambda and the corresponding eigenvector q.
   call dstmev(S, lambda, q)
!   call DSYEV('V','U',4,S,4,q(1),work(1),lwork,info)
!   lambda=q(4)

!   call rotation_matrix(S(:,4), U)
   call rotation_matrix(q, U)
   do iat=1,nat
      rxyzrot(1:3,iat) = matmul(U, rxyz1(1:3,iat) )  ! rotate
   enddo
!using the eigenvalue lambda for caluclating rmsd2
!fail for some configurations
!rmsd2out=max(0.0d0,(x_norm+y_norm)-2.0d0*lambda) ! RMSD=sqrt(error2/n)
!the follwing is not much slower:
!   rmsd2out = sum((coord1-coord2)**2)
        rmsd=0.d0
        do iat=1,nat
        rmsd=rmsd+(rxyz2(1,iat)- rxyzrot(1,iat))**2 +  (rxyz2(2,iat)- rxyzrot(2,iat))**2 + (rxyz2(3,iat)- rxyzrot(3,iat))**2
        enddo
        rmsd=sqrt(rmsd/nat)
end subroutine




!-----------------------------------------------------------------------
subroutine rotation_matrix(q, U)
!-----------------------------------------------------------------------
! This subroutine constructs rotation matrix U from quaternion q.
!-----------------------------------------------------------------------

  real*8, dimension(4) :: q
  real*8, dimension(3,3) :: U
  real*8::q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33

  q0 = q(1)
  q1 = q(2)
  q2 = q(3)
  q3 = q(4)

  b0 = 2.0d0*q0
  b1 = 2.0d0*q1
  b2 = 2.0d0*q2
  b3 = 2.0d0*q3

  q00 = b0*q0-1.0d0
  q01 = b0*q1
  q02 = b0*q2
  q03 = b0*q3

  q11 = b1*q1
  q12 = b1*q2
  q13 = b1*q3

  q22 = b2*q2
  q23 = b2*q3

  q33 = b3*q3

  U(1,1) = q00+q11
  U(1,2) = q12-q03
  U(1,3) = q13+q02

  U(2,1) = q12+q03
  U(2,2) = q00+q22
  U(2,3) = q23-q01

  U(3,1) = q13-q02
  U(3,2) = q23+q01
  U(3,3) = q00+q33

end subroutine rotation_matrix

!-----------------------------------------------------------------------
subroutine DSTMEV(A,lambda,evec)
!-----------------------------------------------------------------------
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through
!     multiplication by V'.
!-----------------------------------------------------------------------
  real*8, dimension(4,4) :: A, T, V, SV
  integer :: i
  integer, dimension(1) :: max_loc ! must be an array
  real*8, dimension(4) :: evec, SW
  real*8, dimension(8) :: rv1
  real*8 :: lambda

  !-----------------------------------------------------------------------
  !(I).   Convert to tridiagonal form, keeping similarity transform
  ! (a product of 3 Givens rotations)
  call givens4(A,T,V)

  !-----------------------------------------------------------------------
  !(II)Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
  lambda=min(T(1,1)-abs(T(1,2)),-abs(T(2,1))+T(2,2)-abs(T(2,3)),&
       -abs(T(3,2))+T(3,3)-abs(T(3,4)),-abs(T(4,3))+T(4,4))
  !-----------------------------------------------------------------------
  !(III). Form positive definite matrix     T <== lambda*I - T
  do i = 1,4
     T(i,i) = T(i,i)-lambda
  enddo
  !-----------------------------------------------------------------------
  !(IV). Compute singular values/vectors of SPD matrix B
  call svdcmp(4,T,4,4,SW,SV,rv1)
  !-----------------------------------------------------------------------
  !(V). Shift spectrum back
  max_loc = maxloc(SW)
  lambda = SW(max_loc(1)) + lambda
  !lambda = SW(1) + lambda
  !-----------------------------------------------------------------------
  !(VI). Convert eigenvector to original coordinates: (V is transposed!)
  evec = matmul(V,SV(:,max_loc(1)))
  !write(*,*)'-----------------------------------------------------------'
  !write(*,*) 'lambda = ', lambda,'  eigenvector:  '
  !write(*,99) evec
  !write(*,*)'-----------------------------------------------------------'
  !-----------------------------------------------------------------------
  !99 format(1x,4(d19.13,1x))

end subroutine dstmev

!-----------------------------------------------------------------------
subroutine givens4(S,T,V)
!-----------------------------------------------------------------------
  real*8, dimension(4,4) :: S
  real*8, dimension(4,4) :: T,V
  real*8  :: c1,c2,c3, s1,s2,s3, r1,r2,r3, c1c2, s1c2
  real*8 :: pythag
  ! external        pythag
  !performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
  T=S; V = 0.d0
  !-----------------------------------------------------------------------
  !Zero out entries T(4,1) and T(1,4)
  ! compute cos and sin of rotation angle in the 3-4 plane
  r1 = pythag(T(3,1),T(4,1))
  if(r1 .ne. 0.d0) then
     c1 = T(3,1)/r1; s1 = T(4,1)/r1
     V(3,3) = c1  ; V(3,4) = s1
     V(4,3) =-s1  ; V(4,4) = c1
     T(3,1) = r1  ; T(4,1) = 0.d0
     T(3:4,2:4) = matmul(V(3:4,3:4),T(3:4,2:4))
     T(1:2,3:4) = transpose(T(3:4,1:2))
     T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
  else
     c1 = 1.d0; s1 = 0.d0
  endif
  !-----------------------------------------------------------------------
  !Zero out entries T(3,1) and T(1,3)
  ! compute cos and sin of rotation angle in the 2-3 plane
  r2 = pythag(T(3,1), T(2,1))
  if(r2 .ne. 0.d0) then
     c2 = T(2,1)/r2; s2 = T(3,1)/r2
     V(2,2) = c2  ; V(2,3) = s2
     V(3,2) =-s2  ; V(3,3) = c2
     T(2,1) = r2  ; T(3,1) = 0.d0
     T(2:3,2:4) = matmul(V(2:3,2:3),T(2:3,2:4))
     T(1,2:3)   = T(2:3,1);  T(4,2:3) = T(2:3,4)
     T(2:3,2:3) = matmul(T(2:3,2:3),transpose(V(2:3,2:3)))
  else
     c2 = 1.d0; s2 = 0.d0
  endif
  !-----------------------------------------------------------------------
  !Zero out entries T(4,2) and T(2,4)
  ! compute cos and sin of rotation angle in the 3-4 plane
  r3 = pythag(T(4,2), T(3,2))
  if(r3 .ne. 0.d0) then
     c3 = T(3,2)/r3; s3 = T(4,2)/r3
     V(3,3) = c3  ; V(3,4) = s3
     V(4,3) =-s3  ; V(4,4) = c3
     T(3,2) = r3  ; T(4,2) = 0.d0
     T(3:4,3:4) = matmul(V(3:4,3:4),T(3:4,3:4))
     T(1:2,3:4) = transpose(T(3:4,1:2))
     T(3:4,3:4) = matmul(T(3:4,3:4),transpose(V(3:4,3:4)))
  else
     c3 = 1.d0; s3 = 0.d0
  endif
  !-----------------------------------------------------------------------
  !Compute net rotation matrix (accumulate similarity for evec. computation)
  ! To save transposing later, This is the transpose!
  V(1,1)=1.d0; V(1,2:4) = 0.d0; V(2:4,1) = 0.d0
  V(2,2) = c2;  V(3,2) = c1*s2 ; V(4,2) = s1*s2; c1c2 = c1*c2; s1c2=s1*c2
  V(2,3) = -s2*c3 ; V(3,3) = c1c2*c3-s1*s3 ; V(4,3) =  s1c2*c3+c1*s3
  V(2,4) =  s2*s3 ; V(3,4) =-c1c2*s3-s1*c3 ; V(4,4) = -s1c2*s3+c1*c3
  !-----------------------------------------------------------------------
  !write(*,*) (V(1:4,i) - W(1:4,i),i=1,4)
end subroutine givens4


!-----------------------------------------------------------------------
SUBROUTINE svdcmp(mmax,a,m,n,w,v,rv1)
!-----------------------------------------------------------------------
  integer :: mmax
  INTEGER :: m,n
  real*8 :: a(mmax,*),v(mmax,*),w(*),rv1(*)
  INTEGER :: i,its,j,jj,k,l,nm
  real*8 :: anorm,c,f,g,h,s,scale,x,y,z,pythag

  g = 0.0d0
  scale = 0.0d0
  anorm = 0.0d0
  do i = 1, n
     l=i+1
     rv1(i)=scale*g
     g=0.0d0
     s=0.0d0
     scale=0.0d0
     if(i.le.m)then
        do k=i,m
           scale=scale+abs(a(k,i))
        end do
        if(scale.ne.0.0d0)then
           do k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
           end do
           f=a(i,i)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,i)=f-g
           do j=l,n
              s=0.0d0
              do k=i,m
                 s=s+a(k,i)*a(k,j)
              end do
              f=s/h
              do k=i,m
                 a(k,j)=a(k,j)+f*a(k,i)
              end do
           end do
           do k=i,m
              a(k,i)=scale*a(k,i)
           end do
        endif
     endif
     w(i)=scale *g
     g=0.0d0
     s=0.0d0
     scale=0.0d0
     if((i.le.m).and.(i.ne.n))then
        do k=l,n
           scale=scale+abs(a(i,k))
        end do
        if(scale.ne.0.0d0)then
           do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
           end do
           f=a(i,l)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,l)=f-g
           do k=l,n
              rv1(k)=a(i,k)/h
           end do
           do j=l,m
              s=0.0d0
              do k=l,n
                 s=s+a(j,k)*a(i,k)
              end do
              do k=l,n
                 a(j,k)=a(j,k)+s*rv1(k)
              end do
           end do
           do k=l,n
              a(i,k)=scale*a(i,k)
           end do
        endif
     endif
     anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
  end do

  do i = n, 1, -1
     if(i .lt. n) then
        if(g.ne.0.0d0)then
           do j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
           end do
           do j=l,n
              s=0.0d0
              do k=l,n
                 s=s+a(i,k)*v(k,j)
              end do
              do k=l,n
                 v(k,j)=v(k,j)+s*v(k,i)
              end do
           end do
        endif
        do j=l,n
           v(i,j)=0.0d0
           v(j,i)=0.0d0
        end do
     endif
     v(i,i)=1.0d0
     g=rv1(i)
     l=i
  end do
!-----------------------------------------------------------------------
  do i = min(m,n), 1, -1
     l=i+1
     g=w(i)
     do j=l,n
        a(i,j)=0.0d0
     end do
     if(g.ne.0.0d0)then
        g=1.0d0/g
        do j=l,n
           s=0.0d0
           do k=l,m
              s=s+a(k,i)*a(k,j)
           end do
           f=(s/a(i,i))*g
           do k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
           end do
        end do
        do j=i,m
           a(j,i)=a(j,i)*g
        end do
     else
        do j= i,m
           a(j,i)=0.0d0
        end do
     endif
     a(i,i)=a(i,i)+1.0d0
  end do
!-----------------------------------------------------------------------
  do k=n,1,-1
     do its=1,30
        do l=k,1,-1
           nm=l-1
           if((abs(rv1(l))+anorm).eq.anorm)  goto 2
           if((abs(w(nm))+anorm).eq.anorm)  goto 1
        end do
1       c=0.0d0
        s=1.0d0
        do i=l,k
           f=s*rv1(i)
           rv1(i)=c*rv1(i)
           if((abs(f)+anorm).eq.anorm) goto 2
           g=w(i)
           h=pythag(f,g)
           w(i)=h
           h=1.0d0/h
           c= (g*h)
           s=-(f*h)
           do j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
           end do
        end do
2       z=w(k)
        if(l .eq. k)then
           if(z.lt.0.0d0)then
              w(k)=-z
              do j=1,n
                 v(j,k)=-v(j,k)
              end do
           endif
           goto 3
        endif
        if(its.eq.30) then
           write(*,*) 'no convergence in svdcmp'
           stop
        endif
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
        g=pythag(f,1.0d0)
        f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c=1.0d0
        s=1.0d0
        do j=l,nm
           i=j+1
           g=rv1(i)
           y=w(i)
           h=s*g
           g=c*g
           z=pythag(f,h)
           rv1(j)=z
           c=f/z
           s=h/z
           f= (x*c)+(g*s)
           g=-(x*s)+(g*c)
           h=y*s
           y=y*c
           do jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
           end do
           z=pythag(f,h)
           w(j)=z
           if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
           endif
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           do jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
           end do
        end do
        rv1(l)=0.0d0
        rv1(k)=f
        w(k)=x
     end do
3    continue
  end do
!-----------------------------------------------------------------------

END subroutine svdcmp



!-----------------------------------------------------------------------
double precision FUNCTION pythag(a,b)
!-----------------------------------------------------------------------

  real*8 ::  a, b
  real*8 :: absa,absb

  absa=abs(a)
  absb=abs(b)
  if(absa.gt.absb)then
     pythag=absa*dsqrt(1.0d0+(absb/absa)**2)
  else
     if(absb.eq.0.0d0)then
        pythag=0.0d0
     else
        pythag=absb*dsqrt(1.0d0+(absa/absb)**2)
     endif
  endif

END function pythag
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


  subroutine readfiles(natx,nat,rxyz,alat,symb)
  implicit real*8 (a-h,o-z)
  dimension rxyz(3, natx),alat(3, 3)
  character(len=2), dimension(natx) :: symb
  character(len=20) :: ustring
   logical :: debug

   !debug=.true.
   debug=.false.



        read(10,*) nat ,ustring
        if (debug) write(*,*)  'nat,ustring=',nat,ustring
        if (nat.gt.natx) stop 'natx'
        read(10,*)


       if (trim(ustring).eq.'angstroemd0' .or. trim(ustring).eq.'angstroem') then
               convert=1.d0/0.52917720859d0
               write(*,*) 'The units are assumed as angstroem'
       else
               convert=1.d0
               !write(*,*) 'The units are assumed as Bohr!!!'
       endif
               !write(*,*) 'The units are assumed as angstroem'
               !convert=1.d0/0.52917720859d0
           if (debug) write(*,*)  ' convert ',convert
        do i = 1, nat
           read(10, *) symb(i),( rxyz(j, i), j = 1, 3)
           if (debug) write(*,*)  symb(i),( rxyz(j, i), j = 1, 3)
           do l=1,3
           rxyz(l,i)=(rxyz(l,i))*convert
           enddo
        end do

        do j=1,3
        do i=1,3
           alat(i,j) = 0.d0
        enddo
        enddo
        !write(*,*) 'WARNING: assuming free boundary conditions'


        end subroutine
