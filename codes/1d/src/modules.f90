module modules
  implicit none

  integer,parameter :: n=100 !< 設計変数の数
  real(8) :: ca !< Armijo条件の定数
  real(8) :: cw !< Wolfe条件の定数

  real(8) :: yng=1.d0 !< ヤング率
  real(8) :: lngth=2.d0/dble(n) !< 各セグメントの長さ

  real(8) :: vol_max !< 体積制約の上限
  
  real(8) :: kmat(n,n) !< 剛性行列
  real(8) :: pvec(n) !< 右辺ベクトル
  real(8) :: u(n) !< 変位
  
  real(8) :: dkmat(n,n,n) !< 剛性行列の微分
  real(8) :: eps !< 直線探索パラメタ(アクセル=1/(1-eps), ブレーキ=1/(1+eps))
  
  integer :: ipiv(n), info !< for Lapack

contains

  !> 体積を計算
  real(8) function cal_vol(x) result(out)
    real(8),intent(in) :: x(n)
    integer :: i
    out=0.d0
    do i=1,n
       out=out+x(i)*lngth
    end do
  end function cal_vol

  !> 初期設定
  subroutine init(x,vol,alp,rho)
    real(8),intent(out) :: alp, vol, x(n),rho

    alp=0.1d0 !< 初期ステップサイズ
    eps=0.01d0 !< 直線探索パラメタ
    rho=3.d0 !< 初期ペナルティを設定
    ca=0.01d0 !< Armijo constantを設定
    cw=0.02d0 !< Wolfe constantを設定
    pvec(:)=1.d0 !< 境界条件(荷重)を設定
    x(:)=1.d0 !< 初期条件を設定
    vol=cal_vol(x) !< 初期体積を計算
    vol_max=vol*0.5d0 !< 体積制約を設定

  end subroutine init

  !> 剛性行列を計算
  subroutine set_kmat(x)
    real(8),intent(in) :: x(n)

    integer :: i
    kmat(:,:)=0.d0
    ! diagonal elements
    do i=1,n-1
       kmat(i,i)=x(i)+x(i+1)
    end do
    kmat(n,n)=x(n)

    ! off-diagonal elements
    do i=1,n-1
       kmat(i,i+1)=-x(i+1)
    end do
    
    do i=1,n-1
       kmat(i+1,i)=-x(i+1)
    end do

    kmat(:,:)=yng/lngth*kmat(:,:)

  end subroutine set_kmat

  !> 右辺ベクトル
  subroutine set_rhs(x)
    real(8),intent(in) :: x(n)
    u(:)=pvec(:)
  end subroutine set_rhs

  real(8) function objective_function(rho,x) result(out)
    real(8) :: rho, x(n), vol
    integer :: i
    
    call set_kmat(x)
    call set_rhs(x)
    call dgesv(n,1,kmat,n,ipiv,u,n,info)
    vol=cal_vol(x)
    out=dot_product(u,pvec)
    if(vol-vol_max>0)then
       out=out+rho*(vol-vol_max)**2
    end if
    
  end function objective_function

  !> 感度解析
  subroutine sensitivity(x,rho,sens)
    real(8),intent(in) :: x(n), rho
    real(8),intent(out) :: sens(n)
    integer :: i
    real(8) :: tmp(n), vol
    
    dkmat(:,:,:)=0.d0
    dkmat(1,1,1)=1.d0
    do i=2,n
       dkmat(i-1,i-1,i)=+1.d0
       dkmat(i  ,i-1,i)=-1.d0
       dkmat(i-1,i  ,i)=-1.d0
       dkmat(i  ,i  ,i)=+1.d0
    end do

    dkmat(:,:,:)=yng/lngth*dkmat(:,:,:)

    do i=1,n
       tmp=matmul(dkmat(:,:,i),u)
       sens(i)=-dot_product(u,tmp)
    end do
    vol=cal_vol(x)
    if(vol-vol_max>0)then
       sens(:)=sens(:)+rho*2.d0*(vol-vol_max)*lngth
    end if
    
  end subroutine sensitivity

  !> armijoの条件
  subroutine armijo(x,alp,fobj,sens,d,rho)
    real(8),intent(in) :: x(n), fobj, sens(n), d(n), rho
    real(8),intent(inout) :: alp
    real(8) :: fobj_new

    fobj_new=objective_function(rho,x(:)+alp*d(:))
    do while(fobj_new>fobj+ca*alp*dot_product(sens,d))
       alp=alp/(1.d0+eps)
       fobj_new=objective_function(rho,x(:)+alp*d(:))
    end do
    
  end subroutine armijo

  !> wolfeの条件
  logical function wolfe(x,alp,fobj,sens,d,rho) result(out)
    real(8),intent(in) :: x(n), fobj, sens(n), d(n), rho
    real(8),intent(inout) :: alp

    real(8) :: sens_new(n)

    call sensitivity(x(:)+alp*d(:),rho,sens_new)
    if(dot_product(sens_new,d).ge.cw*dot_product(sens,d))then
       out=.true.
    else
       out=.false.
    end if
    
  end function wolfe
    
  
  !> 直線探索
  subroutine linear_search(x,alp,rho,fobj,sens,d)
    real(8),intent(in) :: fobj, rho, sens(n), d(n)
    real(8), intent(inout) :: x(n), alp
    logical :: wolfe_cond_is_satisfied
    integer :: i
    
    call armijo(x,alp,fobj,sens,d,rho)
    wolfe_cond_is_satisfied=wolfe(x,alp,fobj,sens,d,rho)
    if(.not.wolfe_cond_is_satisfied)then
       alp=alp/(1.d0-eps)
    end if

    ! 側面制約
    do i=1,n
       if(x(i).ge.3.d0) x(i)=3.d0-epsilon(1.d0)
       if(x(i).le.0.d0) x(i)=epsilon(1.d0)
    end do

    x(:)=x(:)+alp*d(:)

  end subroutine linear_search
  
end module modules
