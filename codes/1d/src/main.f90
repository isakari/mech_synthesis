program main
  use modules
  implicit none

  integer :: i, j
  
  real(8) :: x(n) !< 設計変数
  real(8) :: fobj !< 目的関数
  real(8) :: sens(n) !< 感度
  real(8) :: alp !< ステップサイズ
  real(8) :: d(n) !< 探索方向ベクトル
  real(8) :: vol !< 体積
  real(8) :: rho !< 体積制約のペナルティ
  
  call init(x,vol,alp,rho)
  do i=1,300000
     if(alp.le.epsilon(1.d0)) stop
     fobj=objective_function(rho,x)
     write(1,*) i, fobj
     write(2,*) x(:)
     do j=1,n
        write(10000+i,*) (j-1)*lngth, -x(j)/2.d0
        write(10000+i,*) j*lngth, -x(j)/2.d0
        write(10000+i,*) j*lngth, x(j)/2.d0
        write(10000+i,*) (j-1)*lngth, x(j)/2.d0
        write(10000+i,*) (j-1)*lngth, -x(j)/2.d0     
        write(10000+i,*)
        write(10000+i,*)
     end do
     close(10000+i)
     call sensitivity(x,rho,sens)
     d(:)=-sens(:)  !最急降下法
     call linear_search(x,alp,rho,fobj,sens,d)
     rho=rho*(1.d0+eps)
  end do
    
end program main
