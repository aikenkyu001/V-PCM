! vpcm_fortran.f90
! Version strictly matching Python output format and content
module vpcm_mod
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  
  type :: VPCM_type
     integer :: n
     real(dp) :: lam, noise
     real(dp), allocatable :: A(:,:), K_core(:,:), K_fluct(:,:), K_prev(:,:), H(:,:)
     integer :: t, t_sb
     logical :: gauge_broken
  contains
     procedure :: init => vpcm_init
     procedure :: curvature => vpcm_curvature
     procedure :: temporal_covariant_derivative => vpcm_temporal_covariant_derivative
     procedure :: apply_dynamic_gauge => vpcm_apply_dynamic_gauge
     procedure :: maybe_spontaneous_symmetry_breaking => vpcm_maybe_spontaneous_symmetry_breaking
     procedure :: update => vpcm_update
  end type VPCM_type

  interface
     subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
       import :: dp
       character :: jobu, jobvt
       integer :: m, n, lda, ldu, ldvt, lwork, info
       real(dp) :: a(lda,*), s(*), u(ldu,*), vt(ldvt,*), work(*)
     end subroutine dgesvd

     subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
       import :: dp
       character :: jobz, uplo
       integer :: n, lda, lwork, info
       real(dp) :: a(lda,*), w(*), work(*)
     end subroutine dsyev
  end interface

contains

  subroutine vpcm_init(this, n_in, lam_in, noise_in, t_sb_in)
    class(VPCM_type), intent(inout) :: this
    integer, intent(in) :: n_in, t_sb_in
    real(dp), intent(in) :: lam_in, noise_in
    integer :: i, j
    real(dp) :: sigma
    this%n = n_in; this%lam = lam_in; this%noise = noise_in; this%t_sb = t_sb_in
    this%t = 0; this%gauge_broken = .false.
    if (allocated(this%A)) deallocate(this%A)
    if (allocated(this%K_core)) deallocate(this%K_core)
    if (allocated(this%K_fluct)) deallocate(this%K_fluct)
    if (allocated(this%K_prev)) deallocate(this%K_prev)
    if (allocated(this%H)) deallocate(this%H)
    allocate(this%A(n_in,n_in), this%K_core(n_in,n_in), this%K_fluct(n_in,n_in), this%K_prev(n_in,n_in), this%H(n_in,n_in))
    sigma = 1.2_dp
    do i = 1, n_in; do j = 1, n_in
       this%A(i,j) = exp(-((real(i-j,dp))**2)/(2.0_dp*sigma**2))
    end do; end do
    do i = 1, n_in; this%A(i,:) = this%A(i,:) / sum(this%A(i,:)); end do
    call randn_matrix(this%K_core); this%K_core = 0.5_dp * this%K_core
    call randn_matrix(this%K_fluct); this%K_fluct = 0.1_dp * this%K_fluct
    this%K_prev = this%K_core
    this%H = 0.0_dp; do i = 1, n_in; this%H(i,i) = 1.0_dp; end do
  end subroutine vpcm_init

  subroutine vpcm_curvature(this, F)
    class(VPCM_type), intent(in) :: this
    real(dp), intent(out) :: F(this%n,this%n)
    real(dp) :: omega(this%n,this%n)
    integer :: i
    omega = this%K_core; do i = 1, this%n; omega(i,i) = 0.0_dp; end do
    F = matmul(omega, transpose(omega)) - matmul(transpose(omega), omega)
  end subroutine vpcm_curvature

  subroutine vpcm_temporal_covariant_derivative(this, nabla_tK)
    class(VPCM_type), intent(in) :: this
    real(dp), intent(out) :: nabla_tK(this%n,this%n)
    real(dp) :: dK(this%n,this%n), omega_t(this%n,this%n)
    integer :: i
    dK = this%K_core - this%K_prev; omega_t = dK
    do i = 1, this%n; omega_t(i,i) = 0.0_dp; end do
    nabla_tK = dK + matmul(omega_t, this%K_core) - matmul(this%K_core, omega_t)
  end subroutine vpcm_temporal_covariant_derivative

  subroutine vpcm_apply_dynamic_gauge(this)
    class(VPCM_type), intent(inout) :: this
    real(dp) :: Hloc(this%n,this%n)
    if (.not. this%gauge_broken) then
       call orthogonal_matrix(Hloc)
       this%K_core = matmul(Hloc, matmul(this%K_core, transpose(Hloc)))
       this%K_fluct = matmul(Hloc, matmul(this%K_fluct, transpose(Hloc)))
    end if
  end subroutine vpcm_apply_dynamic_gauge

  subroutine vpcm_maybe_spontaneous_symmetry_breaking(this)
    class(VPCM_type), intent(inout) :: this
    integer :: info, lwork
    real(dp), allocatable :: a(:,:), w(:), work(:)
    if ((.not. this%gauge_broken) .and. (this%t >= this%t_sb)) then
       this%gauge_broken = .true.
       allocate(a(this%n,this%n), w(this%n)); a = this%K_core
       lwork = 3*this%n - 1; allocate(work(lwork))
       call dsyev('V','U', this%n, a, this%n, w, work, lwork, info)
       if (info == 0) then
          this%H = a
          this%K_core = matmul(transpose(this%H), matmul(this%K_core, this%H))
          this%K_fluct = matmul(transpose(this%H), matmul(this%K_fluct, this%H))
       end if
       deallocate(a, w, work)
    end if
  end subroutine vpcm_maybe_spontaneous_symmetry_breaking

  subroutine vpcm_update(this, phi_external, nabla_tK)
    class(VPCM_type), intent(inout) :: this
    real(dp), intent(in), optional :: phi_external(:)
    real(dp), intent(out) :: nabla_tK(this%n,this%n)
    complex(dp) :: K(this%n,this%n), K_new(this%n,this%n), Omega(this%n,this%n), xi(this%n,this%n)
    real(dp) :: phi(this%n)
    integer :: i, j
    call this%apply_dynamic_gauge()
    call this%maybe_spontaneous_symmetry_breaking()
    this%K_prev = this%K_core
    K = cmplx(this%K_core, this%K_fluct, kind=dp)
    if (present(phi_external)) then; phi = phi_external; else; phi = 0.0_dp; end if
    do i = 1, this%n; do j = 1, this%n
       Omega(i,j) = cmplx(exp(-abs(real(i-j,dp))/3.0_dp)*tanh(phi(i)), 0.0_dp, kind=dp)
       xi(i,j) = cmplx(this%noise*randn(), this%noise*randn(), kind=dp)
    end do; end do
    K_new = matmul(this%A, K + (matmul(Omega,K) - matmul(K,Omega)) - this%lam*K + xi)
    K_new = matmul(K_new, transpose(this%A))
    this%K_core = real(K_new); this%K_fluct = aimag(K_new); this%t = this%t + 1
    call this%temporal_covariant_derivative(nabla_tK)
  end subroutine vpcm_update

  subroutine randn_matrix(mat)
    real(dp), intent(out) :: mat(:,:)
    integer :: i, j
    do i = 1, size(mat,1); do j = 1, size(mat,2); mat(i,j) = randn(); end do; end do
  end subroutine randn_matrix

  function randn() result(z)
    real(dp) :: z, u1, u2
    call random_number(u1); call random_number(u2)
    if (u1 <= 1.0e-12_dp) u1 = 1.0e-12_dp
    z = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*acos(-1.0_dp)*u2)
  end function randn

  subroutine orthogonal_matrix(Q)
    real(dp), intent(out) :: Q(:,:)
    integer :: n, j, k; real(dp) :: vnorm
    n = size(Q,1); call randn_matrix(Q)
    do j = 1, n
       do k = 1, j-1; Q(:,j) = Q(:,j) - dot_product(Q(:,k), Q(:,j)) * Q(:,k); end do
       vnorm = sqrt(dot_product(Q(:,j), Q(:,j)))
       if (vnorm > 0.0_dp) Q(:,j) = Q(:,j) / vnorm
    end do
  end subroutine orthogonal_matrix

  function rank_threshold(mat, eps) result(r)
    real(dp), intent(in) :: mat(:,:)
    real(dp), intent(in) :: eps
    integer :: r, m, n, lwork, info
    real(dp), allocatable :: a(:,:), s(:), u(:,:), vt(:,:), work(:)
    m = size(mat,1); n = size(mat,2)
    allocate(a(m,n), s(min(m,n)), u(m,m), vt(n,n)); a = mat
    lwork = max(3*min(m,n)+max(m,n), 5*min(m,n)); allocate(work(lwork))
    call dgesvd('A','A', m, n, a, m, s, u, m, vt, n, work, lwork, info)
    r = count(s > eps)
    deallocate(a,s,u,vt,work)
  end function rank_threshold

  function entropy_vec(phi) result(e)
    real(dp), intent(in) :: phi(:)
    real(dp) :: e, ssum, eps_val; real(dp), allocatable :: p(:)
    integer :: i; eps_val = 1.0e-8_dp; allocate(p(size(phi))); p = abs(phi)
    ssum = sum(p); if (ssum <= eps_val) then; e = 0.0_dp; deallocate(p); return; end if
    p = p / ssum; do i = 1, size(p); if (p(i) < eps_val) p(i) = eps_val; end do
    e = -sum(p * log(p)); deallocate(p)
  end function entropy_vec

  function connected_components_2d(mat, threshold) result(comps)
    real(dp), intent(in) :: mat(:,:)
    real(dp), intent(in) :: threshold
    integer :: comps, n, i, j, head, tail, x, y, nx, ny, k
    logical, allocatable :: visited(:,:)
    integer, allocatable :: qx(:), qy(:)
    integer :: dx(4) = (/1,-1,0,0/), dy(4) = (/0,0,1,-1/)
    n = size(mat,1); allocate(visited(n,n), qx(n*n), qy(n*n)); visited = .false.; comps = 0
    do i = 1, n; do j = 1, n
       if (abs(mat(i,j)) > threshold .and. .not. visited(i,j)) then
          comps = comps + 1; head = 1; tail = 1; qx(1) = i; qy(1) = j; visited(i,j) = .true.
          do while (head <= tail)
             x = qx(head); y = qy(head); head = head + 1
             do k = 1, 4; nx = x+dx(k); ny = y+dy(k)
                if (nx >= 1 .and. nx <= n .and. ny >= 1 .and. ny <= n) then
                   if (abs(mat(nx,ny)) > threshold .and. .not. visited(nx,ny)) then
                      tail = tail + 1; qx(tail) = nx; qy(tail) = ny; visited(nx,ny) = .true.
                   end if
                end if
             end do
          end do
       end if
    end do; end do; deallocate(visited,qx,qy)
  end function connected_components_2d

  function cosine_sim(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c, na, nb, eps
    eps = 1.0e-8_dp
    na = sqrt(dot_product(a,a)) + eps
    nb = sqrt(dot_product(b,b)) + eps
    c = dot_product(a,b) / (na*nb)
  end function cosine_sim

end module vpcm_mod

program main
  use vpcm_mod
  implicit none
  integer, parameter :: n = 16, steps = 400
  type(VPCM_type) :: vpcm, v, v2
  real(dp) :: ranks(steps), ents(steps), comps(steps), curv(steps), nabla_t(steps)
  real(dp) :: phi(n), nt(n,n), F(n,n), evolution_data(steps,6)
  integer :: t, i, j, k, c
  real(dp) :: sigma_pix, f_num, phi_onp(n), finals(20,n), sim_mat(20,20)
  real(dp) :: lams(5), noises(5), res(5,5)
  real(dp) :: centers(20,n), basin(20)
  integer :: cluster_sizes(20)
  logical :: assigned
  character(len=256) :: fmt_str

  ! summary_fortran.txt を初期化
  open(40, file="summary_fortran.txt", status="replace"); close(40)

  print *, "V‑PCM Fortran版シミュレーション開始 (Strict Python Sync)"
  open(40, file="summary_fortran.txt", position="append")
  write(40, '(A)') "V‑PCM 究極完全版シミュレーション開始"
  close(40)
  
  ! 1. 構造進化
  call vpcm%init(n, 0.05_dp, 0.02_dp, 100)
  do t = 1, steps
     phi = 0.0_dp; if (t <= 250) phi(1:8) = 1.0_dp
     call vpcm%update(phi, nt)
     ranks(t) = real(rank_threshold(vpcm%K_core, 1.0e-3_dp), dp)
     ents(t) = entropy_vec(diagonal(vpcm%K_core))
     comps(t) = real(connected_components_2d(vpcm%K_core, 0.2_dp), dp)
     call vpcm%curvature(F); curv(t) = sqrt(sum(F*F))
     nabla_t(t) = sqrt(sum(nt*nt))
     evolution_data(t,:) = (/real(t-1,dp), ranks(t), ents(t), comps(t), curv(t), nabla_t(t)/)
  end do
  open(10, file="evolution_metrics_fortran.csv", status="replace")
  write(10, '(A)') "# step,rank,entropy,components,curvature,nabla_t"
  do t = 1, steps
     write(10, '(F0.6,A,F0.6,A,F0.6,A,F0.6,A,F0.6,A,F0.6)') evolution_data(t,1),',',evolution_data(t,2),',',evolution_data(t,3),',',evolution_data(t,4),',',evolution_data(t,5),',',evolution_data(t,6)
  end do; close(10)

  ! 2. ONP 解析 (Basin計算追加)
  phi_onp = 0.0_dp; phi_onp(1:8) = 1.0_dp
  do i = 1, 20
     call v%init(n, 0.04_dp, 0.01_dp, 100)
     do t = 1, 300; call v%update(phi_onp, nt); end do
     finals(i,:) = diagonal(v%K_core)
  end do
  do i = 1, 20; do j = 1, 20; sim_mat(i,j) = cosine_sim(finals(i,:), finals(j,:)); end do; end do
  open(20, file="onp_similarity_matrix_fortran.csv", status="replace")
  write(20, '(A)') "# pattern similarity matrix"
  do i = 1, 20; do j = 1, 20
     if (j<20) then; write(20, '(F0.6,A)', advance='no') sim_mat(i,j), ','; else; write(20, '(F0.6)') sim_mat(i,j); end if
  end do; end do; close(20)

  ! Clustering for summary
  c = 0; cluster_sizes = 0; centers = 0.0_dp
  do i = 1, 20
     assigned = .false.
     do k = 1, c
        if (cosine_sim(finals(i,:), centers(k,:)) > 0.98_dp) then
           cluster_sizes(k) = cluster_sizes(k) + 1; assigned = .true.; exit
        end if
     end do
     if (.not. assigned) then
        c = c + 1; cluster_sizes(c) = 1; centers(c,:) = finals(i,:)
     end if
  end do
  do k = 1, c; basin(k) = real(cluster_sizes(k),dp) / 20.0_dp; end do

  open(40, file="summary_fortran.txt", position="append")
  write(40, '(A)') ""
  write(40, '(A)') "実験: ONP 解析結果"
  write(40, '(A,I0)') " - Attractors: ", c
  write(40, '(A)', advance='no') " - Basin Sizes: ["
  do k = 1, c
     write(40, '(F5.2)', advance='no') basin(k)
     if (k < c) write(40, '(A)', advance='no') ", "
  end do
  write(40, '(A)') "]"
  close(40)

  ! 3. 相図
  do i = 1, 5; lams(i) = 0.01_dp + (0.12_dp-0.01_dp)*real(i-1,dp)/4.0_dp; noises(i) = 0.0_dp + (0.06_dp-0.0_dp)*real(i-1,dp)/4.0_dp; end do
  do i = 1, 5; do j = 1, 5
     call v2%init(n, lams(i), noises(j), 100)
     do t = 1, 100; call v2%update(phi, nt); end do; res(i,j) = real(rank_threshold(v2%K_core, 1.0e-3_dp), dp)
  end do; end do
  open(30, file="phase_diagram_data_fortran.csv", status="replace")
  write(30, '(A)') "# Phase Diagram: rows=lambda, cols=noise"
  do i = 1, 5; do j = 1, 5
     if (j<5) then; write(30, '(F0.6,A)', advance='no') res(i,j), ','; else; write(30, '(F0.6)') res(i,j); end if
  end do; end do; close(30)

  ! 4. 最終パラメータ出力 (文言をPython版に統一)
  call estimate_sigma_f(vpcm%A, sigma_pix); f_num = (sigma_pix * 5.0_dp) / (1.5_dp * 0.55_dp)
  open(40, file="summary_fortran.txt", position="append")
  write(40, '(A)') ""
  write(40, '(A)') "最終光学設計パラメータ:"
  write(40, '(A,F0.3,A)') " - 推定 PSF sigma: ", sigma_pix, " px"
  write(40, '(A,F0.2)') " - 推奨レンズ F値: ", f_num
  write(40, '(A,F0.2)') " - 散逸係数 lambda: ", vpcm%lam
  write(40, '(A,F0.2)') " - ノイズレベル: ", vpcm%noise
  write(40, '(A)') ""
  write(40, '(A)') "全工程完了。すべてのファイルがルートディレクトリに保存されました。"
  close(40)

  print *, "全工程完了。summary_fortran.txt, evolution_metrics_fortran.csv などを出力しました。"

contains

  function diagonal(mat) result(d)
    real(dp), intent(in) :: mat(:,:)
    real(dp) :: d(size(mat,1))
    integer :: i
    do i = 1, size(mat,1); d(i) = mat(i,i); end do
  end function diagonal

  subroutine estimate_sigma_f(A, sigma)
    real(dp), intent(in) :: A(:,:); real(dp), intent(out) :: sigma
    real(dp) :: row(n), x(n), y(n), ata(2,2), atb(2), det, a_fit
    integer :: i, m, mid; mid = (n+1)/2; row = A(mid,:) / maxval(A(mid,:))
    m = 0; do i = 1, n
       if (row(i) > 1.0e-4_dp) then; m = m + 1; x(m) = real(i-mid,dp); y(m) = log(row(i)); end if
    end do
    ata = 0.0_dp; atb = 0.0_dp
    do i = 1, m
       ata(1,1) = ata(1,1) + x(i)**4; ata(1,2) = ata(1,2) + x(i)**2
       ata(2,1) = ata(2,1) + x(i)**2; ata(2,2) = ata(2,2) + 1.0_dp
       atb(1) = atb(1) + x(i)**2 * y(i); atb(2) = atb(2) + y(i)
    end do
    det = ata(1,1)*ata(2,2) - ata(1,2)*ata(2,1)
    a_fit = (atb(1)*ata(2,2) - atb(2)*ata(1,2)) / det
    sigma = sqrt(-1.0_dp/(2.0_dp*a_fit))
  end subroutine estimate_sigma_f

end program main
