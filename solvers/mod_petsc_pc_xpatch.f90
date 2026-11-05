!> The three logically rectangular patches of a lower single-null X-point grid
!! (core, SOL, private region), for the X-point GMG hierarchy (mod_petsc_pc_gmg).
!!
!! Element lattice (A, B): A along the flux surface (the element's t direction),
!! B across it (s); vertex slots v1 = (A, B), v2 = (A, B+1), v3 = (A+1, B+1),
!! v4 = (A+1, B). The lattice comes from a walk over element%neighbours in which
!! elements are adjacent through the vertex POSITIONS they share: at the X-point
!! four co-located nodes each serve only two of the eight elements around it.
!! Node lattice (I, J) = (B, A). The patches follow the grid's topology, not psi:
!! the core is the shallowest breadth-first layer from the axis that touches the
!! X-point, the private region the flood fill from its two X-point elements that
!! does not cross the leg separatrices, the rest the SOL.
module mod_petsc_pc_xpatch
  implicit none
  private

  public :: xpatch_t, xpatch_build, xpatch_report, xpatch_is_xgrid
  public :: XP_CORE, XP_SOL, XP_PRIV, XP_NPATCH

  integer, parameter :: XP_CORE = 1, XP_SOL = 2, XP_PRIV = 3, XP_NPATCH = 3
  character(len=8), parameter :: XP_NAME(XP_NPATCH) = ["core    ", "SOL     ", "private "]

  !> slot k (1..4) -> (dA, dB) of its vertex in the element lattice
  integer, parameter :: SLOT_A(4) = [0, 0, 1, 1]
  integer, parameter :: SLOT_B(4) = [0, 1, 1, 0]

  type :: xpatch_t
    character(len=12) :: name = ""
    integer :: na = 0, nb = 0                    !< element lattice: A = 0..na-1, B = 0..nb-1
    logical :: periodic = .false.                !< A periodic (the core)
    logical :: axis = .false.                    !< node row I = 0 is the magnetic axis
    integer, allocatable :: eat(:,:)             !< (0:na-1, 0:nb-1) -> element (1-based)
    integer :: ni = 0                            !< node rows I = 0..ni-1
    integer :: nxp = 0                           !< X-point positions in the node lattice
    integer :: xpI(4) = -1, xpJ(4) = -1
  end type xpatch_t

contains

  !> An X-point grid that the closed-grid GMG cannot take.
  logical function xpatch_is_xgrid()
    use phys_module,    only: xpoint, n_flux, n_tht
    use nodes_elements, only: node_list
    xpatch_is_xgrid = xpoint .and. node_list%n_nodes /= n_flux * n_tht
  end function xpatch_is_xgrid

  !> Build the three patches. ok = .false. (and msg) when the grid is not a
  !! lower single-null X-point grid of three logically rectangular patches.
  subroutine xpatch_build(xp, ok, msg)
    use nodes_elements, only: element_list
    use phys_module,    only: xcase, LOWER_XPOINT
    type(xpatch_t), intent(inout) :: xp(XP_NPATCH)
    logical, intent(out)          :: ok
    character(len=*), intent(out) :: msg

    integer :: ne, k, c, pax, pxp
    integer, allocatable :: pid(:), pcnt(:), cls(:)

    ok = .false.; msg = ""
    if (xcase /= LOWER_XPOINT) then
      msg = "only the lower single-null X-point grid is handled (xcase /= LOWER_XPOINT)"
      return
    endif
    ne = element_list%n_elements

    call position_ids(pid, pcnt)
    pax = maxloc(pcnt, 1)                                  ! the axis: most nodes at one position
    pxp = -1
    do k = 1, size(pcnt)
      if (pcnt(k) == 4 .and. k /= pax) then
        if (pxp > 0) then
          msg = "more than one position with four co-located nodes (double null?)"
          return
        endif
        pxp = k
      endif
    enddo
    if (pxp < 0) then
      msg = "no X-point (no position with four co-located nodes)"
      return
    endif
    allocate(cls(ne))
    call classify_topological(pid, pax, pxp, cls, ok, msg)
    if (.not. ok) return

    do c = 1, XP_NPATCH
      call build_patch(xp(c), c, cls, pid, pxp, ok, msg)
      if (.not. ok) return
    enddo
  end subroutine xpatch_build

  !> One patch: the element lattice by the position walk, its orientation and
  !! the X-point's node lattice positions.
  subroutine build_patch(P, c, cls, pid, pxp, ok, msg)
    use nodes_elements, only: node_list, element_list
    type(xpatch_t), intent(inout) :: P
    integer, intent(in)  :: c, cls(:), pid(:), pxp
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg

    integer :: ne, e, nb_, k, n, q, head, tail, nmem, dA, dB, a, b, amin, bmin, amax, bmax, i, j, nj
    integer, allocatable :: ca(:), cb(:), queue(:), npos(:,:), pnodes(:,:,:)
    logical, allocatable :: inP(:)
    logical :: is_core, wrap

    ok = .false.
    ne = element_list%n_elements
    P%name = trim(XP_NAME(c))
    is_core = (c == XP_CORE)
    P%periodic = is_core

    allocate(ca(ne), cb(ne), inP(ne), queue(ne))
    inP = .false.
    ca = huge(0); cb = huge(0)
    nmem = count(cls == c)
    if (nmem == 0) then
      msg = "empty patch "//trim(XP_NAME(c)); return
    endif
    e = findloc(cls, c, 1)
    ca(e) = 0; cb(e) = 0; inP(e) = .true.
    head = 1; tail = 1; queue(1) = e
    wrap = .false.
    do while (head <= tail)
      e = queue(head); head = head + 1
      do k = 1, 4
        nb_ = element_list%element(e)%neighbours(k)
        if (nb_ <= 0 .or. nb_ > ne) cycle
        if (cls(nb_) /= c) cycle
        if (.not. step(e, nb_, dA, dB)) cycle
        if (inP(nb_)) then
          if (cb(nb_) /= cb(e) + dB) then
            msg = "patch "//trim(XP_NAME(c))//" is not logically rectangular (B)"; return
          endif
          if (ca(nb_) /= ca(e) + dA) then
            if (.not. is_core) then
              msg = "patch "//trim(XP_NAME(c))//" is not logically rectangular (A)"; return
            endif
            wrap = .true.                       ! the core's periodic seam
          endif
        else
          ca(nb_) = ca(e) + dA; cb(nb_) = cb(e) + dB
          inP(nb_) = .true.
          tail = tail + 1; queue(tail) = nb_
        endif
      enddo
    enddo
    if (tail /= nmem) then
      msg = "patch "//trim(XP_NAME(c))//" is not connected"; return
    endif
    if (is_core .and. .not. wrap) then
      msg = "the core patch does not close on itself"; return
    endif

    !--- lattice extents
    amin = huge(0); bmin = huge(0); amax = -huge(0); bmax = -huge(0)
    do q = 1, tail
      e = queue(q)
      bmin = min(bmin, cb(e)); bmax = max(bmax, cb(e))
      if (.not. is_core) then
        amin = min(amin, ca(e)); amax = max(amax, ca(e))
      endif
    enddo
    if (is_core) then
      ! the walk numbers A without wrapping; the ring length is the number of
      ! elements of one B row
      P%na = count(cb(queue(1:tail)) == cb(queue(1)))
      amin = 0
      do q = 1, tail
        ca(queue(q)) = modulo(ca(queue(q)), P%na)
      enddo
      ! the seam on the X-point's column: the core's flux surfaces then run
      ! from the X-point to the X-point (J = 0 .. na), as the SOL's middle
      do q = 1, tail
        e = queue(q)
        do k = 1, 4
          if (pid(element_list%element(e)%vertex(k)) == pxp) exit
        enddo
        if (k <= 4) exit
      enddo
      if (q <= tail) then
        a = modulo(ca(e) + SLOT_A(k), P%na)
        do q = 1, tail
          ca(queue(q)) = modulo(ca(queue(q)) - a, P%na)
        enddo
      endif
      amax = P%na - 1
    else
      P%na = amax - amin + 1
    endif
    P%nb = bmax - bmin + 1
    if (P%na * P%nb /= tail) then
      msg = "patch "//trim(P%name)//" is not logically rectangular (element count)"; return
    endif
    allocate(P%eat(0:P%na - 1, 0:P%nb - 1))
    P%eat = 0
    do q = 1, tail
      e = queue(q)
      a = ca(e) - amin; b = cb(e) - bmin
      if (P%eat(a, b) /= 0) then
        msg = "patch "//trim(P%name)//": two elements at one lattice position"; return
      endif
      P%eat(a, b) = e
    enddo

    !--- orientation: the core's axis on B = 0
    P%axis = .false.
    if (is_core) then
      if (any_axis_row(P, 0)) then
        P%axis = .true.
      else if (any_axis_row(P, P%nb - 1)) then
        msg = "core patch: the axis lies on the last B row (orientation not handled)"; return
      else
        msg = "core patch: no axis row found"; return
      endif
    endif

    !--- node lattice: at most four nodes per position; the X-point's positions
    P%ni = P%nb + 1
    nj = P%na + 1
    if (P%periodic) nj = P%na
    allocate(npos(0:P%ni - 1, 0:nj - 1), pnodes(4, 0:P%ni - 1, 0:nj - 1))
    npos = 0; pnodes = 0
    P%nxp = 0
    do b = 0, P%nb - 1
      do a = 0, P%na - 1
        e = P%eat(a, b)
        do k = 1, 4
          n = element_list%element(e)%vertex(k)
          i = b + SLOT_B(k)
          j = a + SLOT_A(k)
          if (P%periodic) j = modulo(j, P%na)
          if (.not. any(pnodes(1:npos(i, j), i, j) == n)) then
            if (npos(i, j) == 4) then
              msg = "patch "//trim(P%name)//": more than four nodes at one lattice position"; return
            endif
            npos(i, j) = npos(i, j) + 1; pnodes(npos(i, j), i, j) = n
          endif
          if (pid(n) == pxp) then
            if (.not. any(P%xpI(1:P%nxp) == i .and. P%xpJ(1:P%nxp) == j)) then
              if (P%nxp < 4) then
                P%nxp = P%nxp + 1; P%xpI(P%nxp) = i; P%xpJ(P%nxp) = j
              endif
            endif
          endif
        enddo
      enddo
    enddo
    ok = .true.

  contains

    !> Lattice step from e to its edge neighbour nb2 through their shared
    !! vertex POSITIONS; .false. unless exactly one edge is shared.
    logical function step(e1, nb2, sA, sB)
      integer, intent(in)  :: e1, nb2
      integer, intent(out) :: sA, sB
      logical :: sh(4)
      integer :: k1, k2
      sh = .false.
      do k1 = 1, 4
        do k2 = 1, 4
          if (pid(element_list%element(e1)%vertex(k1)) == pid(element_list%element(nb2)%vertex(k2))) &
            sh(k1) = .true.
        enddo
      enddo
      step = (count(sh) == 2)
      sA = 0; sB = 0
      if (.not. step) return
      if (sh(1) .and. sh(2)) then
        sA = -1
      else if (sh(3) .and. sh(4)) then
        sA = 1
      else if (sh(2) .and. sh(3)) then
        sB = 1
      else if (sh(1) .and. sh(4)) then
        sB = -1
      else
        step = .false.
      endif
    end function step

    !> Does node row I = 0 (bb = 0) or I = nb (bb = nb-1) hold an axis node?
    logical function any_axis_row(PP, bb)
      type(xpatch_t), intent(in) :: PP
      integer, intent(in) :: bb
      integer :: aa, kk2, want
      any_axis_row = .false.
      want = 0
      if (bb /= 0) want = 1
      do aa = 0, PP%na - 1
        do kk2 = 1, 4
          if (SLOT_B(kk2) /= want) cycle
          if (node_list%node(element_list%element(PP%eat(aa, bb))%vertex(kk2))%axis_node) then
            any_axis_row = .true.; return
          endif
        enddo
      enddo
    end function any_axis_row

  end subroutine build_patch

  !> Position ids of the nodes: co-located nodes (the axis, the X-point)
  !! share one. pcnt(id) = nodes at that position.
  subroutine position_ids(pid, pcnt)
    use nodes_elements, only: node_list
    integer, allocatable, intent(out) :: pid(:), pcnt(:)
    integer(8), allocatable :: kr(:), kz(:)
    integer, allocatable :: ord(:)
    integer :: nn, n, np_
    real*8 :: scl

    nn = node_list%n_nodes
    allocate(kr(nn), kz(nn), ord(nn), pid(nn))
    scl = 0.d0
    do n = 1, nn
      scl = max(scl, maxval(abs(node_list%node(n)%x(1, 1, 1:2))))
    enddo
    scl = max(scl, 1.d-30) * 1.d-10
    do n = 1, nn
      kr(n) = nint(node_list%node(n)%x(1, 1, 1) / scl, 8)
      kz(n) = nint(node_list%node(n)%x(1, 1, 2) / scl, 8)
      ord(n) = n
    enddo
    call msort(ord, kr, kz)
    np_ = 0
    do n = 1, nn
      if (n == 1) then
        np_ = 1
      else if (kr(ord(n)) /= kr(ord(n - 1)) .or. kz(ord(n)) /= kz(ord(n - 1))) then
        np_ = np_ + 1
      endif
      pid(ord(n)) = np_
    enddo
    allocate(pcnt(np_))
    pcnt = 0
    do n = 1, nn
      pcnt(pid(n)) = pcnt(pid(n)) + 1
    enddo
  end subroutine position_ids

  !> cls(e) = XP_CORE / XP_SOL / XP_PRIV from the element adjacency alone (see
  !! the module header). pax / pxp: the axis and X-point position ids. ok =
  !! .false. (and msg) when the grid around the X-point is not the single-null
  !! pattern: eight elements in the cycle core, core, SOL, SOL, private,
  !! private, SOL, SOL.
  subroutine classify_topological(pid, pax, pxp, cls, ok, msg)
    use nodes_elements, only: element_list
    integer, intent(in)           :: pid(:), pax, pxp
    integer, intent(out)          :: cls(:)
    logical, intent(out)          :: ok
    character(len=*), intent(out) :: msg

    integer :: ne, e, k, nb_, q, head, tail, rstar, nx, i, leg, s, p, s2, p2, last, qpos, a, b
    integer :: xel(8), xadj(2, 8), ring(6), cores(2), nc, nsep
    integer, allocatable :: depth(:), queue(:), sep(:,:)
    logical :: found

    ok = .false.
    ne = element_list%n_elements
    allocate(depth(ne), queue(ne), sep(2, ne))

    !--- the core: breadth-first layers from the axis elements
    depth = -1
    tail = 0
    do e = 1, ne
      if (has_pos(e, pax)) then
        depth(e) = 0; tail = tail + 1; queue(tail) = e
      endif
    enddo
    head = 1
    do while (head <= tail)
      e = queue(head); head = head + 1
      do k = 1, 4
        nb_ = element_list%element(e)%neighbours(k)
        if (nb_ <= 0 .or. nb_ > ne) cycle
        if (depth(nb_) >= 0) cycle
        if (.not. edge_of(e, nb_, a, b)) cycle
        depth(nb_) = depth(e) + 1; tail = tail + 1; queue(tail) = nb_
      enddo
    enddo
    nx = 0
    do e = 1, ne
      if (.not. has_pos(e, pxp)) cycle
      nx = nx + 1
      if (nx > 8) exit
      xel(nx) = e
    enddo
    if (nx /= 8) then
      write(msg, '(A,I0,A)') "the X-point has ", nx, " elements around it, not 8"; return
    endif
    rstar = minval(depth(xel))
    if (rstar < 0) then
      msg = "the X-point is not connected to the axis"; return
    endif
    do e = 1, ne
      cls(e) = XP_SOL
      if (depth(e) >= 0 .and. depth(e) <= rstar) cls(e) = XP_CORE
    enddo

    !--- the cycle of the X-point's elements: core, core, s, s, p, p, s, s
    do i = 1, 8
      q = 0
      do k = 1, 8
        if (k == i) cycle
        if (.not. edge_of(xel(i), xel(k), a, b)) cycle
        if (a /= pxp .and. b /= pxp) cycle
        q = q + 1
        if (q <= 2) xadj(q, i) = k
      enddo
      if (q /= 2) then
        msg = "the elements around the X-point do not form a cycle"; return
      endif
    enddo
    nc = 0
    do i = 1, 8
      if (cls(xel(i)) /= XP_CORE) cycle
      nc = nc + 1
      if (nc <= 2) cores(nc) = i
    enddo
    found = (nc == 2)
    if (found) found = any(xadj(:, cores(1)) == cores(2))
    if (.not. found) then
      msg = "the X-point does not have two adjacent core elements"; return
    endif
    last = cores(1); i = cores(2)
    do q = 1, 6
      k = xadj(1, i)
      if (k == last) k = xadj(2, i)
      ring(q) = k; last = i; i = k
      if (cls(xel(k)) == XP_CORE) then
        msg = "more than two core elements at the X-point"; return
      endif
    enddo

    !--- the two leg separatrices, edge by edge from the X-point to the targets
    nsep = 0
    do leg = 1, 2
      if (leg == 1) then
        s = xel(ring(2)); p = xel(ring(3))
      else
        s = xel(ring(5)); p = xel(ring(4))
      endif
      last = pxp
      do
        if (.not. edge_of(s, p, a, b)) then
          msg = "broken leg separatrix"; return
        endif
        nsep = nsep + 1; sep(:, nsep) = [a, b]
        qpos = a
        if (a == last) qpos = b
        found = .false.
        do k = 1, 4
          s2 = element_list%element(s)%neighbours(k)
          if (.not. leg_step(s, s2, p, qpos)) cycle
          do i = 1, 4
            p2 = element_list%element(p)%neighbours(i)
            if (.not. leg_step(p, p2, s, qpos)) cycle
            if (.not. edge_of(s2, p2, a, b)) cycle
            if (a /= qpos .and. b /= qpos) cycle
            found = .true.; exit
          enddo
          if (found) exit
        enddo
        if (.not. found) exit                     ! qpos is on the target
        last = qpos; s = s2; p = p2
        if (nsep >= ne) then
          msg = "the leg separatrix does not end"; return
        endif
      enddo
    enddo

    !--- the private region: flood from its X-point elements, not across the legs
    tail = 0
    do q = 3, 4
      e = xel(ring(q)); cls(e) = XP_PRIV; tail = tail + 1; queue(tail) = e
    enddo
    head = 1
    do while (head <= tail)
      e = queue(head); head = head + 1
      do k = 1, 4
        nb_ = element_list%element(e)%neighbours(k)
        if (nb_ <= 0 .or. nb_ > ne) cycle
        if (cls(nb_) /= XP_SOL) cycle
        if (.not. edge_of(e, nb_, a, b)) cycle
        if (any(sep(1, 1:nsep) == a .and. sep(2, 1:nsep) == b)) cycle
        cls(nb_) = XP_PRIV; tail = tail + 1; queue(tail) = nb_
      enddo
    enddo
    ok = .true.

  contains

    logical function has_pos(e1, ps)
      integer, intent(in) :: e1, ps
      integer :: kk
      has_pos = .false.
      do kk = 1, 4
        if (pid(element_list%element(e1)%vertex(kk)) == ps) has_pos = .true.
      enddo
    end function has_pos

    !> e1 and e2 share an edge (two vertex positions); pa < pb its positions.
    logical function edge_of(e1, e2, pa, pb)
      integer, intent(in)  :: e1, e2
      integer, intent(out) :: pa, pb
      integer :: k1, k2, ns, sp(4), p1
      ns = 0
      do k1 = 1, 4
        p1 = pid(element_list%element(e1)%vertex(k1))
        if (any(sp(1:ns) == p1)) cycle
        do k2 = 1, 4
          if (pid(element_list%element(e2)%vertex(k2)) == p1) then
            ns = ns + 1; sp(ns) = p1; exit
          endif
        enddo
      enddo
      edge_of = (ns == 2)
      pa = 0; pb = 0
      if (edge_of) then
        pa = min(sp(1), sp(2)); pb = max(sp(1), sp(2))
      endif
    end function edge_of

    !> e2, an edge neighbour of e1 other than its partner o1 and not core,
    !! whose shared edge with e1 runs through position ps.
    logical function leg_step(e1, e2, o1, ps)
      integer, intent(in) :: e1, e2, o1, ps
      integer :: pa, pb
      leg_step = .false.
      if (e2 <= 0 .or. e2 > ne .or. e2 == o1) return
      if (cls(e2) == XP_CORE) return
      if (.not. edge_of(e1, e2, pa, pb)) return
      leg_step = (pa == ps .or. pb == ps)
    end function leg_step
  end subroutine classify_topological

  !> Stable merge sort of ord by (k1, k2).
  subroutine msort(ord, k1, k2)
    integer, intent(inout)  :: ord(:)
    integer(8), intent(in)  :: k1(:), k2(:)
    integer, allocatable :: tmp(:)
    integer :: n, w, lo, mid, hi, i, j, k
    n = size(ord)
    allocate(tmp(n))
    w = 1
    do while (w < n)
      lo = 1
      do while (lo <= n)
        mid = min(lo + w - 1, n); hi = min(lo + 2 * w - 1, n)
        i = lo; j = mid + 1; k = lo
        do while (i <= mid .and. j <= hi)
          if (less(ord(j), ord(i))) then
            tmp(k) = ord(j); j = j + 1
          else
            tmp(k) = ord(i); i = i + 1
          endif
          k = k + 1
        enddo
        do while (i <= mid)
          tmp(k) = ord(i); i = i + 1; k = k + 1
        enddo
        do while (j <= hi)
          tmp(k) = ord(j); j = j + 1; k = k + 1
        enddo
        lo = lo + 2 * w
      enddo
      ord = tmp
      w = 2 * w
    enddo
  contains
    logical function less(a, b)
      integer, intent(in) :: a, b
      less = (k1(a) < k1(b)) .or. (k1(a) == k1(b) .and. k2(a) < k2(b))
    end function less
  end subroutine msort

  !> One line per patch, rank 0.
  subroutine xpatch_report(xp, my_id)
    type(xpatch_t), intent(in) :: xp(XP_NPATCH)
    integer, intent(in) :: my_id
    integer :: c, q
    if (my_id /= 0) return
    do c = 1, XP_NPATCH
      write(*,'(A,A,A,I0,A,I0,A,L1,A,L1,A)', advance="no") "[Physics PC]   X-point patch ", trim(xp(c)%name), &
        ": ", xp(c)%na, "x", xp(c)%nb, " elements, periodic ", xp(c)%periodic, ", axis ", xp(c)%axis, &
        ", X-point at (I,J)"
      do q = 1, xp(c)%nxp
        write(*,'(A,I0,A,I0,A)', advance="no") " (", xp(c)%xpI(q), ",", xp(c)%xpJ(q), ")"
      enddo
      write(*,*)
    enddo
  end subroutine xpatch_report

end module mod_petsc_pc_xpatch
