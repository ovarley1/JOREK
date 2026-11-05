/* Galerkin coarse operator C = P^T A P for the GMG's packed harmonic layout.
 *
 * A (MPIAIJ) keeps the bs slots of a fine block, one (field, scalar DOF), in bs
 * consecutive rows; P = Ps x I_bs, Ps (MPIAIJ) one row per fine block, one column per
 * coarse block. C (MPIAIJ, dense bs x bs blocks) takes its block pattern and layout
 * from PETSc's PtAP of A's block graph. Per coarse block row I touched here (one
 * thread each) R = sum_i Ps(i,I) A(i,:), a dense strip, and C(I,:) = R Ps; owned rows
 * go into C's CSR arrays, the others to their owner as dense blocks, added there in
 * rank order at offsets fixed by the symbolic phase. reuse = 1 refills a C made here
 * (same A pattern). *ok = 0 (nothing done, on all ranks) if A, Ps or the layout do not
 * fit. Collective. Requires PETSc >= 3.26.
 */
#if defined(USE_PETSC)
#include <stdint.h>
#include <string.h>
#include <petscmat.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

#define GBLK_KEY "jorek_galerkin_blk"

static PetscLogEvent gblk_ev_sym = -1, gblk_ev_num = -1, gblk_ev_asm = -1;

typedef struct {
  PetscInt  bs, nlb, ngb, ncb, cstart, nlc; /* fine local / ghost blocks, coarse local / first blocks, coarse blocks Ps reaches */
  PetscInt *dcb, *dcn;      /* A's diagonal part: local column -> local fine block, slot */
  PetscInt *ocb, *ocn;      /* A's off-diagonal part: column -> fine block id (nlb + ghost), slot */
  PetscBool *ofb;           /* off-diagonal column starting a whole slot block (bs columns, slots in order) */
  PetscInt *pp, *pl;        /* Ps rows of the fine blocks (local, then ghosts): CSR, local coarse block */
  PetscScalar *pv;
  PetscInt  nt;             /* coarse block rows touched here */
  PetscInt *tI;             /* their global ids */
  PetscInt *cp, *ci;        /* contributors: local fine blocks i ... */
  PetscScalar *cv;          /* ... with Ps(i, I) */
  PetscInt *jp, *jl;        /* fine block columns of R per touched row */
  PetscInt *Jp, *Jl;        /* local coarse block columns of the partial row, ascending */
  PetscInt *Jo;             /* owned rows: offset of block Jl[k] in its row, diagonal part >= 0, off-diagonal part -1 - offset */
  PetscInt *rb;             /* remote rows: offset into rv (-1 for owned), grouped by owner */
  PetscScalar *rv;
  /* the exchange: sends of rv to nsd owners (sdr ranks, sdo offsets into
     rv), receives from nrc ranks (rcr, rco offsets into rcv); received row k:
     its local coarse block row inI, block list inp / inJo (as Jo), values at
     rcv + ino; the rows of one target grouped (igp / igk, k ascending) */
  PetscMPIInt nsd, nrc, *sdr, *rcr, tag;
  PetscInt *sdo, *rco, nin, *inI, *inp, *inJo, *ino, nig, *igp, *igk;
  PetscScalar *rcv;
  PetscInt  mxj, mxJ;       /* largest R / out strip, in blocks */
  /* what a reuse is checked against */
  PetscInt  nA, nAo, nzd, nzo;
  uint64_t  hash;
  PetscObjectId    aid, pid, cid;
  PetscObjectState anz, cnz;
} gblk_ctx;

static PetscErrorCode ctx_destroy(void *pctx)
{
  gblk_ctx *c = *(gblk_ctx **)pctx;

  PetscFunctionBegin;
  PetscCall(PetscFree5(c->dcb, c->dcn, c->ocb, c->ocn, c->ofb));
  PetscCall(PetscFree3(c->pp, c->pl, c->pv));
  PetscCall(PetscFree4(c->tI, c->cp, c->ci, c->cv));
  PetscCall(PetscFree2(c->jp, c->Jp));
  PetscCall(PetscFree3(c->jl, c->Jl, c->Jo));
  PetscCall(PetscFree(c->rb));
  PetscCall(PetscFree(c->rv));
  PetscCall(PetscFree2(c->sdr, c->sdo));
  PetscCall(PetscFree2(c->rcr, c->rco));
  PetscCall(PetscFree4(c->inI, c->inp, c->inJo, c->ino));
  PetscCall(PetscFree2(c->igp, c->igk));
  PetscCall(PetscFree(c->rcv));
  PetscCall(PetscFree(c));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Row pointers and column indices of a SeqAIJ matrix, in host memory. */
static PetscErrorCode csr(Mat M, const PetscInt **i, const PetscInt **j)
{
  PetscMemType mt;

  PetscFunctionBegin;
  PetscCall(MatSeqAIJGetCSRAndMemType(M, i, j, NULL, &mt));
  PetscCheck(PetscMemTypeHost(mt), PETSC_COMM_SELF, PETSC_ERR_SUP, "galerkin_blk: needs host CSR arrays");
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Hash of A's local pattern (both parts, the ghost column map); *sorted:
   every row's columns strictly ascending. */
static PetscErrorCode pattern_scan(Mat A, uint64_t *hash, PetscBool *sorted)
{
  Mat             Ad, Ao;
  const PetscInt *cmap, *ia, *ja, *oia, *oja;
  PetscInt        n, nAo, i, k;
  uint64_t        h = 14695981039346656037ULL;
  PetscBool       s = PETSC_TRUE;

  PetscFunctionBegin;
  PetscCall(MatGetLocalSize(A, &n, NULL));
  PetscCall(MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, &cmap));
  PetscCall(MatGetSize(Ao, NULL, &nAo));
  PetscCall(csr(Ad, &ia, &ja));
  PetscCall(csr(Ao, &oia, &oja));
#define GBLK_MIX(v) (h = (h ^ (uint64_t)(v)) * 1099511628211ULL)
  GBLK_MIX(n); GBLK_MIX(nAo);
  for (i = 0; i < n; i++) {
    GBLK_MIX(ia[i + 1] - ia[i]); GBLK_MIX(oia[i + 1] - oia[i]);
    for (k = ia[i]; k < ia[i + 1]; k++) { GBLK_MIX(ja[k]); if (k > ia[i] && ja[k] <= ja[k - 1]) s = PETSC_FALSE; }
    for (k = oia[i]; k < oia[i + 1]; k++) { GBLK_MIX(oja[k]); if (k > oia[i] && oja[k] <= oja[k - 1]) s = PETSC_FALSE; }
  }
  for (k = 0; k < nAo; k++) GBLK_MIX(cmap[k]);
#undef GBLK_MIX
  *hash = h; *sorted = s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Rank owning global coarse block row Ib (rng: C's row ownership ranges). */
static PetscMPIInt gblk_owner(const PetscInt *rng, PetscMPIInt np, PetscInt Ib, PetscInt bs)
{
  PetscMPIInt lo = 0, hi = np - 1;
  while (lo < hi) {
    PetscMPIInt mid = (lo + hi + 1) / 2;
    if (rng[mid] <= Ib * bs) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/* Remote rows (coarse block rows touched here, owned elsewhere): rv offsets
 * grouped by owner (rb), and the exchange with the owners. Every rank sends
 * each owner its rows' (Ib, nJ, J list); the owner turns each received block
 * J into its offset in its CSR row (as Jo), so the numeric phase only moves
 * values. lcb: global id of each local coarse block. Jpos: work array of Ncb
 * entries, all PETSC_MIN_INT (restored). */
static PetscErrorCode gblk_exchange_setup(Mat C, gblk_ctx *c, const PetscInt *lcb, const PetscInt *cia, const PetscInt *cja,
                                          const PetscInt *coia, const PetscInt *coja, const PetscInt *cgarr,
                                          PetscInt *Jpos, PetscInt *nrv)
{
  MPI_Comm        comm = PetscObjectComm((PetscObject)C);
  const PetscInt  bs = c->bs;
  const PetscInt *rng;
  PetscMPIInt     np, *own, *sc, *rc, *sd, *rd, p;
  PetscInt        t, k, q, rem = 0, *sh, *rh, nsh, nrh, *hcnt, *tcnt, *ord, nrem = 0;

  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(comm, &np));
  PetscCall(PetscObjectGetNewTag((PetscObject)C, &c->tag));
  PetscCall(MatGetOwnershipRanges(C, &rng));
  PetscCall(PetscMalloc1(c->nt + 1, &own));
  PetscCall(PetscCalloc4(np, &sc, np, &rc, np + 1, &sd, np + 1, &rd));
  PetscCall(PetscCalloc2(np, &hcnt, np + 1, &tcnt));
  /* owners, header ints and rows per owner */
  for (t = 0; t < c->nt; t++) {
    own[t] = -1;
    if (c->tI[t] >= c->cstart && c->tI[t] < c->cstart + c->ncb) continue;
    own[t] = gblk_owner(rng, np, c->tI[t], bs);
    hcnt[own[t]] += 2 + (c->Jp[t + 1] - c->Jp[t]);
    tcnt[own[t] + 1]++;
    nrem++;
  }
  for (p = 0; p < np; p++) PetscCall(PetscMPIIntCast(hcnt[p], &sc[p]));
  PetscCallMPI(MPI_Alltoall(sc, 1, MPI_INT, rc, 1, MPI_INT, comm));
  for (p = 0; p < np; p++) {
    PetscCall(PetscMPIIntCast((PetscInt)sd[p] + sc[p], &sd[p + 1]));
    PetscCall(PetscMPIIntCast((PetscInt)rd[p] + rc[p], &rd[p + 1]));
  }
  nsh = sd[np]; nrh = rd[np];
  PetscCall(PetscMalloc2(nsh + 1, &sh, nrh + 1, &rh));
  /* remote rows ordered by owner (counting sort, ascending t within) */
  for (p = 0; p < np; p++) tcnt[p + 1] += tcnt[p];
  PetscCall(PetscMalloc1(nrem + 1, &ord));
  for (t = 0; t < c->nt; t++)
    if (own[t] >= 0) ord[tcnt[own[t]]++] = t;
  /* rv offsets and headers in that order, one send per owner */
  c->nsd = 0;
  for (p = 0; p < np; p++) if (sc[p]) c->nsd++;
  PetscCall(PetscMalloc2(c->nsd, &c->sdr, c->nsd + 1, &c->sdo));
  c->nsd = 0; q = 0;
  for (k = 0; k < nrem; k++) {
    PetscInt nJ, J;
    t = ord[k];
    if (k == 0 || own[t] != own[ord[k - 1]]) { c->sdr[c->nsd] = own[t]; c->sdo[c->nsd] = rem; c->nsd++; }
    nJ = c->Jp[t + 1] - c->Jp[t];
    c->rb[t] = rem; rem += nJ * bs * bs;
    sh[q++] = c->tI[t]; sh[q++] = nJ;
    for (J = c->Jp[t]; J < c->Jp[t + 1]; J++) sh[q++] = lcb[c->Jl[J]];
  }
  c->sdo[c->nsd] = rem;
  PetscCallMPI(MPI_Alltoallv(sh, sc, sd, MPIU_INT, rh, rc, rd, MPIU_INT, comm));
  /* received rows: local block row, block offsets, value offset */
  c->nrc = 0;
  for (p = 0; p < np; p++) if (rc[p]) c->nrc++;
  PetscCall(PetscMalloc2(c->nrc, &c->rcr, c->nrc + 1, &c->rco));
  c->nin = 0;
  for (k = 0; k < nrh; k += 2 + rh[k + 1]) c->nin++;
  PetscCall(PetscMalloc4(c->nin, &c->inI, c->nin + 1, &c->inp, nrh - 2 * c->nin + 1, &c->inJo, c->nin + 1, &c->ino));
  c->nin = 0; c->inp[0] = 0; c->nrc = 0;
  {
    PetscInt nv = 0;
    for (p = 0; p < np; p++) {
      if (!rc[p]) continue;
      c->rcr[c->nrc] = p; c->rco[c->nrc] = nv; c->nrc++;
      for (k = rd[p]; k < rd[p + 1]; k += 2 + rh[k + 1]) {
        PetscInt Ib = rh[k], nJ = rh[k + 1], r = (Ib - c->cstart) * bs, kk;
        PetscCheck(Ib >= c->cstart && Ib < c->cstart + c->ncb, PETSC_COMM_SELF, PETSC_ERR_PLIB, "galerkin_blk: received a row not owned here");
        for (kk = cia[r]; kk < cia[r + 1]; kk += bs) Jpos[c->cstart + cja[kk] / bs] = kk - cia[r];
        for (kk = coia[r]; kk < coia[r + 1]; kk += bs) Jpos[cgarr[coja[kk]] / bs] = -1 - (kk - coia[r]);
        c->inI[c->nin] = Ib - c->cstart;
        for (q = 0; q < nJ; q++) {
          PetscCheck(Jpos[rh[k + 2 + q]] != PETSC_MIN_INT, PETSC_COMM_SELF, PETSC_ERR_PLIB, "galerkin_blk: received block not in C's row");
          c->inJo[c->inp[c->nin] + q] = Jpos[rh[k + 2 + q]];
        }
        for (kk = cia[r]; kk < cia[r + 1]; kk += bs) Jpos[c->cstart + cja[kk] / bs] = PETSC_MIN_INT;
        for (kk = coia[r]; kk < coia[r + 1]; kk += bs) Jpos[cgarr[coja[kk]] / bs] = PETSC_MIN_INT;
        c->inp[c->nin + 1] = c->inp[c->nin] + nJ;
        c->ino[c->nin] = nv;
        c->nin++;
        nv += nJ * bs * bs;
      }
    }
    c->rco[c->nrc] = nv;
    PetscCall(PetscMalloc1(nv + 1, &c->rcv));
  }
  /* received rows grouped by target row (counting sort, ascending k within) */
  {
    PetscInt *pos, r;
    PetscCall(PetscCalloc1(c->ncb + 1, &pos));
    for (k = 0; k < c->nin; k++) pos[c->inI[k] + 1]++;
    c->nig = 0;
    for (r = 0; r < c->ncb; r++) if (pos[r + 1]) c->nig++;
    for (r = 0; r < c->ncb; r++) pos[r + 1] += pos[r];
    PetscCall(PetscMalloc2(c->nig + 1, &c->igp, c->nin + 1, &c->igk));
    c->nig = 0; c->igp[0] = 0;
    for (r = 0; r < c->ncb; r++) if (pos[r + 1] > pos[r]) c->igp[++c->nig] = pos[r + 1];
    for (k = 0; k < c->nin; k++) c->igk[pos[c->inI[k]]++] = k;
    PetscCall(PetscFree(pos));
  }
  PetscCall(PetscFree2(sh, rh));
  PetscCall(PetscFree(ord));
  PetscCall(PetscFree2(hcnt, tcnt));
  PetscCall(PetscFree(own));
  PetscCall(PetscFree4(sc, rc, sd, rd));
  *nrv = rem;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Symbolic phase: the context and C (dense slot blocks, all zero). */
static PetscErrorCode gblk_symbolic(Mat A, Mat Ps, PetscInt bs, Mat *Cout, gblk_ctx **cout)
{
  MPI_Comm        comm = PetscObjectComm((PetscObject)A);
  gblk_ctx       *c;
  Mat             Ad, Ao, Ablk, Acb, C, *Psub;
  const PetscInt *cmap, *ia, *ja, *oia, *oja;
  PetscInt        nA, rstart, nAo, Ncb, n, i, k, m, t, rbs, *gblk, *lcb, *mark, *bgp, *bgj;
  IS              isr, isc;

  PetscFunctionBegin;
  PetscCall(PetscNew(&c));
  c->bs = bs;
  PetscCall(MatGetLocalSize(A, &nA, NULL));
  PetscCall(MatGetOwnershipRange(A, &rstart, NULL));
  c->nlb = nA / bs;
  rbs    = rstart / bs;
  PetscCall(MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, &cmap));
  PetscCall(MatGetSize(Ao, NULL, &nAo));
  PetscCall(MatGetOwnershipRangeColumn(Ps, &c->cstart, NULL));
  PetscCall(MatGetLocalSize(Ps, NULL, &c->ncb));
  PetscCall(MatGetSize(Ps, NULL, &Ncb));
  PetscCall(csr(Ad, &ia, &ja));
  PetscCall(csr(Ao, &oia, &oja));
  c->nA = nA; c->nAo = nAo; c->nzd = ia[nA]; c->nzo = oia[nA];

  /* fine block ids of A's columns; ghosts: the off-diagonal part's blocks,
     ascending (the column map need not be sorted) */
  PetscCall(PetscMalloc5(nA, &c->dcb, nA, &c->dcn, nAo, &c->ocb, nAo, &c->ocn, nAo, &c->ofb));
  for (k = 0; k < nA; k++) { c->dcb[k] = k / bs; c->dcn[k] = k % bs; }
  PetscCall(PetscMalloc1(nAo + 1, &gblk));
  for (k = 0; k < nAo; k++) gblk[k] = cmap[k] / bs;
  c->ngb = nAo;
  PetscCall(PetscSortRemoveDupsInt(&c->ngb, gblk));
  {
    PetscInt run = 0, g; /* consecutive global columns from k on */
    for (k = nAo - 1; k >= 0; k--) {
      run = (k + 1 < nAo && cmap[k + 1] == cmap[k] + 1) ? run + 1 : 1;
      PetscCall(PetscFindInt(cmap[k] / bs, c->ngb, gblk, &g));
      c->ocb[k] = c->nlb + g;
      c->ocn[k] = cmap[k] % bs;
      c->ofb[k] = (PetscBool)(c->ocn[k] == 0 && run >= bs);
    }
  }

  /* Ps rows of the local and ghost fine blocks (Ps is fixed: once); their
     coarse blocks by local id (lcb: the global ids, ascending) */
  {
    PetscInt *rows, nr = c->nlb + c->ngb, *pos, a = 0, b = 0, q = 0;
    PetscCall(PetscMalloc2(nr, &rows, nr, &pos));
    /* merge the local range into the sorted ghosts: sorted row IS */
    while (a < c->nlb || b < c->ngb) {
      if (b >= c->ngb || (a < c->nlb && rbs + a < gblk[b])) { pos[a] = q; rows[q++] = rbs + a; a++; }
      else { pos[c->nlb + b] = q; rows[q++] = gblk[b]; b++; }
    }
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, nr, rows, PETSC_USE_POINTER, &isr));
    PetscCall(ISCreateStride(PETSC_COMM_SELF, Ncb, 0, 1, &isc));
    PetscCall(MatCreateSubMatrices(Ps, 1, &isr, &isc, MAT_INITIAL_MATRIX, &Psub));
    PetscCall(ISDestroy(&isr));
    PetscCall(ISDestroy(&isc));
    {
      const PetscInt    *si, *sj;
      const PetscScalar *sv;
      PetscInt           f, nz = 0;
      PetscCall(csr(Psub[0], &si, &sj));
      PetscCall(MatSeqAIJGetArrayRead(Psub[0], &sv));
      PetscCall(PetscMalloc3(nr + 1, &c->pp, si[nr], &c->pl, si[nr], &c->pv));
      c->pp[0] = 0;
      for (f = 0; f < nr; f++) {
        for (k = si[pos[f]]; k < si[pos[f] + 1]; k++)
          if (sv[k] != 0.0) { c->pl[nz] = sj[k]; c->pv[nz] = sv[k]; nz++; }
        c->pp[f + 1] = nz;
      }
      PetscCall(MatSeqAIJRestoreArrayRead(Psub[0], &sv));
      PetscCall(PetscMalloc1(nz + 1, &lcb));
      PetscCall(PetscArraycpy(lcb, c->pl, nz));
      c->nlc = nz;
      PetscCall(PetscSortRemoveDupsInt(&c->nlc, lcb));
      for (k = 0; k < nz; k++) PetscCall(PetscFindInt(c->pl[k], c->nlc, lcb, &c->pl[k]));
    }
    PetscCall(MatDestroySubMatrices(1, &Psub));
    PetscCall(PetscFree2(rows, pos));
  }

  /* the block graph of A (value 1) and its PtAP: C's block pattern and
     layout; the graph is kept (bgp / bgj: local fine block -> its fine
     block columns, local ids as dcb / ocb) for the touched rows' lists */
  {
    PetscInt    *dnz, *onz, *cols, *gcols, nc, gi;
    PetscScalar *one;
    PetscCall(PetscMalloc6(c->nlb, &dnz, c->nlb, &onz, c->nlb + c->ngb, &mark, c->nlb + c->ngb, &cols, c->nlb + c->ngb, &gcols, c->nlb + c->ngb, &one));
    for (k = 0; k < c->nlb + c->ngb; k++) { mark[k] = -1; one[k] = 1.0; }
    for (i = 0; i < c->nlb; i++) {
      dnz[i] = onz[i] = 0;
      for (m = 0; m < bs; m++) {
        PetscInt r = i * bs + m;
        for (k = ia[r]; k < ia[r + 1]; k++) if (mark[c->dcb[ja[k]]] != i) { mark[c->dcb[ja[k]]] = i; dnz[i]++; }
        for (k = oia[r]; k < oia[r + 1]; k++) if (mark[c->ocb[oja[k]]] != i) { mark[c->ocb[oja[k]]] = i; onz[i]++; }
      }
    }
    PetscCall(PetscMalloc1(c->nlb + 1, &bgp));
    bgp[0] = 0;
    for (i = 0; i < c->nlb; i++) bgp[i + 1] = bgp[i] + dnz[i] + onz[i];
    PetscCall(PetscMalloc1(bgp[c->nlb] + 1, &bgj));
    PetscCall(MatCreate(comm, &Ablk));
    PetscCall(MatSetSizes(Ablk, c->nlb, c->nlb, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetType(Ablk, MATMPIAIJ));
    PetscCall(MatMPIAIJSetPreallocation(Ablk, 0, dnz, 0, onz));
    for (k = 0; k < c->nlb + c->ngb; k++) mark[k] = -1;
    for (i = 0; i < c->nlb; i++) {
      nc = 0;
      for (m = 0; m < bs; m++) {
        PetscInt r = i * bs + m;
        for (k = ia[r]; k < ia[r + 1]; k++) if (mark[c->dcb[ja[k]]] != i) { mark[c->dcb[ja[k]]] = i; cols[nc++] = c->dcb[ja[k]]; }
        for (k = oia[r]; k < oia[r + 1]; k++) if (mark[c->ocb[oja[k]]] != i) { mark[c->ocb[oja[k]]] = i; cols[nc++] = c->ocb[oja[k]]; }
      }
      PetscCall(PetscArraycpy(bgj + bgp[i], cols, nc));
      for (k = 0; k < nc; k++) gcols[k] = cols[k] < c->nlb ? rbs + cols[k] : gblk[cols[k] - c->nlb];
      gi = rbs + i;
      PetscCall(MatSetValues(Ablk, 1, &gi, nc, gcols, one, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(Ablk, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Ablk, MAT_FINAL_ASSEMBLY));
    PetscCall(MatPtAP(Ablk, Ps, MAT_INITIAL_MATRIX, 2.0, &Acb));
    PetscCall(MatDestroy(&Ablk));
    PetscCall(PetscFree6(dnz, onz, mark, cols, gcols, one));
  }

  /* C: Acb's block pattern with dense bs x bs blocks */
  {
    PetscInt *ci_, *cj_, nzb = 0, Ib, q;
    for (Ib = 0; Ib < c->ncb; Ib++) {
      PetscCall(MatGetRow(Acb, c->cstart + Ib, &n, NULL, NULL));
      nzb += n;
      PetscCall(MatRestoreRow(Acb, c->cstart + Ib, &n, NULL, NULL));
    }
    PetscCall(PetscMalloc2(c->ncb * bs + 1, &ci_, nzb * bs * bs, &cj_));
    ci_[0] = 0; q = 0;
    for (Ib = 0; Ib < c->ncb; Ib++) {
      const PetscInt *bj;
      PetscCall(MatGetRow(Acb, c->cstart + Ib, &n, &bj, NULL));
      for (m = 0; m < bs; m++) {
        for (k = 0; k < n; k++) for (i = 0; i < bs; i++) cj_[q++] = bj[k] * bs + i;
        ci_[Ib * bs + m + 1] = q;
      }
      PetscCall(MatRestoreRow(Acb, c->cstart + Ib, &n, &bj, NULL));
    }
    PetscCall(MatCreate(comm, &C));
    PetscCall(MatSetSizes(C, c->ncb * bs, c->ncb * bs, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetBlockSizes(C, bs, bs));
    PetscCall(MatSetType(C, MATMPIAIJ));
    PetscCall(MatMPIAIJSetPreallocationCSR(C, ci_, cj_, NULL));
    PetscCall(MatSetOption(C, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
    PetscCall(PetscFree2(ci_, cj_));
    PetscCall(MatDestroy(&Acb));
  }

  /* touched coarse block rows: Ps^T of the local fine blocks (local ids
     ascend with the global ones) */
  {
    PetscInt *cnt, *Imark, *tpos;
    PetscCall(PetscCalloc1(c->nlc + 1, &Imark));
    for (i = 0; i < c->nlb; i++)
      for (k = c->pp[i]; k < c->pp[i + 1]; k++) Imark[c->pl[k]]++;
    c->nt = 0;
    for (k = 0; k < c->nlc; k++) if (Imark[k]) c->nt++;
    PetscCall(PetscMalloc4(c->nt, &c->tI, c->nt + 1, &c->cp, c->pp[c->nlb], &c->ci, c->pp[c->nlb], &c->cv));
    PetscCall(PetscMalloc2(c->nlc + 1, &tpos, c->nt + 1, &cnt));
    c->cp[0] = 0; t = 0;
    for (k = 0; k < c->nlc; k++) {
      tpos[k] = -1;
      if (Imark[k]) { c->tI[t] = lcb[k]; tpos[k] = t; c->cp[t + 1] = c->cp[t] + Imark[k]; cnt[t] = 0; t++; }
    }
    for (i = 0; i < c->nlb; i++)
      for (k = c->pp[i]; k < c->pp[i + 1]; k++) {
        t = tpos[c->pl[k]];
        c->ci[c->cp[t] + cnt[t]] = i; c->cv[c->cp[t] + cnt[t]] = c->pv[k]; cnt[t]++;
      }
    PetscCall(PetscFree(Imark));
    PetscCall(PetscFree2(tpos, cnt));
  }

  /* per touched row: fine block columns of R, coarse block columns of the
     partial row, and where they go */
  {
    PetscInt *jmark, *Jmark, nj = 0, nJ = 0, cap_j, cap_J, *jbuf, *Jbuf;
    PetscCall(PetscMalloc2(c->nlb + c->ngb + 1, &jmark, c->nlc + 1, &Jmark));
    for (k = 0; k < c->nlb + c->ngb; k++) jmark[k] = -1;
    for (k = 0; k < c->nlc; k++) Jmark[k] = -1;
    cap_j = cap_J = 64 * (c->nt + 1);
    PetscCall(PetscMalloc1(cap_j, &jbuf));  /* grown one by one: separate allocations */
    PetscCall(PetscMalloc1(cap_J, &Jbuf));
    PetscCall(PetscMalloc2(c->nt + 1, &c->jp, c->nt + 1, &c->Jp));
    c->jp[0] = c->Jp[0] = 0; c->mxj = c->mxJ = 0;
#define GBLK_PUSH(buf, cap, len, v) do { \
      if (len == cap) { PetscInt *x_; cap *= 2; PetscCall(PetscMalloc1(cap, &x_)); PetscCall(PetscArraycpy(x_, buf, len)); PetscCall(PetscFree(buf)); buf = x_; } \
      buf[len++] = (v); } while (0)
    for (t = 0; t < c->nt; t++) {
      PetscInt j0 = nj, J0 = nJ, q, kk;
      for (q = c->cp[t]; q < c->cp[t + 1]; q++) {
        i = c->ci[q];
        for (kk = bgp[i]; kk < bgp[i + 1]; kk++) {
          PetscInt b = bgj[kk];
          if (jmark[b] != t) { jmark[b] = t; GBLK_PUSH(jbuf, cap_j, nj, b); }
        }
      }
      for (q = j0; q < nj; q++)
        for (kk = c->pp[jbuf[q]]; kk < c->pp[jbuf[q] + 1]; kk++) {
          PetscInt J = c->pl[kk];
          if (Jmark[J] != t) { Jmark[J] = t; GBLK_PUSH(Jbuf, cap_J, nJ, J); }
        }
      PetscCall(PetscSortInt(nJ - J0, Jbuf + J0));
      c->jp[t + 1] = nj; c->Jp[t + 1] = nJ;
      c->mxj = PetscMax(c->mxj, nj - j0);
      c->mxJ = PetscMax(c->mxJ, nJ - J0);
    }
#undef GBLK_PUSH
    PetscCall(PetscMalloc3(nj, &c->jl, nJ, &c->Jl, nJ, &c->Jo));
    PetscCall(PetscArraycpy(c->jl, jbuf, nj));
    PetscCall(PetscArraycpy(c->Jl, Jbuf, nJ));
    PetscCall(PetscFree(jbuf));
    PetscCall(PetscFree(Jbuf));
    PetscCall(PetscFree2(jmark, Jmark));
  }

  /* owned rows: block offsets in C's diagonal / off-diagonal row (C's
     off-diagonal columns come in whole slot blocks, slots in order) */
  {
    Mat             Cd, Co;
    const PetscInt *cgarr, *cia, *cja, *coia, *coja;
    PetscInt        rem = 0, nCo, *Jpos;
    PetscCall(MatMPIAIJGetSeqAIJ(C, &Cd, &Co, &cgarr));
    PetscCall(MatGetSize(Co, NULL, &nCo));
    for (k = 0; k < nCo; k++)
      PetscCheck(cgarr[k] % bs == k % bs && cgarr[k] - cgarr[k - k % bs] == k % bs, PETSC_COMM_SELF, PETSC_ERR_PLIB, "galerkin_blk: C's ghost columns not in slot blocks");
    PetscCall(csr(Cd, &cia, &cja));
    PetscCall(csr(Co, &coia, &coja));
    PetscCall(PetscMalloc1(c->nt + 1, &c->rb));
    PetscCall(PetscMalloc1(Ncb, &Jpos));
    for (k = 0; k < Ncb; k++) Jpos[k] = PETSC_MIN_INT;
    for (t = 0; t < c->nt; t++) {
      PetscInt Ib = c->tI[t], q;
      if (Ib < c->cstart || Ib >= c->cstart + c->ncb) continue;   /* remote: gblk_exchange_setup */
      c->rb[t] = -1;
      {
        PetscInt r = (Ib - c->cstart) * bs;
        for (k = cia[r]; k < cia[r + 1]; k += bs) Jpos[c->cstart + cja[k] / bs] = k - cia[r];
        for (k = coia[r]; k < coia[r + 1]; k += bs) Jpos[cgarr[coja[k]] / bs] = -1 - (k - coia[r]);
        for (q = c->Jp[t]; q < c->Jp[t + 1]; q++) {
          PetscCheck(Jpos[lcb[c->Jl[q]]] != PETSC_MIN_INT, PETSC_COMM_SELF, PETSC_ERR_PLIB, "galerkin_blk: block %" PetscInt_FMT " not in C's row", lcb[c->Jl[q]]);
          c->Jo[q] = Jpos[lcb[c->Jl[q]]];
        }
        for (k = cia[r]; k < cia[r + 1]; k += bs) Jpos[c->cstart + cja[k] / bs] = PETSC_MIN_INT;
        for (k = coia[r]; k < coia[r + 1]; k += bs) Jpos[cgarr[coja[k]] / bs] = PETSC_MIN_INT;
      }
    }
    PetscCall(gblk_exchange_setup(C, c, lcb, cia, cja, coia, coja, cgarr, Jpos, &rem));
    PetscCall(PetscMalloc1(rem + 1, &c->rv));
    PetscCall(PetscFree(Jpos));
  }

  PetscCall(PetscFree(gblk));
  PetscCall(PetscFree(lcb));
  PetscCall(PetscFree(bgp));
  PetscCall(PetscFree(bgj));
  *Cout = C;
  *cout = c;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Numeric phase: C's values from A's. */
static PetscErrorCode gblk_numeric(Mat A, Mat C, gblk_ctx *c)
{
  Mat                Ad, Ao, Cd, Co;
  const PetscInt    *ia, *ja, *oia, *oja, *cia, *coia;
  const PetscScalar *ad, *ao;
  PetscScalar       *cd, *co;
  const PetscInt     bs = c->bs, ncr = c->ncb * c->bs;
  int                fail = 0;

  PetscFunctionBegin;
  PetscCall(MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, NULL));
  PetscCall(MatMPIAIJGetSeqAIJ(C, &Cd, &Co, NULL));
  PetscCall(csr(Ad, &ia, &ja));
  PetscCall(csr(Ao, &oia, &oja));
  PetscCall(csr(Cd, &cia, NULL));
  PetscCall(csr(Co, &coia, NULL));
  PetscCall(MatSeqAIJGetArrayRead(Ad, &ad));
  PetscCall(MatSeqAIJGetArrayRead(Ao, &ao));
  PetscCall(MatSeqAIJGetArray(Cd, &cd));
  PetscCall(MatSeqAIJGetArray(Co, &co));

#if defined(_OPENMP)
  #pragma omp parallel
#endif
  {
    const PetscInt ldRmax = c->mxj * bs, ldOmax = c->mxJ * bs;
    PetscScalar   *R, *O;
    PetscInt      *jmark, *Jmark, k;
    int            mine;
    /* thread-private work (PETSc's malloc is not thread-safe) */
    R     = (PetscScalar *)malloc(sizeof(PetscScalar) * (size_t)(ldRmax * bs + 1));
    O     = (PetscScalar *)malloc(sizeof(PetscScalar) * (size_t)(ldOmax * bs + 1));
    jmark = (PetscInt *)malloc(sizeof(PetscInt) * (size_t)(c->nlb + c->ngb + 1));
    Jmark = (PetscInt *)malloc(sizeof(PetscInt) * (size_t)(c->nlc + 1));
    mine  = (R && O && jmark && Jmark);
    if (mine) {
      for (k = 0; k < c->nlb + c->ngb; k++) jmark[k] = -1;
      for (k = 0; k < c->nlc; k++) Jmark[k] = -1;
    } else {
#if defined(_OPENMP)
      #pragma omp atomic write
#endif
      fail = 1;
    }
    /* C = 0: the owned partial rows are written below, received ones added */
#if defined(_OPENMP)
    #pragma omp for schedule(static)
#endif
    for (PetscInt r = 0; r < ncr; r++) {
      memset(cd + cia[r], 0, sizeof(PetscScalar) * (size_t)(cia[r + 1] - cia[r]));
      memset(co + coia[r], 0, sizeof(PetscScalar) * (size_t)(coia[r + 1] - coia[r]));
    }
#if defined(_OPENMP)
    #pragma omp for schedule(dynamic, 1)
#endif
    for (PetscInt t = 0; t < c->nt; t++) {
      const PetscInt j0 = c->jp[t], nj = c->jp[t + 1] - j0, J0 = c->Jp[t], nJ = c->Jp[t + 1] - J0;
      const PetscInt ldR = nj * bs, ldO = nJ * bs;
      PetscInt       q, m, kk;
      if (!mine) continue;
      for (q = 0; q < nj; q++) jmark[c->jl[j0 + q]] = q * bs;
      for (q = 0; q < nJ; q++) Jmark[c->Jl[J0 + q]] = q * bs;
      memset(R, 0, sizeof(PetscScalar) * (size_t)(ldR * bs));
      memset(O, 0, sizeof(PetscScalar) * (size_t)(ldO * bs));
      /* R = sum_i Ps(i, I) A(i, :); a row's columns come as whole slot
         blocks (bs consecutive columns of one fine block) wherever all slots
         couple: one lookup and a vector axpy per block, entry by entry
         elsewhere */
      for (q = c->cp[t]; q < c->cp[t + 1]; q++) {
        const PetscInt    i = c->ci[q];
        const PetscScalar p = c->cv[q];
        for (m = 0; m < bs; m++) {
          const PetscInt r  = i * bs + m;
          PetscScalar   *Rm = R + m * ldR;
          for (kk = ia[r]; kk < ia[r + 1];) {
            const PetscInt col = ja[kk];
            if (c->dcn[col] == 0 && kk + bs <= ia[r + 1] && ja[kk + bs - 1] == col + bs - 1) {
              PetscScalar       *dst = Rm + jmark[c->dcb[col]];
              const PetscScalar *src = ad + kk;
              PetscInt           x;
#if defined(_OPENMP)
              #pragma omp simd
#endif
              for (x = 0; x < bs; x++) dst[x] += p * src[x];
              kk += bs;
            } else {
              Rm[jmark[c->dcb[col]] + c->dcn[col]] += p * ad[kk];
              kk++;
            }
          }
          for (kk = oia[r]; kk < oia[r + 1];) {
            const PetscInt col = oja[kk];
            if (c->ofb[col] && kk + bs <= oia[r + 1] && oja[kk + bs - 1] == col + bs - 1) {
              PetscScalar       *dst = Rm + jmark[c->ocb[col]];
              const PetscScalar *src = ao + kk;
              PetscInt           x;
#if defined(_OPENMP)
              #pragma omp simd
#endif
              for (x = 0; x < bs; x++) dst[x] += p * src[x];
              kk += bs;
            } else {
              Rm[jmark[c->ocb[col]] + c->ocn[col]] += p * ao[kk];
              kk++;
            }
          }
        }
      }
      /* O = R Ps, block by block */
      for (q = 0; q < nj; q++) {
        const PetscInt j = c->jl[j0 + q];
        for (kk = c->pp[j]; kk < c->pp[j + 1]; kk++) {
          const PetscScalar w  = c->pv[kk];
          const PetscInt    oo = Jmark[c->pl[kk]];
          for (m = 0; m < bs; m++) {
            const PetscScalar *src = R + m * ldR + q * bs;
            PetscScalar       *dst = O + m * ldO + oo;
            PetscInt           x;
#if defined(_OPENMP)
            #pragma omp simd
#endif
            for (x = 0; x < bs; x++) dst[x] += w * src[x];
          }
        }
      }
      /* out: owned rows into C's arrays, remote ones into the buffer */
      if (c->rb[t] < 0) {
        const PetscInt r0 = (c->tI[t] - c->cstart) * bs;
        for (m = 0; m < bs; m++)
          for (q = 0; q < nJ; q++) {
            const PetscInt off = c->Jo[J0 + q];
            PetscScalar   *dst = off >= 0 ? cd + cia[r0 + m] + off : co + coia[r0 + m] + (-1 - off);
            memcpy(dst, O + m * ldO + q * bs, sizeof(PetscScalar) * (size_t)bs);
          }
      } else {
        memcpy(c->rv + c->rb[t], O, sizeof(PetscScalar) * (size_t)(ldO * bs));
      }
      for (q = 0; q < nj; q++) jmark[c->jl[j0 + q]] = -1;
      for (q = 0; q < nJ; q++) Jmark[c->Jl[J0 + q]] = -1;
    }
    free(R); free(O); free(jmark); free(Jmark);
  }
  PetscCheck(!fail, PETSC_COMM_SELF, PETSC_ERR_MEM, "galerkin_blk: no memory for the thread work arrays");

  /* partial rows owned elsewhere: raw blocks to their owners, added there
     in rank order at the offsets of the symbolic phase */
  PetscCall(PetscLogEventBegin(gblk_ev_asm, A, 0, 0, 0));
  {
    MPI_Comm     comm = PetscObjectComm((PetscObject)C);
    MPI_Request *req;
    PetscMPIInt  i;
    PetscCall(PetscMalloc1(c->nsd + c->nrc + 1, &req));
    for (i = 0; i < c->nrc; i++) {
      PetscMPIInt cnt;
      PetscCall(PetscMPIIntCast(c->rco[i + 1] - c->rco[i], &cnt));
      PetscCallMPI(MPI_Irecv(c->rcv + c->rco[i], cnt, MPIU_SCALAR, c->rcr[i], c->tag, comm, &req[i]));
    }
    for (i = 0; i < c->nsd; i++) {
      PetscMPIInt cnt;
      PetscCall(PetscMPIIntCast(c->sdo[i + 1] - c->sdo[i], &cnt));
      PetscCallMPI(MPI_Isend(c->rv + c->sdo[i], cnt, MPIU_SCALAR, c->sdr[i], c->tag, comm, &req[c->nrc + i]));
    }
    PetscCallMPI(MPI_Waitall(c->nsd + c->nrc, req, MPI_STATUSES_IGNORE));
    PetscCall(PetscFree(req));
  }
  /* one thread per target row, its received rows in rank order */
#if defined(_OPENMP)
  #pragma omp parallel for schedule(dynamic, 4)
#endif
  for (PetscInt g = 0; g < c->nig; g++)
    for (PetscInt kq = c->igp[g]; kq < c->igp[g + 1]; kq++) {
      const PetscInt     k = c->igk[kq], r0 = c->inI[k] * bs, nJ = c->inp[k + 1] - c->inp[k], ldO = nJ * bs;
      const PetscScalar *src = c->rcv + c->ino[k];
      for (PetscInt m = 0; m < bs; m++)
        for (PetscInt q = 0; q < nJ; q++) {
          const PetscInt     off = c->inJo[c->inp[k] + q];
          PetscScalar       *dst = off >= 0 ? cd + cia[r0 + m] + off : co + coia[r0 + m] + (-1 - off);
          const PetscScalar *s_  = src + m * ldO + q * bs;
          for (PetscInt x = 0; x < bs; x++) dst[x] += s_[x];
        }
    }
  PetscCall(MatSeqAIJRestoreArray(Cd, &cd));
  PetscCall(MatSeqAIJRestoreArray(Co, &co));
  PetscCall(MatSeqAIJRestoreArrayRead(Ad, &ad));
  PetscCall(MatSeqAIJRestoreArrayRead(Ao, &ao));
  PetscCall(MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY));
  PetscCall(PetscLogEventEnd(gblk_ev_asm, A, 0, 0, 0));
  PetscCall(MatGetNonzeroState(C, &c->cnz));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* A reuse fits the symbolic phase: the same Ps and C (pattern unchanged), A
   of the same pattern (the same Mat in the same nonzero state, or one of equal
   sizes and pattern hash). O(1) unless A is another Mat. */
static PetscErrorCode gblk_fits(Mat A, Mat Ps, PetscInt bs, Mat C, gblk_ctx *c, PetscBool *fits)
{
  PetscBool        mpi;
  PetscObjectId    aid, pid, cid;
  PetscObjectState anz, cnz;
  Mat              Ad, Ao;
  const PetscInt  *ia, *oia;
  PetscInt         nA, mA, rs, cs, nAo;

  PetscFunctionBegin;
  *fits = PETSC_FALSE;
  PetscCall(PetscObjectTypeCompare((PetscObject)A, MATMPIAIJ, &mpi));
  PetscCall(PetscObjectGetId((PetscObject)Ps, &pid));
  PetscCall(PetscObjectGetId((PetscObject)C, &cid));
  PetscCall(MatGetNonzeroState(C, &cnz));
  if (!mpi || bs != c->bs || pid != c->pid || cid != c->cid || cnz != c->cnz) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(MatGetLocalSize(A, &nA, &mA));
  PetscCall(MatGetOwnershipRange(A, &rs, NULL));
  PetscCall(MatGetOwnershipRangeColumn(A, &cs, NULL));
  PetscCall(MatMPIAIJGetSeqAIJ(A, &Ad, &Ao, NULL));
  PetscCall(MatGetSize(Ao, NULL, &nAo));
  PetscCall(csr(Ad, &ia, NULL));
  PetscCall(csr(Ao, &oia, NULL));
  if (nA != c->nA || mA != nA || cs != rs || nAo != c->nAo || ia[nA] != c->nzd || oia[nA] != c->nzo) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscObjectGetId((PetscObject)A, &aid));
  PetscCall(MatGetNonzeroState(A, &anz));
  if (aid != c->aid || anz != c->anz) {
    uint64_t  h;
    PetscBool sorted;
    PetscCall(pattern_scan(A, &h, &sorted));
    if (h != c->hash) PetscFunctionReturn(PETSC_SUCCESS);
    c->aid = aid; c->anz = anz;
  }
  *fits = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode jorek_galerkin_blk(Mat A, Mat Ps, int bs, int reuse, Mat *C, int *ok)
{
  MPI_Comm  comm = PetscObjectComm((PetscObject)A);
  PetscBool mpi, pmpi, sorted = PETSC_FALSE;
  gblk_ctx *c = NULL;
  uint64_t  hash = 0;
  int       good;

  PetscFunctionBegin;
  *ok = 0;
  if (gblk_ev_sym < 0) {
    PetscClassId cid;
    PetscCall(PetscClassIdRegister("JOREK galerkin", &cid));
    PetscCall(PetscLogEventRegister("GMG_GalSym", cid, &gblk_ev_sym));
    PetscCall(PetscLogEventRegister("GMG_GalNum", cid, &gblk_ev_num));
    PetscCall(PetscLogEventRegister("GMG_GalAsm", cid, &gblk_ev_asm));
  }
  if (reuse) {
    PetscContainer cont;
    PetscBool      fits = PETSC_FALSE;
    int            st   = 0; /* 0 fits, 1 C not made here, 2 A's pattern / Ps / C changed */
    PetscCall(PetscObjectQuery((PetscObject)*C, GBLK_KEY, (PetscObject *)&cont));
    if (!cont) st = 1;
    else {
      PetscCall(PetscContainerGetPointer(cont, (void **)&c));
      PetscCall(gblk_fits(A, Ps, bs, *C, c, &fits));
      if (!fits) st = 2;
    }
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &st, 1, MPI_INT, MPI_MAX, comm));
    if (st == 1) PetscFunctionReturn(PETSC_SUCCESS);
    PetscCheck(st == 0, comm, PETSC_ERR_ARG_INCOMP, "galerkin_blk: reuse, but A's pattern, Ps or C changed since the symbolic phase");
    PetscCall(PetscLogEventBegin(gblk_ev_num, A, 0, 0, 0));
    PetscCall(gblk_numeric(A, *C, c));
    PetscCall(PetscLogEventEnd(gblk_ev_num, A, 0, 0, 0));
    *ok = 1;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  PetscCall(PetscObjectTypeCompare((PetscObject)A, MATMPIAIJ, &mpi));
  PetscCall(PetscObjectTypeCompare((PetscObject)Ps, MATMPIAIJ, &pmpi));
  good = mpi && pmpi && bs > 1;
  if (good) {
    PetscInt nA, mA, rs, cs, np;
    PetscCall(MatGetLocalSize(A, &nA, &mA));
    PetscCall(MatGetOwnershipRange(A, &rs, NULL));
    PetscCall(MatGetOwnershipRangeColumn(A, &cs, NULL));
    PetscCall(MatGetLocalSize(Ps, &np, NULL));
    good = (nA % bs == 0) && (rs % bs == 0) && (np * bs == nA) && (mA == nA) && (cs == rs);
  }
  if (good) {
    PetscCall(pattern_scan(A, &hash, &sorted));
    good = (sorted == PETSC_TRUE);
  }
  PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &good, 1, MPI_INT, MPI_LAND, comm));
  if (!good) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscLogEventBegin(gblk_ev_sym, A, 0, 0, 0));
  PetscCall(gblk_symbolic(A, Ps, bs, C, &c));
  PetscCall(PetscLogEventEnd(gblk_ev_sym, A, 0, 0, 0));
  c->hash = hash;
  PetscCall(PetscObjectGetId((PetscObject)A, &c->aid));
  PetscCall(MatGetNonzeroState(A, &c->anz));
  PetscCall(PetscObjectGetId((PetscObject)Ps, &c->pid));
  PetscCall(PetscObjectGetId((PetscObject)*C, &c->cid));
  PetscCall(PetscObjectContainerCompose((PetscObject)*C, GBLK_KEY, c, ctx_destroy));
  PetscCall(PetscLogEventBegin(gblk_ev_num, A, 0, 0, 0));
  PetscCall(gblk_numeric(A, *C, c));
  PetscCall(PetscLogEventEnd(gblk_ev_num, A, 0, 0, 0));
  *ok = 1;
  PetscFunctionReturn(PETSC_SUCCESS);
}
#endif
