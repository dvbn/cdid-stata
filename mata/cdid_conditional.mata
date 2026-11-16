*! Native streaming engine of the chained conditional pretrend CvM test, current StataNow 19.5.
*! Pinned R cdid 1.1.0.0 (commit ):
*! R/chained_conditional_engine.R (pure computational contract only). No RNG, no grid
*! construction from data, no nuisance refit and no public frontend live here: the caller
*! supplies the full original-F grid design, the cluster map, the cell spec (s, D, sigma,
*! coordinates), bounded extractor pieces (nu, rho, P, xi, blam over possibly repeated
*! original-row indices) and the common B x C Rademacher sign batches.
*!   phi_i(u) = h_u(z_i) A_i - sum_r R_r(u) kappa_r,i + xi_i' L(u)
*!   A_i = sum_r s_r nu_r,i rho_r,i / D_r,  kappa_r,i = s_r nu_r,i / D_r
*!   R_r(u) = H_N[nu_r h_u rho_r] / D_r,  L(u) = H_N[h_u P] - sum_r (s_r/D_r) R_r(u) blam_r
*!   J(u) = sigma sum_r s_r R_r(u);  CvM = (N/M) sum_c sum_m J_c(u_m)^2 (M = N grid rows)
*!   J*_b(u) = N^-1 sum_k eta_bk S_k(u);  CvM_b analogous;  p = mean(CvM_b > CvM)
*! N is the supplied full-frame normalization; nothing is centered.
version 19.5
mata:
/* Same (start,end) rows as panelsetup(V,idcol) for runs of equal keys, vectorized. */
real matrix _cdid_ps_cond(transmorphic matrix V, real scalar idcol)
{
    real scalar n
    real colvector s
    n=rows(V)
    if (n<2 | cols(V)<idcol) return(panelsetup(V,idcol))
    s=1\(1:+vec(selectindex(V[|2,idcol\n,idcol|]:!=V[|1,idcol\n-1,idcol|])))
    if (rows(s)==1) return((1,n))
    return((s,((s[|2\rows(s)|]:-1)\n)))
}
struct cdid_cpt_budget {
    real scalar max_units, max_clusters, max_cells, max_draws
    real scalar max_distinct_points, max_terms, max_coords
    real scalar max_piece_rows, max_events
    real scalar max_sign_batch_rows, max_sign_batch_cells
    real scalar grid_chunk, event_chunk
    real scalar max_workspace_cells, max_workspace_bytes, max_work_ops
}
struct cdid_cpt_grid {
    real scalar M, n_distinct, p
    real matrix zd
    real colvector first, map, mult
    pointer(real matrix) scalar pH
}
struct cdid_cpt_spec {
    string scalar id
    real scalar sigma, R, Dl
    real colvector s, D
    string colvector coords
}
struct cdid_cpt_acc {
    struct cdid_cpt_budget scalar budget
    struct cdid_cpt_spec scalar spec
    real scalar N, C, Md, keep_rows, finished, n_rows_seen, n_pieces, ev_n, ev_cap
    real matrix Vd, Pd, Bl, AggX, AggKu, rows_nu, rows_xi
    real colvector sum_nu, ev_key, ev_val, rows_A
}
struct cdid_cpt_cell {
    string scalar id
    real scalar sigma, M, n_distinct, n_terms, n_coordinates, n_events, mu1, has_rows, A_d_uniq
    real colvector J, mult, R1, L1, A_k, A_d, A_v, A_d_ord, D_supplied, D_accumulated, rows_A, rows_map
    real matrix zd, Rmat, Lmat, AggX, AggK, rows_kappa, rows_xi
    string colvector coords
    pointer(real matrix) scalar pH
}
struct cdid_cpt_mult_state {
    struct cdid_cpt_budget scalar budget
    real scalar N, C, B, n_cells, batch_rows, cursor, finished
    real colvector CvMb
}
struct cdid_cpt_result {
    real scalar N, C, B, n_cells, alpha, CvM, CvMcval, CvMpval, process_nonzero
    real colvector CvMb
    string scalar status, reason
}
void _cdid_cpt_fail(string scalar message, real scalar rc)
{
    errprintf("%s\n", message)
    exit(rc)
}
real colvector _cdid_cpt_index(real matrix predicate)
{
    return(vec(selectindex(predicate)))
}
real scalar _cdid_cpt_whole(real scalar x, real scalar lo)
{
    return(x<. & x==floor(x) & x>=lo & x<=2^53-1)
}
real scalar _cdid_cpt_charge(real scalar a, real scalar b, real scalar cap, real scalar total)
{
    if(a>0 & b>(2^53-1)/a)
        _cdid_cpt_fail("Conditional-test allocation product exceeds exact integer range.",908)
    if(a*b>cap-total)
        _cdid_cpt_fail("Conditional-test numeric storage exceeds the explicit workspace budget.",908)
    return(total+a*b)
}
/* The work cap is opt-in: max_work_ops=. means no counted-work limit. */
void _cdid_cpt_work(struct cdid_cpt_budget scalar b, real scalar ops)
{
    if(b.max_work_ops<. & ops>b.max_work_ops)
        _cdid_cpt_fail("Counted conditional-test work exceeds the explicit opt-in work cap.",908)
}
struct cdid_cpt_budget scalar cdid_cpt_budget_create(real scalar max_units,
    real scalar max_clusters, real scalar max_cells, real scalar max_draws,
    real scalar max_distinct_points, real scalar max_terms, real scalar max_coords,
    real scalar max_piece_rows, real scalar max_events,
    real scalar max_sign_batch_rows, real scalar max_sign_batch_cells,
    real scalar grid_chunk, real scalar event_chunk,
    real scalar max_workspace_cells, real scalar max_workspace_bytes, real scalar max_work_ops)
{
    struct cdid_cpt_budget scalar b
    real rowvector values
    values=(max_units,max_clusters,max_cells,max_draws,max_distinct_points,max_terms,max_coords,
        max_piece_rows,max_events,max_sign_batch_rows,max_sign_batch_cells,grid_chunk,event_chunk,
        max_workspace_cells,max_workspace_bytes)
    if(missing(values) | any(values:<1) | any(values:!=floor(values)) | any(values:>2^53-1))
        _cdid_cpt_fail("Every conditional-test budget is an explicit positive finite integer <=2^53-1.",498)
    if(max_work_ops<. & !_cdid_cpt_whole(max_work_ops,1))
        _cdid_cpt_fail("max_work_ops is either . (no counted-work cap) or an explicit positive integer.",498)
    b.max_units=max_units; b.max_clusters=max_clusters; b.max_cells=max_cells; b.max_draws=max_draws
    b.max_distinct_points=max_distinct_points; b.max_terms=max_terms; b.max_coords=max_coords
    b.max_piece_rows=max_piece_rows; b.max_events=max_events
    b.max_sign_batch_rows=max_sign_batch_rows; b.max_sign_batch_cells=max_sign_batch_cells
    b.grid_chunk=grid_chunk; b.event_chunk=event_chunk
    b.max_workspace_cells=max_workspace_cells; b.max_workspace_bytes=max_workspace_bytes
    b.max_work_ops=max_work_ops
    return(b)
}
/* Sum of the rows of X by group into target rows: target[g,] += sum_{group==g} X. */
real matrix _cdid_cpt_rowsum_into(real matrix target, real matrix X, real colvector group)
{
    real colvector o, g, ids
    real matrix info, S
    if(rows(X)<1 | cols(X)<1) return(target)
    o=order(group,1)
    g=group[o]
    info=_cdid_ps_cond(g,1)
    if(rows(info)==rows(X)) {
        target[group,.]=target[group,.]+X
        return(target)
    }
    S=panelsum(X[o,.],info)
    ids=g[info[.,1]]
    target[ids,.]=target[ids,.]+S
    return(target)
}
/* H[i,m] = 1{zrow_i <= zcol_m componentwise}; all ones when there is no covariate. */
real matrix cdid_cpt_indicator(real matrix zrow, real matrix zcol)
{
    real matrix H
    real scalar j, p, r, c, has_h, rmin, rmax, cmin, cmax
    real colvector zr
    real rowvector zc
    p=cols(zrow)
    if(cols(zcol)!=p) _cdid_cpt_fail("Indicator designs must share their columns.",498)
    r=rows(zrow); c=rows(zcol)
    if(r==0 | c==0 | p==0) return(J(r,c,1))
    has_h=0
    for(j=1;j<=p;j++) {
        zr=zrow[.,j]; zc=zcol[.,j]'
        rmax=colmax(zr); cmin=rowmin(zc)
        if(rmax<. & cmin<. & rmax<=cmin) continue
        rmin=colmin(zr); cmax=rowmax(zc)
        if(rmin<. & cmax<. & rmin>cmax) return(J(r,c,0))
        if(!has_h) {
            H=(zr:<=J(r,1,zc))
            has_h=1
        }
        else H=H:*(zr:<=J(r,1,zc))
    }
    if(!has_h) return(J(r,c,1))
    return(H)
}
pointer(real matrix) scalar _cdid_cpt_alloc_H(real matrix H)
{
    return(&H)
}
/* Q*H for a 0/1 dominance matrix H. Mata products skip zero entries of the left operand,
   so the product is formed as (H'Q')' whenever H carries fewer nonzero products than Q;
   both orientations accumulate every entry over the same ascending inner index. */
real matrix _cdid_cpt_mulH(real matrix Q, real matrix H, real scalar nnzQ)
{
    real scalar r, s, dens
    real rowvector sel
    r=rows(H)
    if(r<2 | rows(Q)<1 | cols(H)<1) return(Q*H)
    s=min((r,64))
    sel=ceil((1..s):*(r/s))
    dens=sum(H[sel,.]:!=0)/(s*cols(H))
    if(1.25*dens*r*rows(Q)<nnzQ) return(cross(H,Q')')
    return(Q*H)
}
void _cdid_cpt_ensure_grid_H(struct cdid_cpt_grid scalar grid, struct cdid_cpt_budget scalar budget)
{
    real scalar Md, cap
    Md=grid.n_distinct
    if(grid.pH==NULL) grid.pH=_cdid_cpt_alloc_H(J(0,0,.))
    if(rows(*grid.pH)==Md & cols(*grid.pH)==Md) return
    cap=min((budget.max_workspace_cells,floor(budget.max_workspace_bytes/8)))
    if(Md>=1 & Md<=4096 & Md*Md<=cap/4) *grid.pH=cdid_cpt_indicator(grid.zd,grid.zd)
}
/* Distinct grid points by exact binary64 row equality (-0 = 0), in first-occurrence
   order of the original units; every unit keeps its map and multiplicity (M rows). */
struct cdid_cpt_grid scalar cdid_cpt_grid_create(real matrix Z, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_grid scalar g
    real scalar M, p, G, cap, j
    real matrix W, info
    real colvector grp, orig, firstg, ord, rank, multg, is_new
    M=rows(Z); p=cols(Z)
    if(M<1 | M>2^31-1 | missing(Z))
        _cdid_cpt_fail("The conditional grid design must be finite on every original unit.",498)
    if(M>budget.max_units)
        _cdid_cpt_fail("Grid rows exceed the explicit max_units budget.",908)
    cap=min((budget.max_workspace_cells,floor(budget.max_workspace_bytes/8)))
    (void) _cdid_cpt_charge(M,3*p+9,cap,0)
    g.M=M; g.p=p; g.pH=NULL
    if(p==0) {
        g.n_distinct=1; g.first=1; g.map=J(M,1,1); g.mult=M; g.zd=J(1,0,0)
        return(g)
    }
    /* Sorting on the design AND the original index makes the first sorted row of
       each run the first occurrence; -0 is mapped onto +0 for the sort keys only. */
    W=sort((Z:+0,(1..M)'),(1..(p+1)))
    orig=W[.,p+1]
    if(M>1) {
        is_new=J(M-1,1,0)
        for(j=1;j<=p;j++) is_new=is_new:|(W[(2..M),j]:!=W[(1..(M-1)),j])
        grp=1\(1:+runningsum(is_new))
    }
    else grp=1
    G=grp[M]
    if(G==M) {
        g.n_distinct=M
        g.first=(1..M)'
        g.mult=J(M,1,1)
        g.map=(1..M)'
        g.zd=Z
        return(g)
    }
    info=_cdid_ps_cond(grp,1)
    firstg=orig[info[.,1]]
    multg=info[.,2]-info[.,1]:+1
    ord=order(firstg,1)
    rank=J(G,1,0); rank[ord]=(1..G)'
    g.n_distinct=G
    g.first=firstg[ord]
    g.mult=multg[ord]
    g.map=J(M,1,0); g.map[orig]=rank[grp]
    g.zd=Z[g.first,.]
    return(g)
}
/* One declared cell: id, orientation sigma (the cell's +/-1 sign; the multiplier draws omit
   it exactly as the source does, so only a unit modulus keeps CvM and CvM_b consistent),
   term signs s, supplied denominators D>0 and the coordinate labels of the nuisance IF
   columns (Dl may be zero). */
struct cdid_cpt_spec scalar cdid_cpt_spec_create(string scalar id, real scalar sigma,
    real colvector s, real colvector D, string colvector coords)
{
    struct cdid_cpt_spec scalar sp
    real scalar R
    R=rows(s)
    if(R<1 | cols(s)!=1 | rows(D)!=R | cols(D)!=1 | missing(s) | missing(D) | any(D:<=0) |
        (sigma!=1 & sigma!=-1) | cols(coords)>1 | R>2^31-1)
        _cdid_cpt_fail("A cell spec needs sigma in {-1,1}, R>=1 finite term signs, R positive finite denominators and a string column of coordinates.",498)
    sp.id=id; sp.sigma=sigma; sp.R=R; sp.Dl=rows(coords)
    sp.s=s; sp.D=D; sp.coords=coords
    return(sp)
}
/* Whole-buffer numeric envelope of one cell build (state, one piece, finalize scratch) and
   its counted work, checked before any allocation. Numeric payload, not RSS or time. */
real scalar _cdid_cpt_preflight_cell(struct cdid_cpt_budget scalar b, real scalar N,
    real scalar C, real scalar Md, real scalar p, real scalar R, real scalar Dl, real scalar keep_rows)
{
    real scalar cap, total, pr, gc
    if(N>b.max_units | C>b.max_clusters | Md>b.max_distinct_points | R>b.max_terms | Dl>b.max_coords)
        _cdid_cpt_fail("Conditional cell dimensions exceed explicit caller limits.",908)
    cap=min((b.max_workspace_cells,floor(b.max_workspace_bytes/8)))
    pr=min((N,b.max_piece_rows)); gc=min((Md,b.grid_chunk))
    total=0
    total=_cdid_cpt_charge(Md,p+4+3*(R+Dl),cap,total)
    total=_cdid_cpt_charge(N,2,cap,total)
    total=_cdid_cpt_charge(C,2*(R+Dl),cap,total)
    total=_cdid_cpt_charge(2*Dl+3,R,cap,total)
    total=_cdid_cpt_charge(b.max_events,6,cap,total)
    total=_cdid_cpt_charge(pr,3*R+2*Dl+8,cap,total)
    total=_cdid_cpt_charge(Md,gc,cap,total)
    total=_cdid_cpt_charge(gc,2*(R+Dl)+1,cap,total)
    if(keep_rows) total=_cdid_cpt_charge(N,2*(1+R+Dl)+1,cap,total)
    _cdid_cpt_work(b,64*N*(R+Dl)+16*Md*Md*(p+R+Dl))
    return(total)
}
/* Grid rows are the full original frame: exactly one grid row per unit (M = N), cluster
   positions 1..C for every unit (zero-IF clusters and empty declared clusters kept). */
struct cdid_cpt_acc scalar cdid_cpt_cell_start(struct cdid_cpt_spec scalar spec,
    struct cdid_cpt_grid scalar grid, real colvector cmap, real scalar C, real scalar N,
    real scalar keep_rows, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_acc scalar a
    real scalar entries, R, Dl, Md
    R=spec.R; Dl=spec.Dl; Md=grid.n_distinct
    if(!_cdid_cpt_whole(N,1) | N>2^31-1 | !_cdid_cpt_whole(C,1) | C>2^31-1 | grid.M!=N |
        rows(cmap)!=N | cols(cmap)!=1 | missing(cmap) | any(cmap:!=floor(cmap)) |
        any(cmap:<1) | any(cmap:>C) | (keep_rows!=0 & keep_rows!=1))
        _cdid_cpt_fail("Require one grid row per original unit, full-F cluster positions in 1..C and a 0/1 keep_rows flag.",498)
    if(C>(2^53-1)/Md)
        _cdid_cpt_fail("Declared clusters times distinct grid points exceed the exact event-key range.",908)
    entries=_cdid_cpt_preflight_cell(budget,N,C,Md,grid.p,R,Dl,keep_rows)
    a.budget=budget; a.spec=spec; a.N=N; a.C=C; a.Md=Md; a.keep_rows=keep_rows
    a.finished=0; a.n_rows_seen=0; a.n_pieces=0; a.ev_n=0; a.ev_cap=0
    a.Vd=J(Md,R,0); a.Pd=J(Md,Dl,0); a.Bl=J(Dl,R,0)
    a.AggX=J(C,Dl,0); a.AggKu=J(C,R,0); a.sum_nu=J(R,1,0)
    a.ev_key=J(0,1,0); a.ev_val=J(0,1,0)
    if(keep_rows) {
        a.rows_A=J(N,1,0); a.rows_nu=J(N,R,0); a.rows_xi=J(N,Dl,0)
    }
    else {
        a.rows_A=J(0,1,0); a.rows_nu=J(0,0,0); a.rows_xi=J(0,0,0)
    }
    return(a)
}
void _cdid_cpt_events_append(struct cdid_cpt_acc scalar a, real colvector keys, real colvector vals)
{
    real scalar need, newcap
    real colvector tmp
    need=a.ev_n+rows(keys)
    if(need>a.budget.max_events)
        _cdid_cpt_fail("Sparse grid-by-cluster events exceed the explicit max_events budget.",908)
    if(need>a.ev_cap) {
        newcap=min((a.budget.max_events,max((2*a.ev_cap,need,1024))))
        tmp=J(newcap,1,0)
        if(a.ev_n) tmp[(1..a.ev_n)]=a.ev_key[(1..a.ev_n)]
        a.ev_key=tmp
        tmp=J(newcap,1,0)
        if(a.ev_n) tmp[(1..a.ev_n)]=a.ev_val[(1..a.ev_n)]
        a.ev_val=tmp
        a.ev_cap=newcap
    }
    a.ev_key[((a.ev_n+1)..need)]=keys
    a.ev_val[((a.ev_n+1)..need)]=vals
    a.ev_n=need
}
/* One extractor piece: rows are original-F positions (repeats allowed, layout rows of one
   unit are summed); nu/rho are rows x R, P/xi rows x Dl, blam Dl x R (piece sums). */
void cdid_cpt_cell_add(struct cdid_cpt_acc scalar a, struct cdid_cpt_grid scalar grid,
    real colvector cmap, real colvector rows, real matrix nu, real matrix rho,
    real matrix P, real matrix xi, real matrix blam)
{
    real scalar n, R, Dl
    real colvector d, k, A, nz, o, g, ids
    real matrix nr, info
    if(a.finished) _cdid_cpt_fail("The cell accumulator is already finalized.",498)
    if(grid.M!=a.N | grid.n_distinct!=a.Md | rows(cmap)!=a.N)
        _cdid_cpt_fail("Piece grid or cluster map does not match the started cell.",498)
    n=rows(rows)
    if(n==0) {
        a.n_pieces=a.n_pieces+1
        return
    }
    R=a.spec.R; Dl=a.spec.Dl
    if(cols(rows)!=1 | missing(rows) | any(rows:!=floor(rows)) | any(rows:<1) | any(rows:>a.N))
        _cdid_cpt_fail("Piece rows must be original-F positions in 1..N.",498)
    if(n>a.budget.max_piece_rows)
        _cdid_cpt_fail("Piece rows exceed the explicit max_piece_rows budget.",908)
    if(rows(nu)!=n | cols(nu)!=R | rows(rho)!=n | cols(rho)!=R | rows(P)!=n | cols(P)!=Dl |
        rows(xi)!=n | cols(xi)!=Dl | rows(blam)!=Dl | cols(blam)!=R |
        missing(nu) | missing(rho) | missing(P) | missing(xi) | missing(blam))
        _cdid_cpt_fail("Process extractor returned malformed or nonfinite rows for cell "+a.spec.id+".",498)
    d=grid.map[rows]; k=cmap[rows]
    nr=nu:*rho
    o=order(d,1); g=d[o]; info=_cdid_ps_cond(g,1)
    if(rows(info)==n) {
        a.Vd[d,.]=a.Vd[d,.]+nr
        if(Dl) a.Pd[d,.]=a.Pd[d,.]+P
    }
    else {
        ids=g[info[.,1]]
        a.Vd[ids,.]=a.Vd[ids,.]+panelsum(nr[o,.],info)
        if(Dl) a.Pd[ids,.]=a.Pd[ids,.]+panelsum(P[o,.],info)
    }
    a.Bl=a.Bl+blam
    o=order(k,1); g=k[o]; info=_cdid_ps_cond(g,1)
    if(rows(info)==n) {
        if(Dl) a.AggX[k,.]=a.AggX[k,.]+xi
        a.AggKu[k,.]=a.AggKu[k,.]+nu
    }
    else {
        ids=g[info[.,1]]
        if(Dl) a.AggX[ids,.]=a.AggX[ids,.]+panelsum(xi[o,.],info)
        a.AggKu[ids,.]=a.AggKu[ids,.]+panelsum(nu[o,.],info)
    }
    a.sum_nu=a.sum_nu+colsum(nu)'
    A=nr*(a.spec.s:/a.spec.D)
    nz=_cdid_cpt_index(A:!=0)
    if(rows(nz)) _cdid_cpt_events_append(a,((k[nz]:-1):*a.Md):+d[nz],A[nz])
    a.n_rows_seen=a.n_rows_seen+n; a.n_pieces=a.n_pieces+1
    if(a.keep_rows) {
        a.rows_A=_cdid_cpt_rowsum_into(a.rows_A,A,rows)
        a.rows_nu=_cdid_cpt_rowsum_into(a.rows_nu,nu,rows)
        if(Dl) a.rows_xi=_cdid_cpt_rowsum_into(a.rows_xi,xi,rows)
    }
}
/* Process values at every distinct grid point in grid chunks, the scalar link (h = 1),
   the merged nonzero events and the compact multiplier state. */
struct cdid_cpt_cell scalar cdid_cpt_cell_finish(struct cdid_cpt_acc scalar a, struct cdid_cpt_grid scalar grid)
{
    struct cdid_cpt_cell scalar c
    real scalar R, Dl, Md, N, gc, st, en, ne, nv, has_H
    real colvector sD, cols, key, val, o, skey, sval, keep, vidx
    real matrix Vn, Vn_v, Pd_vt, zd_v, H, Rc, info
    if(a.finished) _cdid_cpt_fail("The cell accumulator is already finalized.",498)
    if(grid.M!=a.N | grid.n_distinct!=a.Md)
        _cdid_cpt_fail("Finalize grid does not match the started cell.",498)
    R=a.spec.R; Dl=a.spec.Dl; Md=a.Md; N=a.N
    sD=a.spec.s:/a.spec.D
    Vn=a.Vd:/(a.spec.D')
    c.Rmat=J(Md,R,0); c.Lmat=J(Dl,Md,0)
    vidx=_cdid_cpt_index((rowsum(Vn:!=0):>0):|(Dl ? (rowsum(a.Pd:!=0):>0) : J(Md,1,0)))
    nv=rows(vidx)
    if(nv) _cdid_cpt_ensure_grid_H(grid,a.budget)
    c.pH=grid.pH
    has_H=(grid.pH!=NULL ? (rows(*grid.pH)==Md & cols(*grid.pH)==Md) : 0)
    if(nv) {
        Vn_v=(nv==Md ? Vn : Vn[vidx,.])
        if(Dl) Pd_vt=(nv==Md ? a.Pd' : a.Pd[vidx,.]')
        if(!has_H) zd_v=(nv==Md ? grid.zd : grid.zd[vidx,.])
        gc=min((a.budget.grid_chunk,Md))
        if(has_H & gc==Md) {
            if(Dl) {
                H=_cdid_cpt_mulH((Vn'\a.Pd'),*grid.pH,sum((Vn,a.Pd):!=0))
                Rc=(H[(1..R),.]')/N
                c.Rmat=Rc
                c.Lmat=(H[((R+1)..(R+Dl)),.]-a.Bl*(Rc':*sD))/N
            }
            else {
                Rc=(_cdid_cpt_mulH(Vn',*grid.pH,sum(Vn:!=0))')/N
                c.Rmat=Rc
            }
        }
        else {
            for(st=1;st<=Md;st=st+gc) {
                en=min((st+gc-1,Md)); cols=(st..en)'
                if(has_H) H=(nv==Md ? (*grid.pH)[.,cols] : (*grid.pH)[vidx,cols])
                else H=cdid_cpt_indicator(zd_v,grid.zd[cols,.])
                Rc=(H'*Vn_v)/N
                c.Rmat[cols,.]=Rc
                if(Dl) c.Lmat[.,cols]=((Pd_vt*H)-a.Bl*(Rc':*sD))/N
            }
        }
    }
    c.J=a.spec.sigma:*(c.Rmat*a.spec.s)
    c.R1=(colsum(Vn)')/N
    if(Dl) c.L1=((colsum(a.Pd)')-a.Bl*(c.R1:*sD))/N
    else c.L1=J(0,1,0)
    ne=a.ev_n
    skey=J(0,1,0); sval=J(0,1,0)
    if(ne) {
        key=a.ev_key[(1..ne)]; val=a.ev_val[(1..ne)]
        o=order(key,1); key=key[o]; val=val[o]
        info=_cdid_ps_cond(key,1)
        sval=panelsum(val,info)
        skey=key[info[.,1]]
        keep=_cdid_cpt_index(sval:!=0)
        if(rows(keep)) {
            skey=skey[keep]; sval=sval[keep]
        }
        else {
            skey=J(0,1,0); sval=J(0,1,0)
        }
    }
    c.A_k=floor((skey:-1):/Md):+1
    c.A_d=skey-(c.A_k:-1):*Md
    c.A_v=sval
    c.n_events=rows(sval)
    if(c.n_events) {
        c.A_d_ord=order(c.A_d,1)
        c.A_d_uniq=(c.n_events==1 ? 1 : all(c.A_d[c.A_d_ord[(2..c.n_events)]]:!=c.A_d[c.A_d_ord[(1..(c.n_events-1))]]))
    }
    else {
        c.A_d_ord=J(0,1,0)
        c.A_d_uniq=1
    }
    c.id=a.spec.id; c.sigma=a.spec.sigma; c.M=grid.M; c.n_distinct=Md
    c.n_terms=R; c.n_coordinates=Dl; c.coords=a.spec.coords
    c.mult=grid.mult; c.zd=grid.zd
    c.AggX=a.AggX; c.AggK=a.AggKu:*(sD')
    c.mu1=sum(a.spec.s:*c.R1)
    c.D_supplied=a.spec.D; c.D_accumulated=a.sum_nu/N
    c.has_rows=a.keep_rows
    if(a.keep_rows) {
        c.rows_A=a.rows_A; c.rows_kappa=a.rows_nu:*(sD'); c.rows_xi=a.rows_xi; c.rows_map=grid.map
    }
    else {
        c.rows_A=J(0,1,0); c.rows_kappa=J(0,0,0); c.rows_xi=J(0,0,0); c.rows_map=J(0,1,0)
    }
    if(missing(c.Rmat) | missing(c.Lmat) | missing(c.J) | missing(c.R1) | missing(c.L1) | missing(c.AggK))
        _cdid_cpt_fail("Nonfinite conditional process state for cell "+c.id+".",504)
    a.finished=1
    return(c)
}
/* CvM = (N/M) sum_c sum_m J_c(u_m)^2 over all M grid rows (distinct points x multiplicity). */
real scalar cdid_cpt_statistic(struct cdid_cpt_cell colvector cells, real scalar N)
{
    real scalar i, out
    out=0
    for(i=1;i<=rows(cells);i++) out=out+(N/cells[i].M)*sum(cells[i].mult:*(cells[i].J:^2))
    return(out)
}
/* CvM_b contribution of one sign batch E (b x C) for one cell. */
real colvector _cdid_cpt_batch_cell(real matrix E, struct cdid_cpt_cell scalar cl, real scalar N,
    real scalar gchunk, real scalar echunk)
{
    real scalar b, Md, ne, st, en, nq, has_H
    real colvector out, ii, dsub, o, dd, qidx, cols, keep
    real matrix EX, EK, Q, Q_nz, zd_q, Em, Emt, S, info, H, Js
    b=rows(E); Md=cl.n_distinct
    out=J(b,1,0)
    EK=E*cl.AggK
    if(cl.n_coordinates) EX=E*cl.AggX
    ne=cl.n_events
    has_H=(cl.pH!=NULL ? (rows(*cl.pH)==Md & cols(*cl.pH)==Md) : 0)
    nq=0
    if(ne) {
        if(cl.A_d_uniq & rows(cl.A_d_ord)==ne) {
            o=cl.A_d_ord
            qidx=cl.A_d[o]
            nq=ne
            Q_nz=E[.,cl.A_k[o]]:*(cl.A_v[o]')
            if(!has_H) zd_q=cl.zd[qidx,.]
        }
        else if(echunk>=ne & rows(cl.A_d_ord)==ne) {
            o=cl.A_d_ord
            dd=cl.A_d[o]
            info=_cdid_ps_cond(dd,1)
            Emt=(E[.,cl.A_k[o]]:*(cl.A_v[o]'))'
            S=panelsum(Emt,info)
            dd=dd[info[.,1]]
            keep=_cdid_cpt_index(rowsum(S:!=0):>0)
            nq=rows(keep)
            if(nq) {
                if(nq==rows(dd)) {
                    qidx=dd
                    Q_nz=S'
                }
                else {
                    qidx=dd[keep]
                    Q_nz=S[keep,.]'
                }
                if(!has_H) zd_q=cl.zd[qidx,.]
            }
        }
        else {
            Q=J(b,Md,0)
            for(st=1;st<=ne;st=st+echunk) {
                en=min((ne,st+echunk-1)); ii=(st..en)'
                Em=E[.,cl.A_k[ii]]:*(cl.A_v[ii]')
                dsub=cl.A_d[ii]
                o=order(dsub,1)
                dd=dsub[o]
                info=_cdid_ps_cond(dd,1)
                if(rows(info)==rows(dsub)) {
                    Q[.,dsub]=Q[.,dsub]+Em
                }
                else {
                    Emt=Em'
                    S=panelsum(Emt[o,.],info)
                    dd=dd[info[.,1]]
                    Q[.,dd]=Q[.,dd]+S'
                }
            }
            qidx=_cdid_cpt_index(colsum(Q:!=0):>0)
            nq=rows(qidx)
            if(nq) {
                Q_nz=Q[.,qidx]
                if(!has_H) zd_q=cl.zd[qidx,.]
            }
        }
    }
    if(has_H & gchunk>=Md) {
        if(nq) {
            if(nq==Md) Js=_cdid_cpt_mulH(Q_nz,*cl.pH,rows(Q_nz)*cols(Q_nz))-EK*(cl.Rmat')
            else Js=_cdid_cpt_mulH(Q_nz,(*cl.pH)[qidx,.],rows(Q_nz)*cols(Q_nz))-EK*(cl.Rmat')
        }
        else Js=-EK*(cl.Rmat')
        if(cl.n_coordinates) Js=Js+EX*cl.Lmat
        Js=Js/N
        return((N/cl.M)*((Js:*Js)*cl.mult))
    }
    for(st=1;st<=Md;st=st+gchunk) {
        en=min((Md,st+gchunk-1)); cols=(st..en)'
        if(nq) {
            if(has_H) H=(nq==Md ? (*cl.pH)[.,cols] : (*cl.pH)[qidx,cols])
            else H=cdid_cpt_indicator(zd_q,cl.zd[cols,.])
            Js=Q_nz*H-EK*(cl.Rmat[cols,.]')
        }
        else Js=-EK*(cl.Rmat[cols,.]')
        if(cl.n_coordinates) Js=Js+EX*cl.Lmat[.,cols]
        Js=Js/N
        out=out+(N/cl.M)*((Js:*Js)*cl.mult[cols,1])
    }
    return(out)
}
/* Actual cluster process rows streamed over grid and event chunks: a finite sign sample can
   annihilate a nonzero process, so realized draws never decide global degeneracy. */
real scalar cdid_cpt_process_nonzero(struct cdid_cpt_cell colvector cells, real scalar C,
    real scalar gchunk, real scalar echunk)
{
    real scalar i, Md, ne, st, en, es, ee, has_H
    real colvector cols, ii
    real matrix Q, H
    for(i=1;i<=rows(cells);i++) {
        if(any(cells[i].J:!=0)) return(1)
    }
    for(i=1;i<=rows(cells);i++) {
        Md=cells[i].n_distinct; ne=cells[i].n_events
        has_H=(cells[i].pH!=NULL ? (rows(*cells[i].pH)==Md & cols(*cells[i].pH)==Md) : 0)
        for(st=1;st<=Md;st=st+gchunk) {
            en=min((Md,st+gchunk-1)); cols=(st..en)'
            Q=J(C,rows(cols),0)
            for(es=1;es<=ne;es=es+echunk) {
                ee=min((ne,es+echunk-1)); ii=(es..ee)'
                if(has_H) H=(*cells[i].pH)[cells[i].A_d[ii],cols]
                else H=cdid_cpt_indicator(cells[i].zd[cells[i].A_d[ii],.],cells[i].zd[cols,.])
                Q=_cdid_cpt_rowsum_into(Q,H:*cells[i].A_v[ii],cells[i].A_k[ii])
            }
            Q=Q-cells[i].AggK*(cells[i].Rmat[cols,.]')
            if(cells[i].n_coordinates) Q=Q+cells[i].AggX*cells[i].Lmat[.,cols]
            if(any(Q:!=0)) return(1)
        }
    }
    return(0)
}
real scalar _cdid_cpt_preflight_mult(struct cdid_cpt_budget scalar b, struct cdid_cpt_cell colvector cells,
    real scalar N, real scalar C, real scalar B, real scalar br)
{
    real scalar cap, total, i, Md, p, R, Dl, ne, gc, ec, batch, degen, work, one
    if(B>b.max_draws | rows(cells)>b.max_cells)
        _cdid_cpt_fail("Conditional multiplier draws or cells exceed explicit caller limits.",908)
    cap=min((b.max_workspace_cells,floor(b.max_workspace_bytes/8)))
    total=0; batch=0; degen=0; work=0
    for(i=1;i<=rows(cells);i++) {
        Md=cells[i].n_distinct; p=cols(cells[i].zd); R=cells[i].n_terms; Dl=cells[i].n_coordinates
        ne=cells[i].n_events; gc=min((Md,b.grid_chunk)); ec=max((1,min((ne,b.event_chunk))))
        total=_cdid_cpt_charge(Md,p+R+Dl+2,cap,total)
        total=_cdid_cpt_charge(C,R+Dl,cap,total)
        total=_cdid_cpt_charge(ne,3,cap,total)
        total=_cdid_cpt_charge(R+Dl+2,1,cap,total)
        if(cells[i].has_rows) total=_cdid_cpt_charge(N,2+R+Dl,cap,total)
        one=br*(C+R+Dl+Md+3*gc)+Md*gc+ec*(2*br+3)
        batch=max((batch,one))
        one=2*C*gc+ec*gc+3*ec+gc*(R+Dl)
        degen=max((degen,one))
        work=work+64*N*(R+Dl)+16*Md*Md*(p+R+Dl)
        work=work+8*B*(C*(R+Dl)+Md*Md*(p+1)+Md*(R+Dl)+N)
        work=work+8*(C*Md*(R+Dl)+N*Md*(p+1))
    }
    total=_cdid_cpt_charge(B,2,cap,total)
    total=_cdid_cpt_charge(max((batch,degen)),1,cap,total)
    _cdid_cpt_work(b,work)
    return(total)
}
/* Supplied common sign stream: batches in draw order with exact cluster order 1..C, one sign
   matrix shared by every cell and grid point. The caller (later Stata adapter) owns the RNG. */
struct cdid_cpt_mult_state scalar cdid_cpt_mult_start(struct cdid_cpt_cell colvector cells,
    real scalar N, real scalar C, real scalar B, struct cdid_cpt_budget scalar budget)
{
    struct cdid_cpt_mult_state scalar s
    real scalar i, entries, br
    if(rows(cells)<1 | cols(cells)!=1 | !_cdid_cpt_whole(N,1) | !_cdid_cpt_whole(C,1) |
        !_cdid_cpt_whole(B,1) | B>2^31-1)
        _cdid_cpt_fail("Multiplier needs >=1 finalized cell, positive N, C and a positive draw count.",498)
    for(i=1;i<=rows(cells);i++) {
        if(cells[i].M!=N | rows(cells[i].AggK)!=C | cells[i].n_distinct<1 | rows(cells[i].J)!=cells[i].n_distinct)
            _cdid_cpt_fail("Every cell must be finalized on the same full frame and declared clusters.",498)
    }
    br=min((budget.max_sign_batch_rows,floor(budget.max_sign_batch_cells/C),B))
    if(br<1) _cdid_cpt_fail("The sign batch budget admits no complete cluster row.",908)
    entries=_cdid_cpt_preflight_mult(budget,cells,N,C,B,br)
    s.budget=budget; s.N=N; s.C=C; s.B=B; s.n_cells=rows(cells); s.batch_rows=br
    s.cursor=0; s.finished=0; s.CvMb=J(B,1,.)
    return(s)
}
void cdid_cpt_mult_add(struct cdid_cpt_mult_state scalar s, struct cdid_cpt_cell colvector cells,
    real colvector draw_positions, real rowvector cluster_positions, real matrix signs)
{
    real scalar b, i, Md, cap, can_stack, rlo, rhi, nnzQ
    real colvector expected, out, o, dd
    real matrix Q_stack, Js_stack, Js, Emt, S, info
    b=rows(signs)
    if(s.finished | b<1 | cols(signs)!=s.C | rows(draw_positions)!=b | cols(draw_positions)!=1 |
        b>s.B-s.cursor | cols(cluster_positions)!=s.C | rows(cluster_positions)!=1 |
        missing(draw_positions) | missing(cluster_positions) | rows(cells)!=s.n_cells)
        _cdid_cpt_fail("Invalid or completed conditional multiplier sign stream.",498)
    if(b>s.budget.max_sign_batch_rows | b*s.C>s.budget.max_sign_batch_cells)
        _cdid_cpt_fail("Sign batch exceeds explicit row/cell budget.",908)
    expected=((s.cursor+1)..(s.cursor+b))'
    if(any(draw_positions:!=expected) | any(cluster_positions:!=(1..s.C)))
        _cdid_cpt_fail("Supplied signs must preserve exact declared draw and cluster order.",498)
    if(missing(signs) | any((signs:!=1):&(signs:!=-1)))
        _cdid_cpt_fail("Multiplier signs must be exactly -1 or +1.",498)
    out=J(b,1,0)
    Md=cells[1].n_distinct
    cap=min((s.budget.max_workspace_cells,floor(s.budget.max_workspace_bytes/8)))
    can_stack=(s.n_cells>1 & s.budget.grid_chunk>=Md & cells[1].pH!=NULL & s.n_cells*b*Md*2<=cap/2)
    if(can_stack) {
        if(rows(*cells[1].pH)!=Md | cols(*cells[1].pH)!=Md) can_stack=0
        for(i=1;i<=s.n_cells;i++) {
            if(cells[i].n_distinct!=Md | cells[i].pH!=cells[1].pH | s.budget.event_chunk<cells[i].n_events) can_stack=0
        }
    }
    if(can_stack) {
        Q_stack=J(s.n_cells*b,Md,0)
        nnzQ=0
        for(i=1;i<=s.n_cells;i++) {
            if(cells[i].n_events==0) continue
            rlo=(i-1)*b+1; rhi=i*b
            if(cells[i].A_d_uniq) {
                Q_stack[(rlo..rhi),cells[i].A_d]=signs[.,cells[i].A_k]:*(cells[i].A_v')
                nnzQ=nnzQ+b*cells[i].n_events
            }
            else {
                o=cells[i].A_d_ord
                dd=cells[i].A_d[o]
                info=_cdid_ps_cond(dd,1)
                Emt=(signs[.,cells[i].A_k[o]]:*(cells[i].A_v[o]'))'
                S=panelsum(Emt,info)
                Q_stack[(rlo..rhi),dd[info[.,1]]]=S'
                nnzQ=nnzQ+b*rows(info)
            }
        }
        Js_stack=_cdid_cpt_mulH(Q_stack,*cells[1].pH,nnzQ)
        Q_stack=J(0,0,.)
        for(i=1;i<=s.n_cells;i++) {
            rlo=(i-1)*b+1; rhi=i*b
            Js=Js_stack[(rlo..rhi),.]-((signs*cells[i].AggK)*(cells[i].Rmat'))
            if(cells[i].n_coordinates) Js=Js+(signs*cells[i].AggX)*cells[i].Lmat
            Js=Js/s.N
            out=out+(s.N/cells[i].M)*((Js:*Js)*cells[i].mult)
        }
        Js_stack=J(0,0,.)
    }
    else {
        for(i=1;i<=s.n_cells;i++) out=out+_cdid_cpt_batch_cell(signs,cells[i],s.N,s.budget.grid_chunk,s.budget.event_chunk)
    }
    if(missing(out)) _cdid_cpt_fail("Nonfinite conditional multiplier draw.",504)
    s.CvMb[expected]=out
    s.cursor=s.cursor+b
}
/* R stats::quantile type 1 on the sorted draws: j = floor((B p)*(1+fuzz)), x[j + 1{B p > j}],
   with fuzz = 4 eps and the ends clamped to the sample; no interpolation. */
real scalar _cdid_cpt_type1_quantile(real colvector values, real scalar probability)
{
    real colvector ordered
    real scalar B, nppm, fuzz, j, k
    B=rows(values)
    if(B<1 | cols(values)!=1 | missing(values) | probability<0 | probability>1 | probability>=.)
        _cdid_cpt_fail("Type 1 quantile needs finite draws and a probability in [0,1].",498)
    fuzz=4*epsilon(1)
    nppm=B*probability
    j=floor(nppm*(1+fuzz))
    k=j+(nppm>j)
    k=max((1,min((B,k))))
    ordered=sort(values,1)
    return(ordered[k])
}
/* Statistic, cluster-process gate, strict p-value and type 1 critical value. The gate reads
   the actual cluster process, not realized sign outcomes: an annihilated-but-nonzero process
   stays available with its (possibly zero) draws; a globally zero process is unavailable. */
struct cdid_cpt_result scalar cdid_cpt_finish(struct cdid_cpt_mult_state scalar s,
    struct cdid_cpt_cell colvector cells, real scalar alpha)
{
    struct cdid_cpt_result scalar r
    if(s.finished | s.cursor!=s.B | rows(cells)!=s.n_cells)
        _cdid_cpt_fail("The sign stream must contain all draws exactly once before finalization.",498)
    if(alpha<=0 | alpha>=1 | alpha>=.) _cdid_cpt_fail("alpha must lie strictly inside (0,1).",498)
    r.N=s.N; r.C=s.C; r.B=s.B; r.n_cells=s.n_cells; r.alpha=alpha
    r.CvM=cdid_cpt_statistic(cells,s.N)
    if(r.CvM>0 | any(s.CvMb:>0)) r.process_nonzero=1
    else r.process_nonzero=cdid_cpt_process_nonzero(cells,s.C,s.budget.grid_chunk,s.budget.event_chunk)
    if(!r.process_nonzero) {
        r.status="unavailable"
        r.reason="Globally degenerate process: every cluster process row is exactly zero (no estimated process varies across units), so the conditional test is unavailable."
        r.CvMb=J(0,1,0); r.CvMcval=.; r.CvMpval=.
    }
    else {
        r.status="available"; r.reason=""
        r.CvMb=s.CvMb
        r.CvMpval=mean(r.CvMb:>r.CvM)
        r.CvMcval=_cdid_cpt_type1_quantile(r.CvMb,1-alpha)
    }
    if(r.CvM>=.) _cdid_cpt_fail("Nonfinite CvM statistic.",504)
    s.finished=1
    return(r)
}
/* Diagnostics on retained original-unit rows (small inputs): N x n_distinct process rows
   phi_i(u_m) and the scalar link row (h = 1), both oriented by sigma. */
real matrix cdid_cpt_cell_process_rows(struct cdid_cpt_cell scalar cell)
{
    real matrix H, phi
    real scalar Md, has_H
    if(!cell.has_rows) _cdid_cpt_fail("Process rows were not retained for this cell.",498)
    Md=cell.n_distinct
    has_H=(cell.pH!=NULL ? (rows(*cell.pH)==Md & cols(*cell.pH)==Md) : 0)
    if(has_H) H=(*cell.pH)[cell.rows_map,.]
    else H=cdid_cpt_indicator(cell.zd[cell.rows_map,.],cell.zd)
    phi=(H:*cell.rows_A)-cell.rows_kappa*(cell.Rmat')
    if(cell.n_coordinates) phi=phi+cell.rows_xi*cell.Lmat
    return(cell.sigma:*phi)
}
real colvector cdid_cpt_cell_scalar_rows(struct cdid_cpt_cell scalar cell)
{
    real colvector phi
    if(!cell.has_rows) _cdid_cpt_fail("Process rows were not retained for this cell.",498)
    phi=cell.rows_A-cell.rows_kappa*cell.R1
    if(cell.n_coordinates) phi=phi+cell.rows_xi*cell.L1
    return(cell.sigma:*phi)
}
end
