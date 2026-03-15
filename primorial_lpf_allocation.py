"""
primorial_lpf.py
================
Exact least-prime-factor allocation for primorial Goldbach complements.

    N  = p_k#  = product of first k primes
    S  = [p_{k+1}, lambda * p_{k+1})
    C  = { N - p : p in S, p prime }

Every element of C is p_k-rough (Lemma 1.1).  Under the assumption that
no element of C is prime, the LPF partition

    sum_{q > p_k}  C_q  =  pi_S,      C_q = #{p in S : ell(N-p) = q}

must hold exactly.  This script computes the partition, compares it with
the Buchstab-weighted heuristic, and classifies any leftovers.
"""

import math
import sys
import argparse
import time
from collections import defaultdict
from statistics import median

from sympy import primerange, prime, isprime, factorint
from sympy.ntheory.factor_ import pollard_rho


# ============================================================
# Buchstab function
# ============================================================

def buchstab_table(u_max: float = 80.0, h: float = 1e-3):
    """
    Tabulate Buchstab's function omega(u) on [1, u_max].

    Recurrence:
        omega(u)         = 1/u                          for 1 <= u <= 2
        (u omega(u))'    = omega(u-1)                   for u > 2
    Equivalently:
        u omega(u) = 1 + integral_{1}^{u-1} omega(t) dt   for u > 2

    The prefix-integral array is built incrementally so that each
    omega(u) value is computed from a fully-updated integral, avoiding
    the stale-closure issue present in the original code.

    Returns (us, ws): parallel lists of u-values and omega(u) values.
    """
    n = int((u_max - 1.0) / h) + 2
    us = [1.0 + i * h for i in range(n)]
    ws = [0.0] * n
    # prefix[i] = integral_{1}^{us[i]} omega(t) dt  (trapezoidal)
    prefix = [0.0] * n

    # Phase 1: u in [1, 2]
    for i, u in enumerate(us):
        if u > 2.0 + 0.5 * h:
            break
        ws[i] = 1.0 / u
        if i > 0:
            prefix[i] = prefix[i - 1] + 0.5 * h * (ws[i - 1] + ws[i])

    # Phase 2: u > 2  — each step uses the prefix array already written
    for i in range(1, n):
        u = us[i]
        if u <= 2.0 + 0.5 * h:
            continue
        # Interpolate integral up to u - 1
        x = u - 1.0
        if x <= us[0]:
            integ = 0.0
        else:
            j = min(int((x - 1.0) / h), n - 2)
            t = (x - us[j]) / h
            integ = prefix[j] * (1.0 - t) + prefix[j + 1] * t
        ws[i] = (1.0 + integ) / u
        prefix[i] = prefix[i - 1] + 0.5 * h * (ws[i - 1] + ws[i])

    return us, ws


def buchstab_interp(us: list, ws: list, u: float) -> float:
    """Linear interpolation into a Buchstab table."""
    EULER_GAMMA = 0.5772156649015329
    if u <= us[0]:
        return 1.0 / max(u, 1e-12)
    if u >= us[-1]:
        return math.exp(-EULER_GAMMA)
    h = us[1] - us[0]
    j = min(int((u - us[0]) / h), len(us) - 2)
    t = (u - us[j]) / h
    return ws[j] * (1.0 - t) + ws[j + 1] * t


# ============================================================
# Primorial helpers
# ============================================================

def prime_list(k: int) -> list:
    """Return the first k primes."""
    return [prime(i) for i in range(1, k + 1)]


def primorial(primes: list) -> int:
    out = 1
    for p in primes:
        out *= p
    return out


def theta(primes: list) -> float:
    """log-primorial: sum of log p for p in primes."""
    return sum(math.log(p) for p in primes)


def first_in_AP(lo: int, residue: int, mod: int) -> int:
    """Smallest integer >= lo that is congruent to residue (mod mod)."""
    delta = (residue - lo % mod) % mod
    return lo + delta


def _least_prime_factor_recursive(n: int, rho_rounds: int = 8) -> int:
    """Return the least prime factor of n using recursive Pollard-rho splitting."""
    if n % 2 == 0:
        return 2
    if isprime(n):
        return n

    seeds = [2, 3, 5, 7, 11, 13, 17, 19]
    factor = None
    for seed in seeds[:rho_rounds]:
        try:
            factor = pollard_rho(n, seed=seed)
        except Exception:
            factor = None
        if factor not in (None, 1, n):
            break
    if factor in (None, 1, n):
        # Fall back to sympy's factorint if Pollard-rho stalls.
        fac = factorint(n)
        return min(fac.keys())

    a = _least_prime_factor_recursive(factor, rho_rounds=rho_rounds)
    b = _least_prime_factor_recursive(n // factor, rho_rounds=rho_rounds)
    return min(a, b)


def least_prime_factor_above(m: int, threshold: int, trial_limit: int = 10**6) -> int:
    """
    Return the least prime factor of m that exceeds threshold.

    Fast path: trial-divide only by *primes* between threshold and trial_limit.
    Slow path: recursively split the remaining cofactor with Pollard-rho,
    returning only the least prime factor instead of a full factorization.
    """
    if m <= 1:
        raise ValueError(f"m must be > 1, got {m}")

    sqrt_m = int(math.isqrt(m))
    if sqrt_m <= threshold:
        return m if isprime(m) else _least_prime_factor_recursive(m)

    limit = min(trial_limit, sqrt_m)
    start = max(2, threshold + 1)
    for f in primerange(start, limit + 1):
        if m % f == 0:
            return f

    if limit >= sqrt_m:
        return m if isprime(m) else _least_prime_factor_recursive(m)

    return _least_prime_factor_recursive(m)


# ============================================================
# Core LPF allocation
# ============================================================


def batched_sparse_lpf_scan(
    N: int,
    active_ps: set,
    S_lo: int,
    S_hi: int,
    q_lo: int,
    q_hi: int,
    hist_exact: dict,
    assigned: dict,
    progress: bool = False,
    progress_label: str = "",
    progress_every: int = 5000,
):
    """
    Assign LPFs in batch for complements N-p, restricted to the currently active
    prime positions p in S. For each q, we only scan the congruence class
    p ≡ N (mod q) inside [S_lo, S_hi).

    Returns (resolved_count, q_processed).
    """
    if not active_ps or q_hi < q_lo:
        return 0, 0

    resolved = 0
    q_processed = 0
    q_lo = max(2, q_lo)
    active = active_ps

    for q in primerange(q_lo, q_hi + 1):
        q_processed += 1
        if not active:
            break
        r = int(N % q)
        p = first_in_AP(S_lo, r, q)
        while p < S_hi:
            if p in active:
                assigned[p] = q
                hist_exact[q] += 1
                active.remove(p)
                resolved += 1
            p += q
        if progress and q_processed % max(1, progress_every // 2) == 0:
            print(
                f"{progress_label} q={q:,} processed={q_processed:,} resolved={resolved:,} active={len(active):,}",
                file=sys.stderr,
                flush=True,
            )
    return resolved, q_processed


def exact_lpf_allocation(
    k: int,
    lam: float = 2.0,
    q_max: int = None,
    verify_assignments: bool = False,
    buchstab_cache=None,
    full_factorization: bool = False,
    skip_leftover_lpf=False,
    progress: bool = False,
    progress_every: int = 5000,
) -> dict:
    """
    Compute the exact LPF partition for N = p_k# over S = [p_{k+1}, lam*p_{k+1}).

    Optimization for larger k:
      Phase I  : exact q-scan only up to q_small = min(q_max, S_hi - 1)
      Phase II : for any remaining p, compute LPF(m) directly for m = N-p;
                 if LPF(m) <= q_max, assign it immediately.
    For q > S_hi, every occupied residue class has capacity at most 1, so the
    direct LPF route is much cheaper than scanning every prime q up to q_max.
    """
    small_primes = prime_list(k)
    pk  = small_primes[-1]
    pk1 = prime(k + 1)

    N    = primorial(small_primes)
    logN = theta(small_primes)

    S_lo = pk1
    S_hi = int(lam * pk1)

    primes_S     = list(primerange(S_lo, S_hi))
    primes_S_set = set(primes_S)
    pi_S         = len(primes_S)

    if q_max is None:
        q_max = max(200 * pk, S_hi * 100)

    remaining   = set(primes_S)
    hist_exact  = defaultdict(int)
    assigned    = {}

    q_small = min(q_max, max(pk + 1, S_hi - 1))

    # --- Phase I: exact q-scan only in the short-modulus regime -----------
    phase1_count = 0
    t0 = time.time()
    if progress:
        print(f"[run k={k} lam={lam}] start: pi_S={pi_S}, q_max={q_max}, q_small={q_small}", file=sys.stderr, flush=True)

    for q in primerange(pk + 1, q_small + 1):
        phase1_count += 1
        if not remaining:
            break
        r     = int(N % q)
        start = first_in_AP(S_lo, r, q)
        p     = start
        while p < S_hi:
            if p in remaining:
                if verify_assignments:
                    assert (N - p) % q == 0, (
                        f"Assignment error: q={q} does not divide N-p "
                        f"for p={p}"
                    )
                assigned[p]  = q
                remaining.remove(p)
                hist_exact[q] += 1
            p += q
        if progress and phase1_count % progress_every == 0:
            print(
                f"[run k={k} lam={lam}] phase1 q={q:,} processed={phase1_count:,} "
                f"assigned={len(assigned):,} remaining={len(remaining):,} elapsed={time.time()-t0:.1f}s",
                file=sys.stderr,
            )

    # --- Phase II: batched sparse LPF scan for larger q -------------------
    # Instead of factoring each N-p independently, scan only the relevant
    # congruence class p ≡ N (mod q) through the still-active prime positions.
    phase2_trial_limit = min(q_max, 1_000_000)
    direct_assigned = 0
    phase2_q_processed = 0

    if (not skip_leftover_lpf) and q_max > q_small and remaining:
        if progress:
            print(
                f"[run k={k} lam={lam}] phase2 batched LPF on {len(remaining):,} leftovers "
                f"for q in [{q_small + 1:,}, {phase2_trial_limit:,}]",
                file=sys.stderr,
                flush=True,
            )
        direct_assigned, phase2_q_processed = batched_sparse_lpf_scan(
            N=N,
            active_ps=remaining,
            S_lo=S_lo,
            S_hi=S_hi,
            q_lo=q_small + 1,
            q_hi=phase2_trial_limit,
            hist_exact=hist_exact,
            assigned=assigned,
            progress=progress,
            progress_label=f"[run k={k} lam={lam}] phase2",
            progress_every=progress_every,
        )
        if progress:
            print(
                f"[run k={k} lam={lam}] phase2 batch resolved={direct_assigned:,} "
                f"tail={len(remaining):,} q_processed={phase2_q_processed:,}",
                file=sys.stderr,
                flush=True,
            )

    # --- Buchstab heuristic model ----------------------------------------
    if buchstab_cache is None:
        us, ws = buchstab_table()
    else:
        us, ws = buchstab_cache

    def model_Cq(q: int) -> float:
        arg   = (logN - math.log(q)) / math.log(q)
        omega = buchstab_interp(us, ws, arg)
        return pi_S * math.log(pk) / (q * math.log(q)) * omega

    hist_model = {q: model_Cq(q) for q in hist_exact}

    # --- cumulative coverage rows ----------------------------------------
    band_cutoffs = [2, 3, 5, 10, 20, 50, 100, 200, 500, 1000]
    exact_sorted = sorted(hist_exact.items())

    def cov_exact(Q: int) -> int:
        return sum(c for q, c in exact_sorted if q <= Q)

    # Avoid expensive full prime scan for the model when Q is large.
    model_prefix = []
    running_model = 0.0
    for q in primerange(pk + 1, min(q_max, max(mult * pk for mult in band_cutoffs)) + 1):
        running_model += model_Cq(q)
        model_prefix.append((q, running_model))

    def cov_model(Q: int) -> float:
        total = 0.0
        for q, s in model_prefix:
            if q > Q:
                break
            total = s
        return total

    coverage_rows = []
    for mult in band_cutoffs:
        Q = mult * pk
        if Q > q_max:
            continue
        ce = cov_exact(Q)
        cm = cov_model(Q)
        coverage_rows.append({
            "Q":             Q,
            "exact_cov":     ce,
            "exact_frac":    ce / pi_S if pi_S else 0.0,
            "model_cov":     cm,
            "model_frac":    cm / pi_S if pi_S else 0.0,
            "uncovered":     pi_S - ce,
        })

    # --- leftover classification -----------------------------------------
    leftover_records = []
    rem_sorted = sorted(remaining)
    if progress and rem_sorted:
        mode = "primality-only" if (skip_leftover_lpf or len(rem_sorted) > 64) else "with LPF"
        print(
            f"[run k={k} lam={lam}] classifying {len(rem_sorted):,} final leftovers ({mode})",
            file=sys.stderr,
            flush=True,
        )
    classify_lpf = (not skip_leftover_lpf) and len(rem_sorted) <= 64
    for idx, p in enumerate(rem_sorted, 1):
        m   = N - p
        rec = {"p": p, "m": m, "status": None, "lpf": None, "factorization": None}
        if isprime(m):
            rec["status"] = "prime"
        else:
            if not classify_lpf:
                rec["status"] = "composite_unfactored"
            else:
                rec["status"] = "composite"
                try:
                    lpf = least_prime_factor_above(m, threshold=pk, trial_limit=phase2_trial_limit)
                    rec["lpf"] = lpf
                    if full_factorization:
                        rec["factorization"] = factorint(m)
                except Exception as e:
                    rec["lpf"] = None
                    rec["notes"] = str(e)
        leftover_records.append(rec)
        if progress and idx % max(1, progress_every // 10) == 0:
            print(
                f"[run k={k} lam={lam}] leftover classify={idx:,}/{len(rem_sorted):,} elapsed={time.time()-t0:.1f}s",
                file=sys.stderr,
                flush=True,
            )

    if progress:
        print(
            f"[run k={k} lam={lam}] done: phase1_assign={len(assigned)-direct_assigned:,} "
            f"phase2_assign={direct_assigned:,} leftovers={len(remaining):,} total_elapsed={time.time()-t0:.1f}s",
            file=sys.stderr,
        )

    return {
        "k":               k,
        "pk":              pk,
        "pk1":             pk1,
        "N":               N,
        "logN":            logN,
        "S_lo":            S_lo,
        "S_hi":            S_hi,
        "lam":             lam,
        "q_max":           q_max,
        "q_small":         q_small,
        "primes_S":        primes_S,
        "primes_S_set":    primes_S_set,
        "pi_S":            pi_S,
        "assigned":        assigned,
        "hist_exact":      dict(hist_exact),
        "hist_model":      hist_model,
        "coverage_rows":   coverage_rows,
        "leftover_records": leftover_records,
        "remaining_after_qmax": sorted(remaining),
    }


# ============================================================
# Diagnostic summaries
# ============================================================

def lpf_band_summary(data: dict) -> list:
    """Count assigned LPFs by multiplicative band above p_k."""
    pk         = data["pk"]
    hist_exact = data["hist_exact"]
    bands = [
        ("(pk, 2pk]",       pk,         2   * pk),
        ("(2pk, 5pk]",      2   * pk,   5   * pk),
        ("(5pk, 10pk]",     5   * pk,   10  * pk),
        ("(10pk, 20pk]",    10  * pk,   20  * pk),
        ("(20pk, 50pk]",    20  * pk,   50  * pk),
        ("(50pk, 100pk]",   50  * pk,   100 * pk),
        ("(100pk, inf)",    100 * pk,   math.inf),
    ]
    return [
        (label, sum(c for q, c in hist_exact.items() if lo < q <= hi))
        for label, lo, hi in bands
    ]


def leftover_summary(data: dict) -> dict:
    """Aggregate statistics on leftovers (unassigned after q_max)."""
    leftovers    = data["leftover_records"]
    pk           = data["pk"]
    prime_count  = sum(1 for r in leftovers if r["status"] == "prime")
    comp_records = [r for r in leftovers if r["status"] in ("composite", "composite_unfactored")]
    comp_records_with_lpf = [r for r in comp_records if r.get("lpf") is not None]

    comp_lpfs          = [r["lpf"] for r in comp_records_with_lpf]
    logexp_lpf         = []
    lpf_over_pk        = []
    log_lpf_over_log_pk = []

    exp_bands = {"<=m^1/4": 0, "(m^1/4,m^1/3]": 0, "(m^1/3,m^0.4]": 0,
                 "(m^0.4,m^0.49]": 0, ">m^0.49": 0}
    pk_bands  = {"(pk,2pk]": 0, "(2pk,5pk]": 0, "(5pk,10pk]": 0,
                 "(10pk,1e2pk]": 0, "(1e2pk,1e3pk]": 0,
                 "(1e3pk,1e4pk]": 0, ">1e4pk": 0}

    for r in comp_records_with_lpf:
        m   = r["m"]
        lpf = r["lpf"]

        alpha = math.log(lpf) / math.log(m)
        logexp_lpf.append(alpha)

        if   alpha <= 0.25:          exp_bands["<=m^1/4"]       += 1
        elif alpha <= 1.0 / 3.0:     exp_bands["(m^1/4,m^1/3]"] += 1
        elif alpha <= 0.4:           exp_bands["(m^1/3,m^0.4]"] += 1
        elif alpha <= 0.49:          exp_bands["(m^0.4,m^0.49]"]+= 1
        else:                        exp_bands[">m^0.49"]        += 1

        ratio = lpf / pk
        gamma = math.log(lpf) / math.log(pk)
        lpf_over_pk.append(ratio)
        log_lpf_over_log_pk.append(gamma)

        if   lpf <= 2   * pk: pk_bands["(pk,2pk]"]       += 1
        elif lpf <= 5   * pk: pk_bands["(2pk,5pk]"]       += 1
        elif lpf <= 10  * pk: pk_bands["(5pk,10pk]"]      += 1
        elif lpf <= 1e2 * pk: pk_bands["(10pk,1e2pk]"]    += 1
        elif lpf <= 1e3 * pk: pk_bands["(1e2pk,1e3pk]"]   += 1
        elif lpf <= 1e4 * pk: pk_bands["(1e3pk,1e4pk]"]   += 1
        else:                 pk_bands[">1e4pk"]           += 1

    def _stats(xs):
        if not xs:
            return None, None, None
        return median(xs), min(xs), max(xs)

    med_lpf, min_lpf, max_lpf = _stats(comp_lpfs)
    med_exp, min_exp, max_exp = _stats(logexp_lpf)
    med_rat, min_rat, max_rat = _stats(lpf_over_pk)
    med_gam, min_gam, max_gam = _stats(log_lpf_over_log_pk)

    return {
        "leftover_total":           len(leftovers),
        "leftover_primes":          prime_count,
        "leftover_composites":      len(comp_records),
        "leftover_composites_with_lpf": len(comp_records_with_lpf),
        "leftover_prime_frac":      prime_count / len(leftovers) if leftovers else 0.0,

        "median_comp_lpf":          med_lpf,
        "min_comp_lpf":             min_lpf,
        "max_comp_lpf":             max_lpf,

        "median_logexp_lpf":        med_exp,
        "min_logexp_lpf":           min_exp,
        "max_logexp_lpf":           max_exp,

        "median_lpf_over_pk":       med_rat,
        "min_lpf_over_pk":          min_rat,
        "max_lpf_over_pk":          max_rat,

        "median_log_lpf_over_log_pk": med_gam,
        "min_log_lpf_over_log_pk":    min_gam,
        "max_log_lpf_over_log_pk":    max_gam,

        "comp_exp_band_counts":     exp_bands,
        "comp_pk_band_counts":      pk_bands,
    }


def residue_class_capacity(data: dict, q: int) -> int:
    """
    Exact residue-class capacity
        A_q = #{p in S : p ≡ N (mod q), p prime}.

    This is the total number of prime positions in the window S that the
    congruence class attached to q could possibly hit.
    """
    N = data["N"]
    S_lo = data["S_lo"]
    S_hi = data["S_hi"]
    primes_S_set = data["primes_S_set"]

    r = int(N % q)
    start = first_in_AP(S_lo, r, q)
    count = 0
    p = start
    while p < S_hi:
        if p in primes_S_set:
            count += 1
        p += q
    return count


def lpf_capacity_summary(data: dict) -> dict:
    """
    Capacity proxy for composite leftovers based on exact residue-class sizes.

    For each composite leftover with LPF q, let
        L_q = number of composite leftovers with LPF q,
        A_q = #{p in S : p ≡ N (mod q), p prime}.

    The normalized capacity usage is
        sum_q L_q / A_q.
    Classes with A_q = 1 are effectively unit-capacity classes.
    """
    comp_records = [r for r in data["leftover_records"]
                    if r["status"] == "composite" and r.get("lpf") is not None]

    Lq = defaultdict(int)
    for r in comp_records:
        Lq[r["lpf"]] += 1

    Aq = {}
    rows = []
    capacity_used = 0.0
    capacity_waste = 0
    unit_classes = 0
    two_classes = 0
    le3_classes = 0

    for q in sorted(Lq):
        L = Lq[q]
        A = residue_class_capacity(data, q)
        Aq[q] = A
        occ = (L / A) if A > 0 else float('inf')
        capacity_used += occ if math.isfinite(occ) else 0.0
        capacity_waste += max(0, A - L) if A > 0 else 0
        if A == 1:
            unit_classes += 1
        if A == 2:
            two_classes += 1
        if A <= 3:
            le3_classes += 1
        rows.append((q, L, A, occ))

    num_classes = len(Lq)
    total_composites = sum(Lq.values())
    reuse_ratio = (total_composites / num_classes) if num_classes else 0.0
    max_multiplicity = max(Lq.values()) if Lq else 0
    singleton_classes = sum(1 for v in Lq.values() if v == 1)

    return {
        "num_comp_lpf_classes": num_classes,
        "capacity_used": capacity_used,
        "capacity_waste": capacity_waste,
        "unit_classes": unit_classes,
        "two_classes": two_classes,
        "le3_classes": le3_classes,
        "reuse_ratio": reuse_ratio,
        "max_multiplicity": max_multiplicity,
        "singleton_classes": singleton_classes,
        "rows": rows,
    }


# ============================================================
# Printing helpers
# ============================================================

def print_run_summary(data: dict, top_n: int = 20, show_leftovers: bool = True):
    pk       = data["pk"]
    pk1      = data["pk1"]
    pi_S     = data["pi_S"]
    q_max    = data["q_max"]
    hist_e   = data["hist_exact"]
    hist_m   = data["hist_model"]

    print(f"  N = p_{data['k']}#  (p_k = {pk},  p_{{k+1}} = {pk1})")
    print(f"  S = [{data['S_lo']}, {data['S_hi']})   pi_S = {pi_S}   lambda = {data['lam']}")
    print(f"  q_max = {q_max}")

    exact_covered = sum(hist_e.values())
    print()
    print("  === LPF ALLOCATION ===")
    print(f"  Assigned by q <= q_max : {exact_covered}")
    print(f"  Unassigned (leftovers) : {pi_S - exact_covered}")

    print()
    print(f"  {'q':>10} {'exact':>8} {'model':>12} {'ratio':>10}")
    print("  " + "-" * 44)
    for q, c in sorted(hist_e.items(), key=lambda x: (-x[1], x[0]))[:top_n]:
        m   = hist_m.get(q, 0.0)
        rat = f"{c/m:.3f}" if m > 0 else "inf"
        print(f"  {q:10d} {c:8d} {m:12.4f} {rat:>10}")

    print()
    print("  === CUMULATIVE COVERAGE ===")
    print(f"  {'Q':>8} {'cov':>6} {'cov/pi_S':>10} {'model':>10} {'mdl/pi_S':>10} {'left':>6}")
    print("  " + "-" * 58)
    for row in data["coverage_rows"]:
        print(
            f"  {row['Q']:8d} {row['exact_cov']:6d} "
            f"{row['exact_frac']:10.4f} "
            f"{row['model_cov']:10.4f} "
            f"{row['model_frac']:10.4f} "
            f"{row['uncovered']:6d}"
        )

    print()
    print("  === LPF BANDS (all assigned) ===")
    for label, count in lpf_band_summary(data):
        print(f"  {label:>18}  {count:6d}")

    ls = leftover_summary(data)
    print()
    print("  === LEFTOVER SUMMARY ===")
    print(f"  total      = {ls['leftover_total']}")
    print(f"  primes     = {ls['leftover_primes']}  ({ls['leftover_prime_frac']:.4f})")
    print(f"  composites = {ls['leftover_composites']}")
    if ls["leftover_composites"] != ls.get("leftover_composites_with_lpf", ls["leftover_composites"]):
        print(f"  composites with LPF = {ls['leftover_composites_with_lpf']}")

    if ls["median_lpf_over_pk"] is not None:
        print()
        print("  pk-relative LPF stats (composite leftovers):")
        print(f"    median LPF/p_k           = {ls['median_lpf_over_pk']:.4g}")
        print(f"    median log(LPF)/log(p_k) = {ls['median_log_lpf_over_log_pk']:.4f}")
        print()
        print("  pk-relative bands:")
        for label, cnt in ls["comp_pk_band_counts"].items():
            print(f"    {label:>18}  {cnt:6d}")
        print()
        print("  Exponent bands (log LPF / log m):")
        for label, cnt in ls["comp_exp_band_counts"].items():
            print(f"    {label:>18}  {cnt:6d}")

        cap = lpf_capacity_summary(data)
        print()
        print("  === LPF CAPACITY PROXY (composite leftovers) ===")
        print(f"  composite LPF classes    = {cap['num_comp_lpf_classes']}")
        print(f"  normalized capacity used = {cap['capacity_used']:.6f}")
        print(f"  capacity waste           = {cap['capacity_waste']}")
        print(f"  unit-capacity classes    = {cap['unit_classes']}")
        print(f"  two-capacity classes     = {cap['two_classes']}")
        print(f"  classes with A_q <= 3    = {cap['le3_classes']}")
        print(f"  reuse ratio              = {cap['reuse_ratio']:.6f}")
        print(f"  max class multiplicity   = {cap['max_multiplicity']}")
        print(f"  singleton LPF classes    = {cap['singleton_classes']}")

        if cap['rows']:
            print()
            print(f"  {'q':>12} {'L_q':>8} {'A_q':>8} {'L_q/A_q':>12}")
            print("  " + "-" * 44)
            for q, L, A, occ in cap['rows'][:top_n]:
                occ_s = f"{occ:.6f}" if math.isfinite(occ) else "inf"
                print(f"  {q:12d} {L:8d} {A:8d} {occ_s:>12}")

    if show_leftovers and data["leftover_records"]:
        print()
        print("  === LEFTOVER DETAIL ===")
        for rec in data["leftover_records"]:
            p, m = rec["p"], rec["m"]
            if rec["status"] == "prime":
                print(f"    p={p}  N-p={m}  PRIME")
            else:
                lpf = rec["lpf"]
                if lpf is not None:
                    alpha = math.log(lpf) / math.log(m)
                    gamma = math.log(lpf) / math.log(pk)
                    ratio = lpf / pk
                    print(
                        f"    p={p}  N-p={m}  COMPOSITE"
                        f"  lpf={lpf}"
                        f"  LPF/p_k={ratio:.3g}"
                        f"  log(LPF)/log(p_k)={gamma:.4f}"
                        f"  log(LPF)/log(m)={alpha:.4f}"
                        f"  fac={rec['factorization']}"
                    )
                else:
                    print(f"    p={p}  N-p={m}  COMPOSITE  lpf=UNKNOWN  {rec.get('notes','')}")


# ============================================================
# CSV output
# ============================================================

CSV_HEADER = (
    "k,pk,pk1,lambda,pi_S,q_max,"
    "exact_covered,leftovers,leftover_primes,leftover_composites,leftover_prime_frac,"
    "median_comp_lpf,min_comp_lpf,max_comp_lpf,"
    "median_lpf_over_pk,min_lpf_over_pk,max_lpf_over_pk,"
    "median_log_lpf_over_log_pk,min_log_lpf_over_log_pk,max_log_lpf_over_log_pk,"
    "median_logexp_lpf,min_logexp_lpf,max_logexp_lpf,"
    "cap_num_comp_classes,cap_used,cap_waste,cap_unit_classes,cap_two_classes,cap_le3_classes,"
    "cap_reuse_ratio,cap_max_multiplicity,cap_singleton_classes,"
    "pk_band_pk_2pk,pk_band_2pk_5pk,pk_band_5pk_10pk,"
    "pk_band_10pk_1e2,pk_band_1e2_1e3,pk_band_1e3_1e4,pk_band_gt1e4,"
    "exp_band_le14,exp_band_14_13,exp_band_13_04,exp_band_04_049,exp_band_gt049,"
    "alloc_band_pk_2pk,alloc_band_2pk_5pk,alloc_band_5pk_10pk,"
    "alloc_band_10pk_20pk,alloc_band_20pk_50pk,alloc_band_50pk_100pk,alloc_band_gt100pk"
)


def csv_row(data: dict) -> str:
    ls    = leftover_summary(data)
    cap   = lpf_capacity_summary(data)
    bands = dict(lpf_band_summary(data))
    cb_e  = ls["comp_exp_band_counts"]
    cb_p  = ls["comp_pk_band_counts"]
    exact_covered = sum(data["hist_exact"].values())

    _f  = lambda x: "" if x is None else str(x)
    _ff = lambda x: "" if x is None else f"{x:.6f}"

    return ",".join([
        str(data["k"]),
        str(data["pk"]),
        str(data["pk1"]),
        str(data["lam"]),
        str(data["pi_S"]),
        str(data["q_max"]),
        str(exact_covered),
        str(ls["leftover_total"]),
        str(ls["leftover_primes"]),
        str(ls["leftover_composites"]),
        f"{ls['leftover_prime_frac']:.6f}",
        _f(ls["median_comp_lpf"]),
        _f(ls["min_comp_lpf"]),
        _f(ls["max_comp_lpf"]),
        _ff(ls["median_lpf_over_pk"]),
        _ff(ls["min_lpf_over_pk"]),
        _ff(ls["max_lpf_over_pk"]),
        _ff(ls["median_log_lpf_over_log_pk"]),
        _ff(ls["min_log_lpf_over_log_pk"]),
        _ff(ls["max_log_lpf_over_log_pk"]),
        _ff(ls["median_logexp_lpf"]),
        _ff(ls["min_logexp_lpf"]),
        _ff(ls["max_logexp_lpf"]),
        str(cap["num_comp_lpf_classes"]),
        f"{cap['capacity_used']:.6f}",
        str(cap["capacity_waste"]),
        str(cap["unit_classes"]),
        str(cap["two_classes"]),
        str(cap["le3_classes"]),
        f"{cap['reuse_ratio']:.6f}",
        str(cap["max_multiplicity"]),
        str(cap["singleton_classes"]),
        str(cb_p["(pk,2pk]"]),
        str(cb_p["(2pk,5pk]"]),
        str(cb_p["(5pk,10pk]"]),
        str(cb_p["(10pk,1e2pk]"]),
        str(cb_p["(1e2pk,1e3pk]"]),
        str(cb_p["(1e3pk,1e4pk]"]),
        str(cb_p[">1e4pk"]),
        str(cb_e["<=m^1/4"]),
        str(cb_e["(m^1/4,m^1/3]"]),
        str(cb_e["(m^1/3,m^0.4]"]),
        str(cb_e["(m^0.4,m^0.49]"]),
        str(cb_e[">m^0.49"]),
        str(bands["(pk, 2pk]"]),
        str(bands["(2pk, 5pk]"]),
        str(bands["(5pk, 10pk]"]),
        str(bands["(10pk, 20pk]"]),
        str(bands["(20pk, 50pk]"]),
        str(bands["(50pk, 100pk]"]),
        str(bands["(100pk, inf)"]),
    ])


# ============================================================
# Batch driver
# ============================================================

def run_batch(
    k_values,
    lam_values,
    q_max=None,
    detailed=False,
    csv_out=None,
    verify=False,
    full_factorization=False,
    skip_leftover_lpf=False,
    progress=False,
    progress_every=100000,
):
    """
    Run exact_lpf_allocation for every (k, lambda) pair.

    Parameters
    ----------
    k_values    : iterable of k values
    lam_values  : iterable of lambda values
    q_max       : upper bound on q search (None = auto)
    detailed    : if True, print per-run diagnostics
    verify      : if True, assert q | (N-p) for all assignments
    csv_out     : if given, write CSV to this file path (in addition to stdout)
    """
    buchstab_cache = buchstab_table()

    csv_lines = [CSV_HEADER]
    print(CSV_HEADER)

    for lam in lam_values:
        for k in k_values:
            if progress:
                print(f"[batch] starting k={k}, lambda={lam}", file=sys.stderr, flush=True)
            data = exact_lpf_allocation(
                k=k,
                lam=lam,
                q_max=q_max,
                verify_assignments=verify,
                buchstab_cache=buchstab_cache,
                full_factorization=full_factorization,
                skip_leftover_lpf=skip_leftover_lpf,
                progress=progress,
                progress_every=progress_every,
            )
            row = csv_row(data)
            csv_lines.append(row)
            print(row)

            if detailed:
                print()
                print("=" * 70)
                print_run_summary(data, top_n=20, show_leftovers=True)
                print("=" * 70)
                print()

    if csv_out:
        with open(csv_out, "w") as fh:
            fh.write("\n".join(csv_lines) + "\n")
        print(f"\nCSV written to {csv_out}", file=sys.stderr, flush=True)


# ============================================================
# CLI entry point
# ============================================================

def _parse_args():
    p = argparse.ArgumentParser(
        description="LPF allocation diagnostics for primorial Goldbach",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--k-min",    type=int,   default=33)
    p.add_argument("--k-max",    type=int,   default=34)
    p.add_argument("--lam",      type=float, nargs="+", default=[5.0])
    p.add_argument("--q-max",    type=int,   default=2_000_000)
    p.add_argument("--detailed", action="store_true",
                   help="Print full per-run diagnostics")
    p.add_argument("--verify",   action="store_true",
                   help="Assert q | (N-p) for every assignment (slow for large k)")
    p.add_argument("--csv",      type=str,   default=None,
                   help="Write CSV summary to this file")
    p.add_argument("--full-factorization", action="store_true",
                   help="Factor leftover composites completely (much slower for k >= 34)")
    p.add_argument("--skip-leftover-lpf", action="store_true",
                   help="Only test leftover primality; skip LPF extraction for composite leftovers")
    p.add_argument("--progress", action="store_true",
                   help="Print progress updates to stderr during long runs")
    p.add_argument("--progress-every", type=int, default=5000,
                   help="Progress cadence for phase-I prime scan")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    run_batch(
        k_values    = range(args.k_min, args.k_max + 1),
        lam_values  = args.lam,
        q_max       = args.q_max,
        detailed    = args.detailed,
        verify      = args.verify,
        csv_out     = args.csv,
        full_factorization = args.full_factorization,
        skip_leftover_lpf = args.skip_leftover_lpf,
        progress = args.progress,
        progress_every = args.progress_every,
    )
