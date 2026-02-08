"""
goldbach_primorial_counts.py

Compute exact Goldbach counts G(N) for primorials N = p_k# for small k:

    G(N) = #{ p <= N/2 : p prime and (N - p) prime }

This is feasible only for relatively small N because it uses a sieve up to N.
Default range k=3..8 is safe (N up to 9,699,690).

Usage:
  python goldbach_primorial_counts.py
  python goldbach_primorial_counts.py --kmin 5 --kmax 8
  python goldbach_primorial_counts.py --k 5 6 7 8
  python goldbach_primorial_counts.py --maxN 20000000

Output: table (k, primorial N, G(N))
"""

from __future__ import annotations

import argparse
import math
from typing import List, Tuple


def sieve_is_prime(n: int) -> bytearray:
    """Return bytearray is_prime[0..n] for n >= 0."""
    if n < 1:
        return bytearray(b"\x00") * (n + 1)
    is_prime = bytearray(b"\x01") * (n + 1)
    is_prime[0:2] = b"\x00\x00"
    # eliminate evens > 2
    for i in range(4, n + 1, 2):
        is_prime[i] = 0
    limit = int(math.isqrt(n))
    for p in range(3, limit + 1, 2):
        if is_prime[p]:
            step = p
            start = p * p
            for x in range(start, n + 1, step):
                is_prime[x] = 0
    is_prime[2] = 1 if n >= 2 else 0
    return is_prime


def primes_up_to(n: int, is_prime: bytearray | None = None) -> List[int]:
    """List primes <= n."""
    if n < 2:
        return []
    if is_prime is None:
        is_prime = sieve_is_prime(n)
    return [2] + [i for i in range(3, n + 1, 2) if is_prime[i]]


def first_k_primes(k: int) -> List[int]:
    """Return [p1, ..., pk]. For small k we can oversieve safely."""
    if k <= 0:
        return []
    # crude upper bound for the kth prime (valid for k>=6): k(log k + log log k)
    if k < 6:
        bound = 15
    else:
        bound = int(k * (math.log(k) + math.log(math.log(k)))) + 10
    # increase bound until we have k primes
    while True:
        is_prime = sieve_is_prime(bound)
        ps = primes_up_to(bound, is_prime)
        if len(ps) >= k:
            return ps[:k]
        bound *= 2


def primorial_from_primes(ps: List[int]) -> int:
    n = 1
    for p in ps:
        n *= p
    return n


def goldbach_count(N: int, is_prime: bytearray) -> int:
    """Compute G(N) = #{p <= N/2 prime : N-p prime}, exact."""
    if N < 4 or (N % 2) == 1:
        return 0
    half = N // 2
    count = 0
    # handle p=2 separately, then odd primes
    if half >= 2 and is_prime[2] and is_prime[N - 2]:
        count += 1
    for p in range(3, half + 1, 2):
        if is_prime[p] and is_prime[N - p]:
            count += 1
    return count


def fmt_int(n: int) -> str:
    return f"{n:,}"


def compute_table(k_values: List[int], maxN: int | None) -> List[Tuple[int, int, int]]:
    """Return list of (k, N, G(N)) for requested k."""
    kmax = max(k_values)
    ps = first_k_primes(kmax)

    rows: List[Tuple[int, int, int]] = []
    for k in k_values:
        N = primorial_from_primes(ps[:k])
        if maxN is not None and N > maxN:
            raise ValueError(f"Primorial p_{k}# = {N} exceeds --maxN={maxN}.")
        # sieve up to N for exact primality on [0..N]
        is_prime = sieve_is_prime(N)
        g = goldbach_count(N, is_prime)
        rows.append((k, N, g))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmin", type=int, default=3, help="Minimum k (default: 3)")
    ap.add_argument("--kmax", type=int, default=8, help="Maximum k (default: 8)")
    ap.add_argument("--k", type=int, nargs="*", default=None, help="Explicit k values (overrides kmin/kmax)")
    ap.add_argument("--maxN", type=int, default=None, help="Abort if primorial exceeds this N")
    args = ap.parse_args()

    if args.k is not None and len(args.k) > 0:
        k_values = sorted(set(args.k))
    else:
        if args.kmin > args.kmax:
            raise SystemExit("--kmin must be <= --kmax")
        k_values = list(range(args.kmin, args.kmax + 1))

    # sanity: exact G(N) requires sieve up to N; warn if huge
    # (user can override; but memory/time will explode)
    rows = compute_table(k_values, args.maxN)

    print("| k | Primorial N = p_k# | Goldbach count G(N) |")
    print("|---:|-------------------:|-------------------:|")
    for k, N, g in rows:
        print(f"| {k} | {fmt_int(N)} | {fmt_int(g)} |")


if __name__ == "__main__":
    main()
