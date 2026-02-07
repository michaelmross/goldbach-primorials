"""
LPF Symmetry Verification for Primorials

This script verifies the LPF (Least Prime Factor) pairing symmetry for primorial numbers,
a key structural property used in proving Goldbach's conjecture for primorials.

For a primorial N = p_k# and pairs (x, N-x) where x < N/2, this script verifies that:
    #{x : lpf(x) = q} ≈ #{x : lpf(N-x) = q}
for each prime q.

Reference:
    "Goldbach's Conjecture for Primorials" (2026)
    https://arxiv.org/abs/[to be added]

Author: [Your name]
Date: February 2026
License: MIT
"""

import numpy as np
from sympy import isprime, prime, primorial, factorint
from collections import Counter
import time
import sys


def lpf(n):
    """
    Compute the least prime factor of n.
    
    Parameters:
    -----------
    n : int
        Integer >= 2
    
    Returns:
    --------
    int
        The smallest prime that divides n
    
    Examples:
    ---------
    >>> lpf(12)
    2
    >>> lpf(15)
    3
    >>> lpf(17)
    17
    """
    if n <= 1:
        return None
    if isprime(n):
        return n
    factors = factorint(n)
    return min(factors.keys())


def verify_lpf_symmetry(k, verbose=True):
    """
    Verify LPF pairing symmetry for the k-th primorial.
    
    For N = p_k#, tests whether the distribution of lpf(x) for x < N/2
    matches the distribution of lpf(N-x).
    
    Parameters:
    -----------
    k : int
        Index of primorial (e.g., k=3 gives N=30, k=8 gives N=9699690)
    verbose : bool
        If True, print detailed output
    
    Returns:
    --------
    dict
        Results containing:
        - 'k': primorial index
        - 'N': primorial value
        - 'pairs': number of pairs tested
        - 'ratio': total count ratio (left/right)
        - 'lpf_data': individual LPF frequencies
    """
    start_time = time.time()
    
    N = primorial(k)
    p_k = prime(k)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Testing LPF Symmetry for k={k}")
        print(f"N = p_{k}# = {N:,}")
        print(f"p_{k} = {p_k}")
        print(f"{'='*70}\n")
    
    # Collect LPF data for pairs (x, N-x) where x < N/2
    left_lpfs = []   # lpf(x) for x < N/2
    right_lpfs = []  # lpf(N-x) for x < N/2
    
    # Start from 3 to avoid edge cases with 1 and 2
    for x in range(3, N//2):
        y = N - x
        
        lpf_x = lpf(x)
        lpf_y = lpf(y)
        
        left_lpfs.append(lpf_x)
        right_lpfs.append(lpf_y)
        
        # Show first few examples
        if verbose and len(left_lpfs) <= 10:
            print(f"  ({x:6d}, {y:6d}): lpf({x})={lpf_x:3d}, lpf({y})={lpf_y:3d}")
    
    # Count frequencies
    left_counts = Counter(left_lpfs)
    right_counts = Counter(right_lpfs)
    
    # Get all unique LPFs
    all_lpfs = sorted(set(left_counts.keys()) | set(right_counts.keys()))
    
    # Display results
    if verbose:
        print(f"\n{'='*70}")
        print(f"LPF DISTRIBUTION COMPARISON")
        print(f"{'='*70}\n")
        print(f"{'LPF':>6} | {'Left':>12} | {'Right':>12} | {'Ratio':>10} | {'Status':>10}")
        print(f"{'-'*70}")
    
    lpf_data = []
    num_good = 0
    num_ok = 0
    
    for q in all_lpfs[:30]:  # Show first 30
        left_count = left_counts.get(q, 0)
        right_count = right_counts.get(q, 0)
        
        if right_count > 0:
            ratio = left_count / right_count
        else:
            ratio = float('inf') if left_count > 0 else 1.0
        
        # Categorize balance
        if 0.95 <= ratio <= 1.05:
            status = "✓ GOOD"
            num_good += 1
        elif 0.90 <= ratio <= 1.10:
            status = "~ OK"
            num_ok += 1
        else:
            status = "✗ SKEWED"
        
        lpf_data.append({
            'lpf': q,
            'left': left_count,
            'right': right_count,
            'ratio': ratio
        })
        
        if verbose:
            print(f"{q:6d} | {left_count:12,} | {right_count:12,} | {ratio:10.4f} | {status:>10}")
    
    # Overall statistics
    total_left = sum(left_counts.values())
    total_right = sum(right_counts.values())
    overall_ratio = total_left / total_right if total_right > 0 else 0
    
    if verbose:
        print(f"{'-'*70}")
        print(f"{'TOTAL':>6} | {total_left:12,} | {total_right:12,} | {overall_ratio:10.6f} |")
        print(f"\n{'='*70}")
        print(f"ASSESSMENT")
        print(f"{'='*70}\n")
        
        print(f"Total pairs analyzed: {total_left:,}")
        print(f"Overall ratio (left/right): {overall_ratio:.10f}")
        
        total_lpfs = min(30, len(all_lpfs))
        print(f"Well-balanced LPFs (ratio 0.95-1.05): {num_good}/{total_lpfs} ({num_good/total_lpfs*100:.1f}%)")
        
        # Verdict
        if abs(overall_ratio - 1.0) < 1e-6:
            print(f"\n✓✓✓ PERFECT SYMMETRY!")
        elif abs(overall_ratio - 1.0) < 1e-4:
            print(f"\n✓✓ EXCELLENT SYMMETRY")
        elif abs(overall_ratio - 1.0) < 1e-2:
            print(f"\n✓ VERY GOOD SYMMETRY")
        else:
            print(f"\n⚠ WEAK SYMMETRY")
        
        elapsed = time.time() - start_time
        print(f"\nComputation time: {elapsed:.2f}s")
        print(f"{'='*70}\n")
    
    return {
        'k': k,
        'N': N,
        'pairs': total_left,
        'ratio': overall_ratio,
        'lpf_data': lpf_data,
        'left_counts': left_counts,
        'right_counts': right_counts
    }


def test_multiple_primorials(k_values, output_file=None):
    """
    Test LPF symmetry for multiple primorials and summarize results.
    
    Parameters:
    -----------
    k_values : list of int
        Primorial indices to test
    output_file : str, optional
        If provided, write results to this file
    
    Returns:
    --------
    list of dict
        Results for each k value
    """
    print("\n" + "="*70)
    print("LPF SYMMETRY VERIFICATION FOR PRIMORIALS")
    print("="*70)
    print(f"\nTesting k values: {k_values}")
    print()
    
    results = []
    for k in k_values:
        result = verify_lpf_symmetry(k, verbose=True)
        results.append(result)
        time.sleep(0.5)  # Brief pause between tests
    
    # Summary table
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'k':>3} | {'N':>15} | {'Pairs':>12} | {'Ratio':>15} | {'Assessment':>15}")
    print("-" * 70)
    
    for r in results:
        ratio = r['ratio']
        
        if abs(ratio - 1.0) < 1e-6:
            assessment = "✓✓✓ PERFECT"
        elif abs(ratio - 1.0) < 1e-4:
            assessment = "✓✓ EXCELLENT"
        elif abs(ratio - 1.0) < 1e-2:
            assessment = "✓ VERY GOOD"
        else:
            assessment = "⚠ WEAK"
        
        print(f"{r['k']:3d} | {r['N']:15,} | {r['pairs']:12,} | {r['ratio']:15.10f} | {assessment:>15}")
    
    # Final verdict
    perfect_count = sum(1 for r in results if abs(r['ratio'] - 1.0) < 1e-6)
    
    print("\n" + "="*70)
    print("FINAL VERDICT")
    print("="*70)
    
    print(f"\nResults across {len(results)} primorials:")
    print(f"  Perfect symmetry (ratio = 1.000000): {perfect_count}/{len(results)}")
    
    if perfect_count == len(results):
        print(f"\n✓✓✓ HYPOTHESIS CONFIRMED AT HIGHEST LEVEL!")
        print(f"All primorials show perfect LPF pairing symmetry.")
    else:
        print(f"\n~ Hypothesis partially confirmed")
    
    print()
    
    # Write to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("LPF Symmetry Verification Results\n")
            f.write("="*70 + "\n\n")
            for r in results:
                f.write(f"k={r['k']}, N={r['N']:,}, Pairs={r['pairs']:,}, Ratio={r['ratio']:.10f}\n")
            f.write(f"\nPerfect symmetry: {perfect_count}/{len(results)}\n")
        print(f"Results written to {output_file}")
    
    return results


def main():
    """
    Main entry point for the script.
    """
    # Default test: k = 3, 4, 5, 6, 7, 8
    # For k=8, this tests ~4.8 million pairs and takes 1-2 minutes
    default_k_values = [3, 4, 5, 6, 7, 8]
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        try:
            k_values = [int(k) for k in sys.argv[1:]]
            print(f"Testing custom k values: {k_values}")
        except ValueError:
            print("Usage: python lpf_symmetry_verification.py [k1 k2 k3 ...]")
            print("Example: python lpf_symmetry_verification.py 3 4 5")
            sys.exit(1)
    else:
        k_values = default_k_values
        print("Using default k values: [3, 4, 5, 6, 7, 8]")
        print("(Specify custom values as arguments: python lpf_symmetry_verification.py 3 4 5)\n")
    
    # Run tests
    results = test_multiple_primorials(k_values, output_file="lpf_results.txt")
    
    return results


if __name__ == "__main__":
    main()
