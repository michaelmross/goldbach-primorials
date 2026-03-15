# `Primorial LPF Allocation`

`primorial_lpf_allocation.py` computes least-prime-factor (LPF) allocation data for
**primorial Goldbach complements**.

Given

- `N = p_k#` (the product of the first `k` primes), and
- `S = [p_{k+1}, lambda * p_{k+1})`,

it studies the complement set

- `C = { N - p : p in S, p prime }`.

Under the hypothesis that no element of `C` is prime, every complement must have
some least prime factor `q > p_k`, and the exact LPF partition is

- `C_q = #{ p in S : lpf(N-p) = q }`.

The script computes that partition up to a user-chosen cutoff `q_max`, compares
it with a Buchstab-weighted heuristic, and summarizes what remains unresolved.

---

## What the script outputs

For each `(k, lambda)` run, the script prints one CSV row with summary statistics.
It can also print a detailed human-readable diagnostic report.

Typical uses:

- measure how much of the complement set is exactly covered by LPFs `<= q_max`
- inspect which LPF bands dominate
- measure the size and composition of the unresolved tail
- estimate whether leftovers are mostly prime or composite
- study LPF-capacity / reuse behavior among leftover composites

---

## Basic usage

Run one small test:

```bash
primorial_lpf_allocation.py --k-min 20 --k-max 20 --lam 5 --q-max 200000
```

Run a short range and write CSV:

```bash
primorial_lpf_allocation.py --k-min 20 --k-max 30 --lam 5 10 20 --q-max 1000000 --csv lpf_runs.csv
```

Show detailed per-run diagnostics:

```bash
primorial_lpf_allocation.py --k-min 24 --k-max 24 --lam 20 --q-max 2000000 --detailed
```

Show progress for long runs:

```bash
primorial_lpf_allocation.py --k-min 34 --k-max 40 --lam 20 --q-max 2000000 --progress --csv lpf_k34_40.csv
```

---

## Command-line options

### Range and window

- `--k-min INT`  
  Smallest `k` to run.

- `--k-max INT`  
  Largest `k` to run.

- `--lam FLOAT [FLOAT ...]`  
  One or more `lambda` values. The prime window is
  `S = [p_{k+1}, lambda * p_{k+1})`.

### LPF search bound

- `--q-max INT`  
  Largest prime `q` considered in the exact LPF scan.
  LPFs above this bound remain in the leftover tail.

### Output control

- `--csv PATH`  
  Save the CSV summary to a file.

- `--detailed`  
  Print a full per-run report, including top LPFs, coverage table,
  LPF bands, leftover summary, and detailed leftover records.

- `--progress`  
  Print progress updates to `stderr` during long runs.

- `--progress-every INT`  
  Controls progress print cadence.

### Correctness / expensive options

- `--verify`  
  Intended to verify assignments, but in the current script the verification flag
  is passed through and not actively used inside `exact_lpf_allocation()`. Treat
  it as reserved for now.

- `--full-factorization`  
  For composite leftovers, also compute full factorizations with `factorint()`.
  This becomes very slow once `k` gets moderately large.

- `--skip-leftover-lpf`  
  Only classify final leftovers as prime vs composite; do **not** try to compute
  LPFs for composite leftovers in the final tail. This is the fastest setting for
  large `k`, but it leaves several leftover-composite columns blank in the CSV.

---

## How the algorithm works

The script has three main stages.

### 1. Build the primorial and the prime window

For each `k`:

- compute the first `k` primes
- form `N = p_k#`
- let `pk = p_k`, `pk1 = p_{k+1}`
- generate primes in `S = [pk1, lambda * pk1)`

The number of such primes is `pi_S`.

### 2. Exact LPF allocation up to `q_max`

The core routine is `exact_lpf_allocation()`. It performs:

- **Phase I:** exact `q`-scan up to `q_small = min(q_max, S_hi - 1)`
- **Phase II:** for remaining unresolved complements, direct LPF extraction via
  `least_prime_factor_above()` when useful

The implementation uses a batched sparse congruence scan:

- for each prime `q`, only numbers with `p ≡ N (mod q)` can contribute
- this is much faster than independently factoring every `N-p`

### 3. Leftover classification

After all assignments with `q <= q_max`, the script classifies unassigned
complements as leftovers.

Each leftover is tagged as one of:

- `prime`
- `composite`
- `composite_unfactored` (when `--skip-leftover-lpf` is used)

---

## Understanding the CSV columns

The CSV header is generated directly by the script and includes these groups.

### Core run parameters

- `k` — primorial index
- `pk` — `p_k`
- `pk1` — `p_{k+1}`
- `lambda` — window multiplier
- `pi_S` — number of primes in `S`
- `q_max` — LPF cutoff used in this run

### Coverage counts

- `exact_covered` — number of complements assigned an exact LPF `<= q_max`
- `leftovers` — number of unresolved complements after the exact scan
- `leftover_primes` — leftovers that are prime
- `leftover_composites` — leftovers that are composite
- `leftover_prime_frac` — `leftover_primes / leftovers`

### Composite-leftover LPF statistics

These are computed only for leftover composites whose LPF is explicitly found.
If LPFs are skipped or remain unresolved, these cells are blank.

- `median_comp_lpf`, `min_comp_lpf`, `max_comp_lpf`
- `median_lpf_over_pk`, `min_lpf_over_pk`, `max_lpf_over_pk`
- `median_log_lpf_over_log_pk`, `min_log_lpf_over_log_pk`, `max_log_lpf_over_log_pk`
- `median_logexp_lpf`, `min_logexp_lpf`, `max_logexp_lpf`

Interpretation:

- `lpf_over_pk` measures LPF size relative to the primorial threshold `p_k`
- `log_lpf_over_log_pk` measures LPF size on a log scale relative to `p_k`
- `logexp_lpf = log(lpf)/log(m)` measures LPF size relative to the complement `m = N-p`

### Capacity / reuse summary for leftover composites

These are computed from LPF classes among leftover composites with known LPF.

- `cap_num_comp_classes` — number of distinct LPF classes among composite leftovers
- `cap_used` — normalized capacity usage
- `cap_waste` — unused capacity across those classes
- `cap_unit_classes` — classes with exact residue capacity `A_q = 1`
- `cap_two_classes` — classes with `A_q = 2`
- `cap_le3_classes` — classes with `A_q <= 3`
- `cap_reuse_ratio` — average multiplicity per LPF class
- `cap_max_multiplicity` — largest `L_q` among composite-leftover LPF classes
- `cap_singleton_classes` — number of classes with multiplicity `1`

### PK-relative LPF bands for leftover composites

These bands classify leftover-composite LPFs by size relative to `p_k`:

- `pk_band_pk_2pk`
- `pk_band_2pk_5pk`
- `pk_band_5pk_10pk`
- `pk_band_10pk_1e2`
- `pk_band_1e2_1e3`
- `pk_band_1e3_1e4`
- `pk_band_gt1e4`

### Exponent bands for leftover composites

These classify LPFs by `log(lpf)/log(m)`:

- `exp_band_le14` = `<= m^(1/4)`
- `exp_band_14_13` = `(m^(1/4), m^(1/3)]`
- `exp_band_13_04` = `(m^(1/3), m^0.4]`
- `exp_band_04_049` = `(m^0.4, m^0.49]`
- `exp_band_gt049` = `> m^0.49`

### Allocation bands for all assigned LPFs

These summarize the exact LPF allocation `hist_exact` for all complements that
were assigned with `q <= q_max`:

- `alloc_band_pk_2pk`
- `alloc_band_2pk_5pk`
- `alloc_band_5pk_10pk`
- `alloc_band_10pk_20pk`
- `alloc_band_20pk_50pk`
- `alloc_band_50pk_100pk`
- `alloc_band_gt100pk`

These are often the most useful columns for the main LPF-allocation story.

---

## How to populate **all** columns for smaller `k`

For small and moderate `k`, you can usually fill every CSV column by forcing the
script to compute LPFs for the final composite leftovers instead of leaving them
unfactored.

### Recommended settings

Use all of the following:

- **do not** pass `--skip-leftover-lpf`
- optionally pass `--full-factorization` if you want explicit complete factors in
  detailed output
- choose `k` small enough that leftover LPF extraction is still feasible
- choose a `q_max` large enough that the unresolved tail is small

A good practical range is usually:

- `k <= 28` : easy
- `k <= 32` : often feasible
- `k >= 34` : increasingly expensive, depending on `lambda` and `q_max`

### Example: populate all columns for one run

```bash
primorial_lpf_allocation.py --k-min 24 --k-max 24 --lam 20 --q-max 2000000 --csv k24_full.csv
```

This will usually populate:

- leftover prime/composite counts
- leftover composite LPF stats
- capacity columns
- PK-relative leftover LPF bands
- exponent bands
- allocation bands

### Example: small sweep with all columns populated

```bash
primorial_lpf_allocation.py --k-min 20 --k-max 28 --lam 20 --q-max 2000000 --csv lpf_small_full.csv
```

### Example: detailed one-off run

```bash
primorial_lpf_allocation.py --k-min 26 --k-max 26 --lam 20 --q-max 2000000 --detailed
```

### When columns go blank

The following columns become blank when the script cannot or does not compute LPFs
for leftover composites:

- `median_comp_lpf`, `min_comp_lpf`, `max_comp_lpf`
- `median_lpf_over_pk`, `min_lpf_over_pk`, `max_lpf_over_pk`
- `median_log_lpf_over_log_pk`, `min_log_lpf_over_log_pk`, `max_log_lpf_over_log_pk`
- `median_logexp_lpf`, `min_logexp_lpf`, `max_logexp_lpf`
- all `cap_*` columns that depend on known leftover composite LPFs
- `pk_band_*` columns for leftover composites
- `exp_band_*` columns

This usually happens when either:

- you pass `--skip-leftover-lpf`, or
- `k` is large enough that leftover composite LPFs are too expensive to extract

### Strategy for a complete small-`k` dataset

For a fully populated calibration table, use something like:

```bash
primorial_lpf_allocation.py \
  --k-min 20 --k-max 30 \
  --lam 5 10 20 \
  --q-max 2000000 \
  --csv lpf_calibration_full.csv
```

This is a good regime for producing paper-quality fully populated rows.

---

## Recommended workflows

### Fast large-`k` exploratory sweep

Use this when you mainly care about exact-covered counts, allocation bands, and
leftover prime fractions.

```bash
primorial_lpf_allocation.py \
  --k-min 34 --k-max 50 \
  --lam 20 \
  --q-max 2000000 \
  --skip-leftover-lpf \
  --progress \
  --csv lpf_large_fast.csv
```

Pros:

- much faster
- avoids getting stuck on final tail LPFs

Cons:

- many leftover-composite columns will be blank

### Fully populated small-`k` calibration sweep

```bash
primorial_lpf_allocation.py \
  --k-min 20 --k-max 30 \
  --lam 20 \
  --q-max 2000000 \
  --csv lpf_small_complete.csv
```

Pros:

- fills essentially all CSV columns
- useful for tables and diagnostics in the paper

Cons:

- slower than fast mode

### Deep inspection of one case

```bash
primorial_lpf_allocation.py \
  --k-min 28 --k-max 28 \
  --lam 20 \
  --q-max 2000000 \
  --detailed \
  --full-factorization
```

Use sparingly. Full factorization becomes expensive quickly.

---

## Notes and caveats

1. **The script name and module docstring differ.**  
   The file is named `primorial_lpf10.py`, but its top docstring still says
   `primorial_lpf.py`.

2. **`--verify` is currently not active.**  
   The CLI exposes it and `run_batch()` passes it through, but the current
   `exact_lpf_allocation()` implementation does not consume `verify_assignments`.

3. **Large `k` tail extraction is the expensive step.**  
   The hardest part is not the exact `q`-scan but the final LPF extraction for
   the remaining composites.

4. **Blank columns are expected in fast mode.**  
   If you choose speed over full leftover LPF extraction, many leftover-composite
   diagnostics are intentionally left blank.

---

## Minimal example set

Populate nearly everything for a manageable run:

```bash
primorial_lpf_allocation.py --k-min 22 --k-max 22 --lam 20 --q-max 2000000 --csv demo_full.csv
```

Run a paper-style calibration sweep:

```bash
primorial_lpf_allocation.py --k-min 20 --k-max 30 --lam 20 --q-max 2000000 --csv calibration.csv
```

Run a large-`k` exploratory sweep quickly:

```bash
primorial_lpf_allocation.py --k-min 34 --k-max 50 --lam 20 --q-max 2000000 --skip-leftover-lpf --progress --csv lpf_large.csv
```

