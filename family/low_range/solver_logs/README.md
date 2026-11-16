# Solver logs (16 November 2026)

| log | command (run from `family/`) | result |
|---|---|---|
| `rho_175.log`, `rho_210.log` | `rho_cpsat.py 175 46 20 600 0,1,4,5 0,1,2 18`, `rho_cpsat.py 210 54 30 600 0,1,6 0,1 18` | LS₅, LS₆ found in < 1 s (structure of Prop. 5.4 fixed) |
| `rho_245.log` | `rho_cpsat.py 245 64 42 900 0,4 0 18` | LS₇ INFEASIBLE, 849 s |
| `rho_245_confirm.log` | `WORKERS=8 SEED=7 rho_cpsat.py 245 64 42 3000 1,5 2 18` (equivalent normalization) | LS₇ INFEASIBLE, 125 s |
| `rho_525.log`, `rho_490.log`, `rho_455.log` | `rho_cpsat.py k λ μ 600 J0 I0 36` with the J0, I0 of Table 3 | UNKNOWN after 600 s (these PDSs exist) |
| `sigma_280.log`, `sigma_315.log`, `sigma_350.log` | `sigma_cpsat.py 280 76 56 1800 - 27 18`, `… 315 90 72 1800 0 27 18`, `… 350 106 90 1800 0,1 27 18` | LS₈, LS₉, LS₁₀ UNKNOWN after 1800 s |
| `low_175.log`, `low_210.log`, `low_245.log` | `divcut.py k λ μ 0 900 3` (unrestricted C5 model in (Z/6)⁴ coordinates; this folder) | LS₅ found 15 s, LS₆ found 70 s, LS₇ UNKNOWN 900 s |
| `run525_cut0.log`, `run525_cut1.log` | `divcut.py 525 216 210 {0,1} 300 6` | UNKNOWN 300 s without and with divisibility cuts (Section 3.3) |

The non-existence of a ρ-invariant LS₇(36) PDS rests on the solver-free `../exhaustive_m7.out`. These CP-SAT runs only
confirm it.
