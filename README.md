# Differentially Private Formation Control: Privacy and Network Co-Design

Code for the paper "Differentially Private Formation Control: Privacy and Network Co-Design".

This implements a budgeted biconvex optimization approach for jointly designing network topology and privacy parameters in formation control. The solver optimizes edge weights and local differential privacy noise levels subject to a budget constraint.

## What's in here

- `budgeted_codesign_solver.py` - Main solver using alternating convex search (ACS)
- `budgeted_biconvex_comparison.py` - Compare results across different budget values
- `budgeted_error_bound_comparison.py` - Compare results across different error bounds
- `tikz_utils.py` - Shared utilities for generating TikZ visualizations

## Running it

Just run either comparison script:

```bash
python budgeted_biconvex_comparison.py
python budgeted_error_bound_comparison.py
```

The scripts will run the optimization, show you matplotlib plots, and save TikZ files to the respective directories.

## What it does

The solver uses an alternating convex search approach:
1. Fix privacy parameters, optimize network weights + connectivity
2. Fix network weights, optimize privacy parameters
3. Repeat until convergence

Key constraints:
- Budget: Total edge weight sum (Tr(L)) ≤ 2B
- Error: Privacy-error tradeoff bound
- Connectivity: Auxiliary variable y ≤ λ₂ (algebraic connectivity)

## Dependencies

- numpy
- scipy
- matplotlib
- networkx
- matplot2tikz

## Output

- TikZ `.tex` files for figures
- Interactive matplotlib plots during execution

