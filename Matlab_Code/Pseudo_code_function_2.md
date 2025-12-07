# Exact Quadratic Optimization — Optimization_function_2 

## Problem

Minimize:
- f(x) = 0.5 * xᵀ P x + qᵀ x + s

Subject to:
- Aeq x = beq
- Aineq x ≤ bineq

Inputs:
- P, q, s: quadratic cost parameters
- Aeq, beq: equality constraints
- Aineq, bineq: inequality constraints

Outputs:
- x_global: best (lowest-cost) feasible solution found
- f_global: objective value at x_global
- X_term: all feasible terminal solutions (each corresponds to one active set)
- F_term: objective values of terminal solutions
- I_term: indices of active inequalities used to produce each terminal solution

---

## Pseudocode (Readable)


Constants (numerical tolerances):
- TOL_RANK := 1e-12  (detect null-space / near-zero rows)
- TOL_KKT  := 1e-10  (KKT stationarity / reduced Hessian check)
- TOL_FEAS := 1e-8   (inequality feasibility slack tolerance)

1) Initialize
- Ensure q is a column vector.
- If beq is empty, set beq := zero-length column.
- If bineq is empty, set bineq := zero-length column.
- nine := number of inequalities (rows in Aineq); if Aineq is empty, nine := 0.
- X_term := empty (stores candidate solutions as columns)
- F_term := empty (stores objective values)
- I_term := empty (stores active-set indices per solution)

2) Enumerate active sets (from no active inequalities to all)
- For k from 0 to nine:
  - Generate all subsets `active_idx` of size k from {1..nine}.
    - If k == 0, the only subset is the empty set.

  For each active_idx:
    a) Build current equality system
       - If active_idx is empty:
         - A := Aeq
         - b := beq
       - Else:
         - A := stackRows(Aeq, Aineq(active_idx, :))
         - b := stackRows(beq, bineq(active_idx))

    b) Check feasibility of A * x = b
       - If A is non-empty and rank(A) < rank([A | b]), the system is inconsistent; skip this active set.

    c) Parameterize the feasible affine space
       - If A is empty:
         - Nullspace basis N := identity (full space is feasible)
         - Particular solution xp := zero vector
       - Else:
         - Compute xp := a particular solution to A * x = b (e.g., via pseudoinverse).
         - Compute a nullspace basis N of A (directions where A * x = 0).

       - Represent feasible x as: x = xp + N * y.

    d) Minimize the quadratic over y (reduced problem)
       - Reduced Hessian:    M1 := Nᵀ P N
       - Reduced gradient:   M2 := Nᵀ (q + P xp)
       - Solve M1 * y = -M2
         - If ||M1|| > TOL_KKT (has curvature):
           - y := -pinv(M1) * M2
           - x := xp + N * y
         - Else (flat along some directions):
           - If ||M2|| > TOL_KKT:
             - No stationary point in the flat directions; skip active set.
           - Else:
             - x := xp (minimum-norm particular solution)

    e) Check all inequalities (not just the active ones)
       - If Aineq is empty OR Aineq * x - bineq ≤ TOL_FEAS elementwise:
         - Accept x as feasible.
         - Compute fval := 0.5 * xᵀ P x + qᵀ x + s
         - Append:
           - x to X_term (as a new column)
           - fval to F_term
           - active_idx to I_term

3) Select the best solution
- If F_term is empty: error “No feasible solution found.”
- Otherwise:
  - idx := index of minimum value in F_term
  - x_global := X_term[:, idx]
  - f_global := F_term[idx]

4) Return x_global, f_global, X_term, F_term, I_term

---

## Notes

- Pseudoinverse (pinv) is used for robustness in singular or semi-definite cases (P or constraint systems).
- Null-space computation can be done via SVD; any reliable method works.
- Tolerances (TOL_RANK, TOL_KKT, TOL_FEAS) control numerical decisions for feasibility, stationarity, and inequality satisfaction.
- Active-set enumeration is exhaustive: try every subset of inequalities as tight, solve the reduced equality problem, then verify all inequalities.