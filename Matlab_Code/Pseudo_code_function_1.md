# Exact Optimization 1
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

## Pseudocode


1) Initialize
- n := length(q)
- nine := number of rows in Aineq (set to 0 if Aineq is empty)
- X_term := empty list/matrix (store solutions as columns)
- F_term := empty list/vector
- I_term := empty list (store active-set indices)

2) Equality-only case (no inequalities)
- If nine == 0:
  - Find a particular solution xp satisfying Aeq * xp = beq.
  - Find a nullspace basis N of Aeq (directions that keep Aeq * x = beq).
  - Represent all feasible x as: x = xp + N * y.
  - Minimize the quadratic in y:
    - If the reduced quadratic (Hessian along y) is nonzero, solve for optimal y and compute x.
    - Otherwise, the cost is flat along N; choose xp.
  - Compute f = 0.5 * xᵀ P x + qᵀ x + s.
  - Set:
    - x_global := x
    - f_global := f
    - X_term := [x]
    - F_term := [f]
    - I_term := [empty set]
  - Return outputs.

3) Inequalities present (nine > 0)
- Enumerate all possible active sets (subsets of {1..nine}):
  - For each subset `active`:
    a) Treat active inequalities as equalities:
       - Build combined system A * x = b by stacking:
         - A := [Aeq; Aineq(active, :)]
         - b := [beq; bineq(active)]
    b) Check feasibility of A * x = b:
       - If inconsistent, skip this active set.
    c) Parameterize feasible x:
       - Find a particular solution xp with A * xp = b.
       - Find a nullspace basis N of A.
       - Feasible x are: x = xp + N * y.
    d) Minimize the quadratic over y:
       - Build reduced quadratic in y.
       - If reduced Hessian is nonzero: solve for optimal y and compute x.
       - Else: require reduced linear term to be zero; if not, skip. If yes, take x = xp.
    e) Check all inactive inequalities:
       - If any inactive inequality violates Aineq[i] * x ≤ bineq[i], skip.
    f) If feasible:
       - Compute f = 0.5 * xᵀ P x + qᵀ x + s.
       - Append:
         - x to X_term (as a new column),
         - f to F_term,
         - active indices to I_term.

4) Select the best solution
- If X_term is empty: report error “No feasible terminal points found.”
- Otherwise:
  - Find index `idx` of the smallest value in F_term.
  - Set:
    - x_global := X_term[:, idx]
    - f_global := F_term[idx]

5) Return x_global, f_global, X_term, F_term, I_term

---

## Notes

- “Particular solution xp” can be obtained with a pseudoinverse or any method that solves A * x = b.
- “Nullspace basis N” can be computed via SVD or any nullspace routine.
- Numerical tolerances (e.g., small thresholds) help decide if matrices are “effectively zero” and to accept tiny violations due to floating-point error.
- Active-set enumeration is exhaustive: it tries every subset of inequalities as potentially tight, solves the equality-reduced problem, then verifies the remaining inequalities.