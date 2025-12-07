function [x_global,f_global,X_term,F_term,I_term] = Optimization_function(P,q,s,Aeq,beq,Aineq,bineq)
% OPTIMIZATION_FUNCTION Exact Quadratic Optimization (Algorithm 3.1 Lin 2023)
% Solves: min 1/2 x'Px + q'x + s  s.t. Aeq x = beq, Aineq x <= bineq
%
% Robustness: Uses pseudoinverse for singular Hessians (semi-definite P).
% Stability : Standardized tolerances for rank and KKT checks.

    %% 1. Initialization and Constants
    TOL_RANK = 1e-12; % Tolerance for Null-space detection
    TOL_KKT  = 1e-10; % Tolerance for KKT stationarity
    TOL_FEAS = 1e-8;  % Tolerance for inequality feasibility

    % Ensure consistent shapes
    q = q(:); 
    if isempty(beq), beq = zeros(0,1); end
    if isempty(bineq), bineq = zeros(0,1); end
    
    nine = size(Aineq, 1);
    X_term = []; F_term = []; I_term = {};

    %% 2. Active Set Enumeration
    % Iterate k from 0 (no inequality active) to nine (all active)
    for k = 0:nine
        
        % Generate combinations of active inequalities
        if k == 0, combs = zeros(1,0); else, combs = nchoosek(1:nine, k); end

        for r = 1:size(combs, 1)
            active_idx = combs(r, :);
            
            % 2.1 Build current Equality System (Original Eq + Active Ineq)
            if isempty(active_idx)
                A = Aeq; b = beq;
            else
                A = [Aeq; Aineq(active_idx, :)];
                b = [beq; bineq(active_idx)];
            end
            
            % 2.2 Solve Equality-Constrained Sub-Problem
            % System: [P A'; A 0] * [x; lambda] = [-q; b]
            % We use the Null-Space method for numerical stability.
            
            Ap = pinv(A); 
            [~, ~, V] = svd(A);
            
            % Feasibility Check: Does Ax=b have ANY solution?
            % If rank(A) < rank([A b]), it's inconsistent.
            if ~isempty(A) && rank(A) < rank([A, b])
                continue; 
            end
            
            % Identify Null Space Basis (directions where A*x = 0)
            if isempty(A)
                V2 = eye(length(q)); % No constraints, full space is null space
                x_part = zeros(length(q),1);
            else
                nullMask = sum(abs(A*V) < TOL_RANK, 1) == size(A,1);
                V2 = V(:, nullMask);
                x_part = Ap*b; % Particular solution
            end

            % Project P and q onto the null space
            M1 = V2' * P * V2;            % Projected Hessian
            M2 = V2' * (q + P * x_part);  % Projected Gradient

            % Solve Reduced KKT: M1 * y = -M2
            if norm(M1) > TOL_KKT
                % Standard case: Curvature exists, unique minimum in null space
                y = -pinv(M1) * M2;
                x = x_part + V2*y;
            else
                % Singular case (Flat valley):
                % Solution exists only if gradient is zero along flat directions
                if norm(M2) > TOL_KKT
                    continue; % Unbounded or no stationary point
                end
                x = x_part; % Pick minimum norm particular solution
            end

            % 2.3 Check Inactive Inequalities
            % The solution must satisfy ALL inequalities, not just active ones
            if isempty(Aineq) || all(Aineq*x - bineq <= TOL_FEAS)
                fval = 0.5*x'*P*x + q'*x + s;
                
                % Store result
                X_term = [X_term, x];     %#ok<*AGROW>
                F_term = [F_term, fval];
                I_term{end+1} = active_idx;
            end
        end
    end

    %% 3. Final Selection
    if isempty(F_term)
        error('Optimization_function:NoFeasiblePoint', 'No feasible solution found.');
    end
    [f_global, idx] = min(F_term);
    x_global = X_term(:, idx);
end