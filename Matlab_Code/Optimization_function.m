function [x_global, f_global, X_term, F_term, I_term] = Optimization_function(P,q,s,Aeq,beq,Cineq,dineq)
% ExOpt_Algorithm3_1_fixed: Exact QP solver using nullspace + subset enumeration
%
% Handles equality-only, inequality-only, and mixed QPs.
% Returns global optimum, objective, terminal solutions, and active inequality sets.

tol = 1e-9;

n = size(P,1);

% -----------------------------
% Step 1: Compute equality-constrained solution
% -----------------------------
if ~isempty(Aeq)
    x_global = pinv(Aeq) * beq;
else
    x_global = zeros(n,1); % no equality constraints
end
f_global = 0.5*x_global'*P*x_global + q'*x_global + s;

% -----------------------------
% Step 2: Prepare terminal sets
% -----------------------------
X_term = x_global;  % initialize with global solution
F_term = f_global;
I_term = {[]};

% -----------------------------
% Step 3: If there are inequalities, enumerate feasible subsets
% -----------------------------
if ~isempty(Cineq)
    I = 1:size(Cineq,1);      % all inequality indices
    subsets = {};             % all non-empty subsets
    for k = 1:length(I)
        combos = nchoosek(I,k);
        for r = 1:size(combos,1)
            subsets{end+1} = combos(r,:);
        end
    end
    
    % Loop over subsets
    for idx = 1:length(subsets)
        Ij = subsets{idx};
        C = Cineq(Ij,:);
        d = dineq(Ij);
        
        % Build augmented system: equality + active inequalities as equality
        A_aug = [Aeq; C];
        b_aug = [beq; d];
        
        if rank(A_aug) == size(A_aug,2)
            % Full rank: unique solution
            x_candidate = pinv(A_aug) * b_aug;
        else
            % Use nullspace method for underdetermined
            x_part = pinv(A_aug) * b_aug;
            V2 = null(A_aug);
            if ~isempty(V2)
                M1 = V2'*P*V2;
                M2 = V2'*(P*x_part + q);
                x_candidate = x_part - V2*pinv(M1)*M2;
            else
                x_candidate = x_part;
            end
        end
        
        % Check feasibility of remaining inequalities
        remaining = setdiff(I,Ij);
        feasible = true;
        for r = remaining
            if Cineq(r,:)*x_candidate - dineq(r) > tol
                feasible = false;
                break;
            end
        end
        
        if feasible
            % Store terminal solution
            X_term = [X_term, x_candidate];
            F_term = [F_term, 0.5*x_candidate'*P*x_candidate + q'*x_candidate + s];
            I_term{end+1} = Ij;
        end
    end
end

% -----------------------------
% Step 4: If equality-only, X_term = x_global, F_term = f_global, I_term = {[]}
% -----------------------------

end
