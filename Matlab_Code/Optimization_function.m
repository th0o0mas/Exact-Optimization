function [x_global,f_global,X_term,F_term,I_term] = ExOpt_Algorithm3_1(P,q,s,Aeq,beq,Aineq,bineq)

% ==============================================================
% Exact Optimization using Algorithm 3.1 (Lin 2023)
% Full Output Version (compact + stable)
%
% Inputs:
%   P,q,s        Quadratic cost  (1/2 x'Px + q'x + s)
%   Aeq, beq     Equality constraints Aeq x = beq
%   Aineq,bineq  Inequality constraints Aineq x <= bineq
%
% Outputs:
%   x_global     Best (minimum-cost) optimal solution
%   f_global     Corresponding objective value
%   X_term       All terminal optimal points (one per active set)
%   F_term       Their objective values
%   I_term       Active-set indices for each terminal solution
% ==============================================================

%% Preprocessing
n = length(q);
neq = size(Aeq,1);
nine = size(Aineq,1);

if isempty(Aineq)
    nine = 0;
end

X_term = [];
F_term = [];
I_term = {};

%% Case 1 — Equality-only (no inequalities)
if nine == 0
    % Projection-based solution
    A = Aeq; b = beq;
    Ap = pinv(A);
    
    % Nullspace basis
    [~,~,V] = svd(A);
    nullMask = sum(abs(A*V) < 1e-12,1)==neq;
    V2 = V(:,nullMask);

    % Reduced KKT (Theorem 3.2)
    M1 = V2' * P * V2;
    M2 = V2' * (q + P*Ap*b);

    if norm(M1) > 1e-10
        y = -pinv(M1) * M2;
        x = Ap*b + V2*y;
    else
        % constant solution set
        x = Ap*b;
    end

    f = 0.5*x'*P*x + q'*x + s;

    % Outputs
    x_global = x;
    f_global = f;
    X_term = x;
    F_term = f;
    I_term = {[]};
    return;
end

%% Case 2 — Inequalities present
% We enumerate all active sets (subset of inequalities)
for k = 0:nine
    combs = nchoosek(1:nine,k);

    for r = 1:size(combs,1)

        active = combs(r,:);
        
        % Build combined equality system
        A = [Aeq; Aineq(active,:)]; 
        b = [beq; bineq(active)];

        % Check if feasible (rank consistency)
        if rank(A) < rank([A b])
            continue;  % no solution to this active set
        end
        
        % Build particular solution + nullspace
        Ap = pinv(A);
        [~,~,V] = svd(A);
        nullMask = sum(abs(A*V) < 1e-12,1)==size(A,1);
        V2 = V(:,nullMask);

        % Reduced matrices
        M1 = V2'*P*V2;
        M2 = V2'*(q + P*Ap*b);

        if norm(M1) > 1e-10
            y = -pinv(M1)*M2;
            x = Ap*b + V2*y;
        else
            % Only feasible if M2 == 0
            if norm(M2) > 1e-10
                continue;
            end
            x = Ap*b;
        end

        % Check inactive inequalities
        if any(Aineq*x - bineq > 1e-8)
            continue; 
        end

        % Save terminal solution
        fval = 0.5*x'*P*x + q'*x + s;
        X_term = [X_term x];
        F_term = [F_term fval];
        I_term{end+1} = active;
    end
end

%% Final selection of global minimum
if isempty(F_term)
    error('No feasible terminal points found.');
end

[f_global, idx] = min(F_term);
x_global = X_term(:,idx);

end
