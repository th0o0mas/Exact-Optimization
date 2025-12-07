clc; clear;

%% --- Settings ---
n_iter = 50;  % number of iterations
options = optimoptions('quadprog','Display','off');

%% -------------------- Example 4.5: Original 3-variable equality QP --------------------
disp('================ Example 4.5: 3-variable equality QP ================');
P1 = [1 -1 1; -1 2 -2; 1 -2 4];
q1 = [-7;-12;-15];
s1 = 0;
Aeq1 = [1 1 1];
beq1 = 3;
Cineq1 = [];
dineq1 = [];

% quadprog
tic;
for i=1:n_iter
    [x_qp1,fval_qp1] = quadprog(P1,q1,[],[],Aeq1,beq1,[],[],[],options);
end
time_qp1 = toc;

% ExOpt
tic;
for i=1:n_iter
    [x_global1,f_global1,X_term1,F_term1,I_term1] = Optimization_function(P1,q1,s1,Aeq1,beq1,Cineq1,dineq1);
end
time_ex1 = toc;

disp('Quadprog solution:'); disp(x_qp1);
disp('Quadprog objective:'); disp(fval_qp1);
disp('ExOpt global solution:'); disp(x_global1);
disp('ExOpt global objective:'); disp(f_global1);
disp('ExOpt terminal solutions (columns):'); disp(X_term1);
disp('ExOpt terminal objectives:'); disp(F_term1);
disp('Active inequality sets:'); disp(I_term1);
disp(['Time quadprog: ', num2str(time_qp1),' s']);
disp(['Time ExOpt: ', num2str(time_ex1),' s']);
disp('------------------------------------------------------------');

%% -------------------- Example 2: Example 4.1 (inequalities only) --------------------
disp('================ Example 2: Example 4.1 (inequalities only) ================');
P2 = [4 1; 1 2];
q2 = [-12;-10];
s2 = 0;
Aeq2 = [];
beq2 = [];
Cineq2 = [1 1; -1 0; 0 -1];
dineq2 = [4;0;0];

% quadprog
tic;
for i=1:n_iter
    [x_qp2,fval_qp2] = quadprog(P2,q2, Cineq2,dineq2,[],[],[],[],[],options);
end
time_qp2 = toc;

% ExOpt
tic;
for i=1:n_iter
    [x_global2,f_global2,X_term2,F_term2,I_term2] = Optimization_function(P2,q2,s2,Aeq2,beq2,Cineq2,dineq2);
end
time_ex2 = toc;

disp('Quadprog solution:'); disp(x_qp2);
disp('Quadprog objective:'); disp(fval_qp2);
disp('ExOpt global solution:'); disp(x_global2);
disp('ExOpt global objective:'); disp(f_global2);
disp('ExOpt terminal solutions (columns):'); disp(X_term2);
disp('ExOpt terminal objectives:'); disp(F_term2);
disp('Active inequality sets:'); disp(I_term2);
disp(['Time quadprog: ', num2str(time_qp2),' s']);
disp(['Time ExOpt: ', num2str(time_ex2),' s']);
disp('------------------------------------------------------------');

%% -------------------- Example 4.2: 3-variable singular P, equality + inequalities --------------------
disp('================ Example 3: Singular P with equality + inequalities ================');
P3 = diag([1 0 0]);
q3 = [0;0;0];
s3 = 0;
Aeq3 = [0 0 1];
beq3 = 2;
Cineq3 = [1 0 0; 0 1 0];
dineq3 = [4;4];

% quadprog
tic;
for i=1:n_iter
    [x_qp3,fval_qp3] = quadprog(P3,q3, Cineq3,dineq3,Aeq3,beq3,[],[],[],options);
end
time_qp3 = toc;

% ExOpt
tic;
for i=1:n_iter
    [x_global3,f_global3,X_term3,F_term3,I_term3] = Optimization_function(P3,q3,s3,Aeq3,beq3,Cineq3,dineq3);
end
time_ex3 = toc;

disp('Quadprog solution:'); disp(x_qp3);
disp('Quadprog objective:'); disp(fval_qp3);
disp('ExOpt global solution:'); disp(x_global3);
disp('ExOpt global objective:'); disp(f_global3);
disp('ExOpt terminal solutions (columns):'); disp(X_term3);
disp('ExOpt terminal objectives:'); disp(F_term3);
disp('Active inequality sets:'); disp(I_term3);
disp(['Time quadprog: ', num2str(time_qp3),' s']);
disp(['Time ExOpt: ', num2str(time_ex3),' s']);
disp('------------------------------------------------------------');
