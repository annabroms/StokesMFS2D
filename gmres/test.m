clear;
close all;
rng(6);

%From Dhairya Malhotra Sept 2024

%%
N = 400;

S = zeros(N, N);
for i = 1:N
  S(i,i) = eps^((i-1)/N);
end
A = rand(N)*S*rand(N) + eye(N);
A_op = @(x) A*x;



% create preconditioner object
precond = KrylovPrecond();

% solve without preconditioner and initialize the preconditioner
x0 = rand(N,1);
b =  A*x0;
[x, e, precond] = precond_gmres(A_op, b, x0*0, N, 1e-12, precond);
fprintf("iterations = %d\n", length(e))
fprintf("error = %e\n\n", max(x-x0))

% solve using preconditioner
x0 = rand(N,1);
b =  A*x0;
[x, e, precond] = precond_gmres(A_op, b, x0*0, N, 1e-10, precond);
fprintf("iterations = %d\n", length(e))
fprintf("error = %e\n\n", max(x-x0))

