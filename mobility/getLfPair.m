function Lf = getLfPair(K1,K2)
%GETLFPAIR Constructs pair projection matrix for fine Stokeslet sources.
%
% Syntax:
%   Lf = getLfPair(K1,K2)
%   getLfPair()    % run built-in self-test
%
% Description:
%   Builds the pair projection matrix using per-particle blocks
%   K*((K'*K)\K') and assembles them in the pair ordering
%   [x_1; x_2; y_1; y_2], matching getLrPair.
%
% Anna Broms, Feb 2026

if nargin == 0
    run_self_test();
    if nargout > 0
        Lf = [];
    end
    return;
end

Lf1 = K1*((K1'*K1)\K1');
Lf2 = K2*((K2'*K2)\K2');

xx1 = Lf1(1:end/2,1:end/2);
xx2 = Lf2(1:end/2,1:end/2);
yy1 = Lf1(end/2+1:end,end/2+1:end);
yy2 = Lf2(end/2+1:end,end/2+1:end);
yx1 = Lf1(end/2+1:end,1:end/2);
yx2 = Lf2(end/2+1:end,1:end/2);
xy1 = Lf1(1:end/2,end/2+1:end);
xy2 = Lf2(1:end/2,end/2+1:end);

Lfxx = [xx1 zeros(size(xx1,1),size(xx2,2)); zeros(size(xx2,1),size(xx1,2)) xx2];
Lfyy = [yy1 zeros(size(yy1,1),size(yy2,2)); zeros(size(yy2,1),size(yy1,2)) yy2];
Lfxy = [xy1 zeros(size(xy1,1),size(xy2,2)); zeros(size(xy2,1),size(xy1,2)) xy2];
Lfyx = [yx1 zeros(size(yx1,1),size(yx2,2)); zeros(size(yx2,1),size(yx1,2)) yx2];

Lf = [Lfxx Lfxy; Lfyx Lfyy];

end

function run_self_test()
% Checks:
% 1) Lf is an (approximately) orthogonal projector.
% 2) (I-Lf) removes rigid-body force/torque moments particle-wise.

rng(1);
tol = 1e-11;

q1 = 0.0 + 0.0i;
q2 = 2.3 + 0.2i;

Nf1 = 24;
Nf2 = 28;
t1 = linspace(0,2*pi,Nf1+1)'; t1(end) = [];
t2 = linspace(0,2*pi,Nf2+1)'; t2(end) = [];

r1 = q1 + 1.20*(cos(t1) + 1i*sin(t1));
r2 = q2 + 1.05*(cos(t2) + 1i*sin(t2));

K1 = getKmat2D(r1,q1);
K2 = getKmat2D(r2,q2);
Lf = getLfPair(K1,K2);

n1 = size(r1,1);
n2 = size(r2,1);
n = 2*(n1+n2);

assert(isequal(size(Lf),[n n]),'getLfPair returned wrong size.');

proj_err = norm(Lf*Lf - Lf,inf)/max(1,norm(Lf,inf));
sym_err = norm(Lf - Lf',inf)/max(1,norm(Lf,inf));

% Random pair-ordered source vector: [x1; x2; y1; y2].
beta = randn(n,1);
beta_proj = beta - Lf*beta;

ind_x1 = 1:n1;
ind_x2 = n1+1:n1+n2;
ind_y1 = n1+n2+1:n1+n2+n1;
ind_y2 = n1+n2+n1+1:n;

b1 = [beta_proj(ind_x1); beta_proj(ind_y1)];
b2 = [beta_proj(ind_x2); beta_proj(ind_y2)];

mom1 = norm(K1'*b1,inf);
mom2 = norm(K2'*b2,inf);

fprintf('getLfPair self-test:\n');
fprintf('  idempotence rel. err  : %.3e\n',proj_err);
fprintf('  symmetry rel. err     : %.3e\n',sym_err);
fprintf('  particle 1 moment err : %.3e\n',mom1);
fprintf('  particle 2 moment err : %.3e\n',mom2);

assert(proj_err < tol,'Idempotence check failed.');
assert(sym_err < tol,'Symmetry check failed.');
assert(mom1 < tol,'Particle 1 projection check failed.');
assert(mom2 < tol,'Particle 2 projection check failed.');

fprintf('  PASSED\n');
end
