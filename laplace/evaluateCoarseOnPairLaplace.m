function Npair = evaluateCoarseOnPairLaplace(q,rbase_in_c,rout_f)
%EVALUATECOARSEONPAIRLAPLACE Cross-evaluation map for scalar pair correction.
%
% Syntax:
%   Npair = evaluateCoarseOnPairLaplace(q,rbase_in_c,rout_f)
%   evaluateCoarseOnPairLaplace() % self-test
%
% Anna Broms, Mar 2026

if nargin==0
    test_evaluateCoarseOnPairLaplace;
    return
end

q = q(:);
rbase_in_c = rbase_in_c(:);
rout_f = rout_f(:);

assert(numel(q)==2,'q must contain two centers.');
assert(mod(numel(rout_f),2)==0,'rout_f must have even length.');

nhalf = numel(rout_f)/2;
rout_1 = rout_f(1:nhalf);
rout_2 = rout_f(nhalf+1:end);

N21 = lapSLPmat(rbase_in_c+q(1),rout_2); % body1 -> body2
N12 = lapSLPmat(rbase_in_c+q(2),rout_1); % body2 -> body1

Npair = [zeros(size(N12)) N12; N21 zeros(size(N12))];

end

function test_evaluateCoarseOnPairLaplace
fprintf('evaluateCoarseOnPairLaplace self-test\n');

q = [0;2.01];
Nc = 60;
Nt = ceil(1.2*Nc);

tc = linspace(0,2*pi,Nc+1)'; tc(end)=[];
t1 = linspace(0,2*pi,Nt+1)'; t1(end)=[];

rbase_in_c = 0.55*(cos(tc)+1i*sin(tc));
rout_1 = q(1)+cos(t1)+1i*sin(t1);
rout_2 = q(2)+cos(t1)+1i*sin(t1);
rout_f = [rout_1; rout_2];

Npair = evaluateCoarseOnPairLaplace(q,rbase_in_c,rout_f);

lam1 = randn(Nc,1);
lam2 = randn(Nc,1);
rhs = -Npair*[lam2; lam1];

u1 = lapSLPmat(rbase_in_c+q(2),rout_1)*lam1;
u2 = lapSLPmat(rbase_in_c+q(1),rout_2)*lam2;
ref = -[u1;u2];

rel = norm(rhs-ref,inf)/max(1,norm(ref,inf));
fprintf('  rel. err: %.3e\n',rel);

end
