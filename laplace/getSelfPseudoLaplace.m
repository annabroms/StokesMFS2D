function [Uii,Yii] = getSelfPseudoLaplace(P,rin_vec,rout_vec,pair_points,project_charge)
%GETSELFPSEUDOLAPLACE Build one-body pseudoinverse factors for Laplace SLP.
%
% Syntax:
%   [Uii,Yii] = getSelfPseudoLaplace(P,rin_vec,rout_vec,pair_points)
%   [Uii,Yii] = getSelfPseudoLaplace(P,rin_vec,rout_vec,pair_points,project_charge)
%
% Inputs:
%   P          - Number of particles.
%   rin_vec    - Source nodes (stacked by particle).
%   rout_vec   - Collocation nodes (stacked by particle).
%   pair_points- P x 2 bookkeeping matrix. Uses column 2 for colloc counts.
%
% Outputs:
%   Uii, Yii   - Cell arrays with pseudoinverse factors.
%
% See also: lapSLPmat, getPairBlockLaplace.
%
% Anna Broms, Mar 2026

if nargin < 4 || isempty(pair_points)
    pair_points = [zeros(P,1) repmat(length(rout_vec)/P,P,1)];
end
if nargin < 5 || isempty(project_charge)
    project_charge = false;
else
    project_charge = logical(project_charge);
end

tol = 1e-14;

Uii = cell(P,1);
Yii = cell(P,1);

N = length(rin_vec)/P;
start_colloc = 0;

for i = 1:P
    M = pair_points(i,2);
    rin = rin_vec((i-1)*N+1:i*N);
    rout = rout_vec(start_colloc+1:start_colloc+M);

    A = lapSLPmat(rin,rout);
    if project_charge
        A = apply_charge_projection_with_closure(A,M,N);
    end
    [Y,U] = getPseudoFactors(A,tol,0);

    Uii{i} = U';
    Yii{i} = Y;

    start_colloc = start_colloc + M;
end

function Aproj = apply_charge_projection_with_closure(A,m,n)
% Apply A*(I-Lq)+Lr for Kq=ones(n,1), Kout=ones(m,1).

Aproj = A;

if n > 0
    row_mean = sum(Aproj,2)/n;
    Aproj = Aproj - row_mean*ones(1,n);
end

if m > 0 && n > 0
    Aproj = Aproj + 1;
end

end

end
