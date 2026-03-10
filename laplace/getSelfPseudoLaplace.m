function [Uii,Yii] = getSelfPseudoLaplace(P,rin_vec,rout_vec,pair_points)
%GETSELFPSEUDOLAPLACE Build one-body pseudoinverse factors for Laplace SLP.
%
% Syntax:
%   [Uii,Yii] = getSelfPseudoLaplace(P,rin_vec,rout_vec,pair_points)
%
% Inputs:
%   P          - Number of particles.
%   rin_vec    - Source nodes (stacked by particle).
%   rout_vec   - Collocation nodes (stacked by particle).
%   pair_points- P x 2 bookkeeping matrix. Uses column 2 for colloc counts.
%
% Outputs:
%   Uii, Yii   - Cell arrays with pseudoinverse factors.

if nargin < 4 || isempty(pair_points)
    pair_points = [zeros(P,1) repmat(length(rout_vec)/P,P,1)];
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
    [Y,U] = getPseudoFactors(A,tol,0);

    Uii{i} = U';
    Yii{i} = Y;

    start_colloc = start_colloc + M;
end

end
