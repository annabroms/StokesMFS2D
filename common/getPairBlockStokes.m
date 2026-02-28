function [U,Y,Atot] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair)
%GETPAIRBLOCKSTOKES Build pseudoinverse factors for a 2-body Stokes pair block.
%
% Syntax:
%   [U,Y] = getPairBlockStokes(rin_pair,rout_f,Lf_pair,Lr_pair)
%
% Inputs:
%   rin_pair  - complex source locations for the pair fine basis (body 1
%               followed by body 2, including any enhancement points).
%   rout_f    - complex collocation locations for the pair solve
%               (body 1 followed by body 2).
%   Lf_pair   - pair projection matrix for fine sources (same ordering as
%               rin_pair); used to remove force/torque-producing modes.
%   Lr_pair   - pair rigid-motion closure matrix (can be empty). If empty,
%               no projection/closure is applied and the raw single-layer
%               operator is factorized.
%
% Outputs:
%   U, Y  - pseudoinverse factors returned by getPseudoFactors for the
%           assembled pair operator Ntot, such that the stabilized solve
%           can be applied as x = Y*(U'*rhs).
% 
% Operator assembled in this routine:
%   S    = singleLayer(rin_pair, rout_f, 1)
%   Atot = S - S*Lf_pair + Lr_pair   (when Lr_pair is non-empty)
%   Atot = S                         (otherwise)
%
% See also: getPseudoFactors, singleLayer, getLrPair, getLfPair
%
% Anna Broms Feb 12, 2026

mu = 1; 

%debug: visualise grid for the pair
debug = 0;
if debug
    figure(8)
    clf; 
    plot(real(rin_pair),imag(rin_pair),'r+')
    hold on
    plot(real(rout_f),imag(rout_f),'bo')
end

if isempty(Lr_pair)
    project = 0; 
else
    project = 1; 
end

% Compute the dense target from source matrix.
S = singleLayer(rin_pair,rout_f,mu); 

% In a mobility problem, need to project so that fine sources don't
% contribute to net force/torque on the particles
if project
    Atot = S-S*Lf_pair+Lr_pair;
else
    Atot = S;
end 

tol = 1e-14; %SVD truncation level

[Y,U] = getPseudoFactors(Atot,tol,0); 


end
