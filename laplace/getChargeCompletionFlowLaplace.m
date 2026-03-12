function [lambda0,u_rhs] = getChargeCompletionFlowLaplace(rsrc,rtar,source_ind,Q_body,use_fmm)
%GETCHARGECOMPLETIONFLOWLAPLACE Build completion sources for prescribed net charges.
%
% Syntax:
%   [lambda0,u_rhs] = getChargeCompletionFlowLaplace(rsrc,rtar,source_ind,Q_body,use_fmm)
%
% Inputs:
%   rsrc       - Complex source locations.
%   rtar       - Complex target locations where the completion field is evaluated.
%   source_ind - Cell array with source indices for each body.
%   Q_body     - Prescribed net charge per body.
%   use_fmm    - Logical flag for lapSLPfield.
%
% Outputs:
%   lambda0    - Source strengths with per-body sums equal to Q_body.
%   u_rhs      - Completion potential at targets, u_rhs = -SLP(lambda0).
%
% See also: solve_elast_1B, solve_elast_2B, solve_elast_peanut.
%
% Anna Broms, Mar 2026

if nargin < 5 || isempty(use_fmm)
    use_fmm = true;
end

Q_body = Q_body(:);
P = numel(source_ind);

assert(numel(Q_body)==P,'Q_body must have one entry per particle.');

lambda0 = zeros(numel(rsrc),1);

for k = 1:P
    idx = source_ind{k};
    lambda0(idx) = Q_body(k)/numel(idx);
end

u_rhs = -lapSLPfield(rsrc,rtar,lambda0,use_fmm);

end
