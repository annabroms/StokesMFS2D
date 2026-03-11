function [U,Y,A] = getPairBlockLaplace(rin_pair,rout_f,proj)
%GETPAIRBLOCKLAPLACE Build pseudoinverse factors for Laplace pair block in 2D.
%
% Syntax:
%   [U,Y,A] = getPairBlockLaplace(rin_pair,rout_f)
%   [U,Y,A] = getPairBlockLaplace(rin_pair,rout_f,proj)
%
% Inputs:
%   rin_pair - Pair source nodes [body i; body j].
%   rout_f   - Pair collocation nodes [body i; body j].
%   proj     - Optional struct for elastance projection/closure:
%              proj.project_charge (logical)
%              proj.nsrc = [n1 n2]
%              proj.ntar = [m1 m2]
%
% Notes:
%   If proj.project_charge is true, this routine factorizes the modified
%   operator
%       A = S*(I-Lq_pair) + Lr_pair,
%   but applies it via block means/sums (rank-1 updates), avoiding explicit
%   dense products with Lq/Lr.
%
% See also: getPairBasisLaplace, getSelfPseudoLaplace.
%
% Anna Broms, Mar 2026

if nargin < 3 || isempty(proj)
    proj = struct();
end

if isfield(proj,'project_charge') && ~isempty(proj.project_charge)
    project_charge = logical(proj.project_charge);
else
    project_charge = false;
end

A = lapSLPmat(rin_pair,rout_f);

if project_charge

    n1 = proj.nsrc(1);
    n2 = proj.nsrc(2);
    m1 = proj.ntar(1);
    m2 = proj.ntar(2);

    if size(A,1)~=m1+m2
        error('getPairBlockLaplace:TargetSizeMismatch', ...
            'Target split does not match matrix size.');
    end
    if size(A,2)~=n1+n2
        error('getPairBlockLaplace:SourceSizeMismatch', ...
            'Source split does not match matrix size.');
    end

    A = apply_pair_charge_projection_with_closure(A,n1,n2,m1,m2);
end

[Y,U] = getPseudoFactors(A,1e-14,0);

end

function Aproj = apply_pair_charge_projection_with_closure(A,n1,n2,m1,m2)
% Apply A*(I-Lq_pair) + Lr_pair with scalar per-body charge constraints.

Aproj = A;

% Body 1: Project off contribution to charge 
cols1 = 1:n1;
row_mean1 = sum(Aproj(:,cols1),2)/n1;
Aproj(:,cols1) = Aproj(:,cols1) - row_mean1*ones(1,n1);

% Body 2: Project off contribution to charge
cols2 = n1+1:n1+n2;
row_mean2 = sum(Aproj(:,cols2),2)/n2;
Aproj(:,cols2) = Aproj(:,cols2) - row_mean2*ones(1,n2);

% Add ones to diagonal blocks. This is the Lr_pair matrix
Aproj(1:m1,1:n1) = Aproj(1:m1,1:n1) + 1;
Aproj(m1+1:m1+m2,n1+1:n1+n2) = Aproj(m1+1:m1+m2,n1+1:n1+n2) + 1;


end
