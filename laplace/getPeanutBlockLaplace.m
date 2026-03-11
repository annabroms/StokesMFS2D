function [DC,Y] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut,proj)
%GETPEANUTBLOCKLAPLACE Factorization for scalar peanut compression.
%
% Syntax:
%   [DC,Y] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut)
%   [DC,Y] = getPeanutBlockLaplace(rin_pair_c,rin_pair_f,rout_peanut,proj)
%
% Inputs:
%   proj - Optional struct for elastance projection:
%          proj.project_charge (logical)
%          proj.nsrc_c = [n1c n2c]
%          proj.nsrc_f = [n1f n2f]
%
% See also: getPairBasisLaplace, transform_laplace_peanut, getPairBlockLaplace.
%
% Anna Broms, Mar 2026

if nargin < 4 || isempty(proj)
    proj = struct();
end

if isfield(proj,'project_charge') && ~isempty(proj.project_charge)
    project_charge = logical(proj.project_charge);
else
    project_charge = false;
end

Nf = lapSLPmat(rin_pair_f,rout_peanut);
Nc = lapSLPmat(rin_pair_c,rout_peanut);

if project_charge
    if ~(isfield(proj,'nsrc_c') && numel(proj.nsrc_c)==2)
        error('getPeanutBlockLaplace:BadProjNsrcC', ...
            'proj.nsrc_c must be [n1 n2].');
    end
    if ~(isfield(proj,'nsrc_f') && numel(proj.nsrc_f)==2)
        error('getPeanutBlockLaplace:BadProjNsrcF', ...
            'proj.nsrc_f must be [n1 n2].');
    end

    Nc = apply_pair_charge_projection_columns(Nc,proj.nsrc_c(1),proj.nsrc_c(2));
    Nf = apply_pair_charge_projection_columns(Nf,proj.nsrc_f(1),proj.nsrc_f(2));
end

[Y,U] = getPseudoFactors(Nc,1e-14,0);
DC = U'*Nf;

end

function Aproj = apply_pair_charge_projection_columns(A,n1,n2)
% Apply A*(I-Lq_pair) with scalar per-body charge constraints.

Aproj = A;

% Body 1
if n1 > 0
    cols1 = 1:n1;
    row_mean1 = sum(Aproj(:,cols1),2)/n1;
    Aproj(:,cols1) = Aproj(:,cols1) - row_mean1*ones(1,n1);
end

% Body 2
if n2 > 0
    cols2 = n1+1:n1+n2;
    row_mean2 = sum(Aproj(:,cols2),2)/n2;
    Aproj(:,cols2) = Aproj(:,cols2) - row_mean2*ones(1,n2);
end

end
