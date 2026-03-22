function res = matvecStokes1BEnhancedLeft(lambda,geom,basis)
%MATVECSTOKES1BENHANCEDLEFT Apply the left-preconditioned enhanced 1B matvec.
%
% Syntax:
%   res = matvecStokes1BEnhancedLeft(lambda,geom,basis)
%
% Inputs:
%   lambda - Stacked source strengths [lambda_x; lambda_y].
%   geom   - Geometry struct from prepareStokes1BEnhanced.
%   basis  - Basis struct from prepareStokes1BEnhanced.

lambda = lambda(:);
ns_tot = geom.total_source_count;
if numel(lambda) ~= 2*ns_tot
    error('matvecStokes1BEnhancedLeft:badInputSize', ...
        'lambda must have length 2*geom.total_source_count.');
end

lambda_x = lambda(1:ns_tot);
lambda_y = lambda(ns_tot+1:end);

u = getVelocityField(geom.rvec_in,geom.rout,lambda_x,lambda_y, ...
    [],[],[],[],[],[],[],geom.use_fmm);
res = applyStokes1BEnhancedPseudo(u,geom,basis);

end
