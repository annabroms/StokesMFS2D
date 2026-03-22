function res = matvecStokes1BEnhanced(tau,geom,basis)
%MATVECSTOKES1BENHANCED Apply the Stokeslet-only enhanced 1B matvec.
%
% Syntax:
%   res = matvecStokes1BEnhanced(tau,geom,basis)

[lambda_x_raw,lambda_y_raw,lambda_body,lambda_px,lambda_py] = mapStokes1BEnhancedSources(tau,geom,basis);

if geom.solve_resistance
    lambda_x = lambda_x_raw;
    lambda_y = lambda_y_raw;
else
    lambda_x = lambda_px;
    lambda_y = lambda_py;
end

res = getVelocityField(geom.rvec_in,geom.rout,lambda_x,lambda_y,[],[],[],[],[],[],[],geom.use_fmm);

for k = 1:numel(geom.q)
    self_k = basis.Aii{k}*lambda_body{k};
    nt = geom.target_count(k);
    res(geom.target_x_ind{k}) = res(geom.target_x_ind{k}) - self_k(1:nt);
    res(geom.target_y_ind{k}) = res(geom.target_y_ind{k}) - self_k(nt+1:end);
end

res = res + tau(:);

end
