function [U,Y,Atot] = getTripletBlockStokes(rin_trip,rout_trip,Lf_trip,Lr_trip,svd_opts)
%GETTRIPLETBLOCKSTOKES Build pseudoinverse factors for a 3-body Stokes block.

if nargin < 5 || isempty(svd_opts)
    svd_opts = struct();
end

S = stokSLPmat(rin_trip,rout_trip,1);
if isempty(Lr_trip)
    Atot = S;
else
    Atot = S - S*Lf_trip + Lr_trip;
end

tol = 1e-11;
tol = 1e-9;
[Y,U] = getPseudoFactors(Atot,tol,0,svd_opts);

end
