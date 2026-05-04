function res = getVelocityField(rvec_in,rcheck,stok_x,stok_y,rimage,nimage,rot,stress_x,stress_y,pot_x,pot_y,use_fmm)
%GETVELOCITYFIELD Evaluate the 2D Stokes velocity field from configured source types.
%
% Syntax:
%   res = getVelocityField(rvec_in,rcheck,stok_x,stok_y)
%   res = getVelocityField(rvec_in,rcheck,stok_x,stok_y,use_fmm)
%   res = getVelocityField(rvec_in,rcheck,stok_x,stok_y,rimage,nimage,...
%                          rot,stress_x,stress_y,pot_x,pot_y)
%   res = getVelocityField(...,use_fmm)
%
% Inputs:
%   rvec_in   - Complex source locations for Stokeslets.
%   rcheck    - Complex target locations.
%   stok_x    - x-strengths for Stokeslets at rvec_in.
%   stok_y    - y-strengths for Stokeslets at rvec_in.
%   rimage    - Complex image-source locations (optional).
%   nimage    - Complex image normals/directions (optional).
%   rot       - Rotlet strengths at image sources (optional).
%   stress_x  - x-strengths for stresslets at image sources (optional).
%   stress_y  - y-strengths for stresslets at image sources (optional).
%   pot_x     - x-strengths for potential dipoles at image sources (optional).
%   pot_y     - y-strengths for potential dipoles at image sources (optional).
%   use_fmm   - Logical flag for Stokeslet evaluation:
%               true  -> use `stfmm2d` (default)
%               false -> use `stokSLPdirect`
%
% Output:
%   res       - Stacked velocity [u_x; u_y] at target points.
%
% Notes:
%   The `use_fmm` switch applies only to direct Stokeslet sources
%   (`rvec_in`, `stok_x`, `stok_y`). Image-source contributions are evaluated
%   with their existing direct routines.

% Allow shorthand: getVelocityField(rvec_in,rcheck,stok_x,stok_y,use_fmm)
if nargin == 5 && isscalar(rimage) && (islogical(rimage) || isnumeric(rimage))
    use_fmm = logical(rimage);
    rimage = [];
    nimage = [];
    rot = [];
    stress_x = [];
    stress_y = [];
    pot_x = [];
    pot_y = [];
elseif nargin < 12 || isempty(use_fmm)
    use_fmm = true;
end

rvec_in = rvec_in(:);
rcheck = rcheck(:);
stok_x = stok_x(:);
stok_y = stok_y(:);

n_tar = numel(rcheck);
n_src = numel(rvec_in);

if n_src > 0
    if use_fmm
        %fmm_eps = 1e-10; %old setting
        fmm_eps = 1e-9;
        targ = [real(rcheck)'; imag(rcheck)'];
        srcinfo.sources = [real(rvec_in)'; imag(rvec_in)'];
        srcinfo.stoklet = [stok_x'; stok_y'];
        fmm_timer = tic;
        U = stfmm2d(fmm_eps, srcinfo, 0, targ, 1);
        manageSolveTimeMeasurement('add_fmm',toc(fmm_timer));
        res = [U.pottarg(1,:)'; U.pottarg(2,:)']/2/pi;
    else
        [udirect,vdirect] = stokSLPdirect(real(rvec_in),imag(rvec_in),...
            real(rcheck),imag(rcheck),stok_x,stok_y,n_src);
        res = [udirect; vdirect];
    end
else
    res = zeros(2*n_tar,1);
end

if nargin > 4 && ~isempty(rimage)
    u_rot = getRotlets(rot,rimage,rcheck);
    u_stress = getStresslets(stress_x,stress_y,rimage,rcheck,real(nimage),imag(nimage));
    u_pot = getPotdip(pot_x,pot_y,rimage,rcheck);
    res = res + u_rot + u_stress + u_pot;
end

end
