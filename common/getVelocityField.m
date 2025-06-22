function res = getVelocityField(rvec_in,rcheck,stok_x,stok_y,rimage,nimage,rot,stress_x,stress_y,pot_x,pot_y)
%GETVELOCITYFIELD Evaluates the 2D Stokes velocity field from various source types
%
% Syntax:
%   res = getVelocityField(rvec_in, rcheck, stok_x, stok_y, rimage, nimage, ...
%                          rot, stress_x, stress_y,pot_x, pot_y,)
%
% Inputs:
%   rvec_in   - Complex vector of source points for Stokeslets 
%   rcheck    - Complex vector of target points where the velocity field is evaluated
%   stok_x    - Real vector of x-force strengths for Stokeslets (same length as rvec_in)
%   stok_y    - Real vector of y-force strengths for Stokeslets
%   rimage    - Complex vector of image source locations 
%   nimage    - Complex vector encoding stresslet "directions" at image points (stored as x + iy)
%   rot       - Real vector of torques (strengths for rotlets)
%   stress_x  - Real vector of x-strengths for stresslet dipoles
%   stress_y  - Real vector of y-strengths for stresslet dipoles
%   pot_x     - Real vector of x-strengths for potential dipoles
%   pot_y     - Real vector of y-strengths for potential dipoles

%
% Output:
%   res       - 2N vector containing the evaluated velocity field at the N target points 
%               (x-velocity components; y-velocity components)
%
% Description:
%   Computes the 2D Stokes velocity field at given target points due to:
%   - Direct Stokeslets at `rvec_in, using fmm2d
%   - Stresslets and potential dipoles at `rimage` 
%   Multiple source contributions are combined into a total velocity field.
%
% See also:
%   getPotdip, getStresslets, 
%
% Anna Broms, April 25, 2025

if size(rvec_in,1)
    eps = 1e-10;
    
    targ = [real(rcheck)';imag(rcheck)'];
    srcinfo.sources = [real(rvec_in)'; imag(rvec_in)'];
    srcinfo.stoklet = [stok_x'; stok_y'];
    U = stfmm2d(eps, srcinfo, 0, targ, 1); % might want to change eps here
    res = [U.pottarg(1,:)';U.pottarg(2,:)']/2/pi;
else
    res = zeros(size(rcheck,1)*2,1);
end
   

% With the old FMM from Lukas Bystricky
% [ufmm,vfmm] = stokesSLPfmm(stok_x,stok_y,real(rvec_in),imag(rvec_in),real(rcheck),imag(rcheck),...
%        0,5); % 5 sets high accuracy
%res = [ufmm; vfmm];

if nargin > 4
    %Not computed with FMM... to be replaced

    u_rot = getRotlets(rot,rimage,rcheck); 
    u_stress = getStresslets(stress_x,stress_y,rimage,rcheck,real(nimage),imag(nimage));
    u_pot = getPotdip(pot_x, pot_y,rimage,rcheck);

    res = res+u_rot+u_stress+u_pot; 
end

end