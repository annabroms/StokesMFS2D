function [rout_fine_other,tout] = getFineOther(q,refine,i,p2,a,N_f)
%GETFINEOTHER Build the fine target grid on body p2 for the pair (i,p2).
%
% Syntax:
%   [rout_fine_other,tout] = getFineOther(q,refine,i,p2,a,N_f)
%
% Inputs:
%   q      - particle centers
%   refine - pair-local clustered boundary nodes or their angular values
%   i,p2   - pair indices; returns targets on body p2 associated with body i
%   a,N_f  - fine-grid upsampling parameters
%
% Notes:
%   This helper is currently used only in the Stokes circular-body code
%   paths, where the body radius is normalized to 1.

rad = 1;

nout = ceil(a*N_f); % number of uniformly distributed collocation points on fine grid.
t = linspace(0,2*pi,nout+1);
t = t(1:end-1)';
rout_base = rad*(cos(t)+1i*sin(t));

t2 = refine{p2,i}; % clustered points or their angular parameterisation
if ~isreal(t2)
    fine_2 = t2(:);
    t2 = angle(fine_2-q(p2));
else
    fine_2 = q(p2)+rad*(cos(t2)+1i*sin(t2));
end

rout_fine_other = [q(p2)+rout_base; fine_2];

tout = [t; t2];


end
