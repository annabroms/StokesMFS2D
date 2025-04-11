function [U,Y] = getPairBlock(rin_pair_f,rout_f,rimage,nimage,s,Lf_pair,Lr,project)
%getPairBlock computes factorisation of the pair correction for two
%circular particles in Stokes flow solved with MFS
%
% Syntax: [U,Y] = getPairBlock(rin_pair,rout_f,nimage,rimage,s,L,Lr,project)
%
% Input: 
% rin_pair_f    - vector of complex valued coordinates for the fine grid of
%                 proxy sources on both particles in the pair: body1, body2
% rout_f        - vector of complex valued coordinates on the boundaries
%                 discretised with a fine grid: body1, body2
% rimage        - vector of complex valued locations for images sources
%                 for body 1, body 2
% nimage        - vector containing x y "directions" stored as x+iy (generated at random) for image stresslets for 
%                 body 1, body 2
% s             - Boolean vector describing combination of source types at image
%                 locations [S R T D]
% Lf_pair       - Projection matrix onto the null space of the force
%                 /torque constraint matrix using the fine grid on the two particles
% Lc            - Projection matrix onto the null space of the force
%                 /torque constraint matrix using the coarse grid on one
%                 particle.
%
% Output: 
% DC            - Matrix representing the product of two matrices: the
%                 matrix of left singular vectors for the coarse
%                 grid on the peanut, and the fine evaluation matrix for
%                 the pair.
% Y             - Matrix formed by VS^+, with S^⁺ a diagonal matrix
%                 containing 1/sigma for each singular value sigma of the matrix
%                 representing coarse sources evaluatated on the peanut
%                 boundary. Only singular values above an internal
%                 tolerance tol are kept and V is the matrix formed by
%                 right singular vectors.
% 
% Note: DC and Y builds the mapping lambda_coarse_effective <- beta backward 
% stably with beta the fine source strengths, as lambda_coarse_effective = Y*(DC*beta).
%
% Anna Broms April 4, 2025

mu = 1; 

%Visualise grid for the pair
% figure(8)
% clf; 
% plot(real(rin_pair),imag(rin_pair),'r+')
% hold on
% plot(real(rout_f),imag(rout_f),'bo')
% plot(real(rimage),imag(rimage),'k*')


if size(rimage,1)
    Nimage = getImageKernels2D(rimage,nimage,rout_f,mu,s);
else
    Nimage = [];
end
N = singleLayer(rin_pair_f,rout_f,mu); 

% In a mobility problem, need to project so that fine sources don't
% contribute to net force/torque on the particles
if project
    N = [N*Lf_pair+Lr Nimage];
else
    N = [N Nimage]; 
end 


%tol = 1e-6;
%tol = 1e-14; %Good level of truncation it seems, for close contacts. Tried with a line of circles.
tol = 1e-11; %4 nov: seems to work well but results in peaks in close to touching points
tol = 1e-14;
%tol = eps; % I don't think this is reasonable. Just to understand what happens now. 

[Y,U] = getPseudoFactors(N,tol,0); 


end