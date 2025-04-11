function [DC, Y] = getPeanutBlock(rin_pair_c, rin_pair_f, rout_peanut, nimage, rimage, s, Lc_pair, Lf_pair)
%GETPEANUTBLOCK Computes factorization for peanut compression for two circular particles
% in Stokes flow using the method of fundamental solutions (MFS)
%
% Syntax:
%   [DC, Y] = getPeanutBlock(rin_pair_c, rin_pair_f, rout_peanut, nimage, rimage, s, Lc_pair, Lf_pair)
%
% Inputs:
%   rin_pair_c  - Complex vector of proxy source coordinates on both particles (coarse grid)
%   rin_pair_f  - Complex vector of proxy source coordinates on both particles (fine grid)
%   rout_peanut - Complex vector of target points on the peanut interface (separation surface)
%   nimage      - Complex vector encoding directions (as x + iy) for image stresslets on body 1 and 2
%   rimage      - Complex vector of image source locations for body 1 and 2
%   s           - Logical vector indicating source types at image locations [S R T D]
%   Lc_pair     - Projection matrix onto the null space of the force/torque constraint
%                 for coarse sources on both particles
%   Lf_pair     - Projection matrix onto the null space of the force/torque constraint
%                 for fine sources on both particles
%
% Outputs:
%   DC - Product of:
%          1) left singular vectors matrix from an SVD of the evaluation of the coarse-grid on the peanut interface, and
%          2) the fine-grid evaluation matrix for the pair on the peanut
%   Y  - Matrix given by V * S⁺ from the same SVD, where:
%          - S⁺ is a diagonal matrix with 1/σ for retained singular values σ
%          - V is the matrix of corresponding right singular vectors
%
% Description:
%   This function enables a backward-stable compression mapping from fine source strengths
%   (β) to effective coarse proxy strengths (λ_coarse_effective), via:
%
%       λ_coarse_effective = Y * (DC * β)
%
%   The factorization is used to accelerate and stabilize close interaction corrections between
%   two near-touching particles in Stokes flow.
%
% Notes:
%   - Only singular values above an internal tolerance are retained
%   - Designed for use in 2-body preconditioning with peanut compression
%
% See also:
%   getPseudoFactors
%
% Anna Broms, April 4, 2025

mu = 1; 

%% Get fine representation: fine proxy grid of Stokeslets + image
% enhancement, evaluated on peanut boundary 
if size(rimage,1)
    Nimage = getImageKernels2D(rimage,nimage,rout_peanut,mu,s);
else
    Nimage =  [];
end
Nf = singleLayer(rin_pair_f,rout_peanut,mu);

%% If mobility, need to project so that the fine sources don't contribute to force and torque
if ~isempty(Lf_pair)
    Ntot = [Nf*Lf_pair Nimage];
    project = 1; 
else
    Ntot = [Nf Nimage]; %matrix in the rhs of the matching LSQ problem representing the fine grid
    project = 0; 
end

%% Get coarse representation 
% %i.e. coarse proxy grid of Stokeslets evaluated on peanut boundary
Npeanut = singleLayer(rin_pair_c,rout_peanut,mu);

% If we solve a mobility problem, the coarse sources equivalent to the fine
% sources should not contribute to a force/torque on the particle pair.
% Project off any such contribution.
if project
    Npeanut = Npeanut*Lc_pair;
end

%tol = 1e-6;
tol = 1e-12;
tol = 1e-14;
%tol = eps; 
%tol = 1e-8; 

%Determine factors for the pseudoinverse of the coarse evaluation
[Y,U]  = getPseudoFactors(Npeanut,tol);
%Get mapping to determine (coarse sources) lambda <- beta (fine sources) 
DC = U'*Ntot;

end