function [Uii,Yii] = getSelfPseudo(P,rin_vec,rout_vec,rimage,nimage,pair_points,s)
%getSelfPseudo computes two blocks per 2D particle to build the MFS single body pseudoinverse 
% based on the SVD factorisaiton of the single body Stokes target-from-source matrix.          
%
% Syntax: [Uii,Yii] = getSelfPseudo(P,q,rin_vec,rout_vec,rimage,nimage,pair_points,s)
%
% Input: 
% P         - Number of particles to compute pseudoinverse for 
% rin_vec   - Complex valued vector of source points on everybody: body1,
%             body 2, ..., body P
% rout_vec  - Complex valued vector of collocation points on everybody:
%             body 1, body2,..., body P
% rimage    - vector containing complex valued locations for images sources
%             for body 1, body 2, ... , body P
% nimage    - vector containing x y "directions" stored as x+iy (generated at random) for image stresslets for 
%             body 1, body 2, ..., body P
% pair_points Matrix of size P x 2 containing the total number of image and
%             collocation locations for particle j on row j
% s         - Boolean vector describing the singularity types active at image
%             points [S R T D]
%
% Output: 
% Uii         - Cell array of size P x 1, with each cell containing the matrix of
%             left singular vectors (only the ones we wish to keep, up to the numerical
%             rank set by a tolerance tol set inside this function
% Yii         - Cell array of size P x 1, with each cell containing the matrix formed
%             by VS^+, with S^⁺ a diagonal matrix containing
%             1/sigma for each sigma that is above tol and V the matrix
%             formed by right singular vectors.
%
% Note: If the 1-body basis is used with 2-body precond, the 1-body basis
% from this function is determined with P = 1 and rimage, nimage empty
%
% See also: getSelfPseudoMobility
%
% Anna Broms April 4, 2025

if nargin < 4
    rimage = [];
    nimage = [];
    P = 1; 
    pair_points = [0,length(rout_vec)]; %no images
    s = 0; 
    solve_xy = 1; 
end

solve_xy = 1; %order unknowns as x, then y globally

mu = 1; 
visualise = 1; 

N = size(rin_vec,1)/P;

if isempty(rimage) %compute just one factor, for a single particle
    tol = 1e-10;
    tol = 1e-12;
   % tol = 1e-14; %best performance if this is reasonably small, but should
   % not be too small as it leads to slow fourier decay on close neighbour
   % and large MFS coefficients
    %tol = 1e-8; 
    tol = 1e-10; 
    %should not choose eps here as it might be that we then don't remove the last singval 
    tol = 1e-14; 
else  
    tol = 1e-9;
    tol = 1e-14; 
    tol = 1e-9;
    tol = 1e-12; 
   % tol = 1e-12;
end

%struct to store pseudoinverses
Uii = cell(P,1); %One block for each particle
Yii = cell(P,1); 

start_im = 0; 
start_colloc = 0; 
for i = 1:P
    %Get image data
    if isempty(rimage)
        rim = [];
        nim = [];
    else
        rim = rimage(start_im+1:start_im+pair_points(i,1));
        nim = nimage(start_im+1:start_im+pair_points(i,1));
        start_im = start_im+pair_points(i,1); 
    end

    %Get collocation points
    M = pair_points(i,2); 
    rout = rout_vec(start_colloc+1:start_colloc+M);

    %...and source points
    rin = rin_vec(N*(i-1)+1:N*i); 

    if s(1)
        rin = [rin; rim];
    end

    %Get single layer
    if solve_xy
        Nio = singleLayer(rin,rout,mu); 
    else
        Nio = singleLayerXY(rin,rout,mu); %Different ordering of unknowns. Not in use. 
    end
    
    %Get image kernel
    if isempty(nim)
        Nimage = [];
    else
        Nimage = getImageKernels2D(rim,nim,rout,mu,s); %get the right singularities
    end
    Nii = [Nio Nimage];

    %Compute SVD
    [Y,Bi]  = getPseudoFactors(Nii,tol,visualise);

    %Store factors
    Uii{i} = Bi';
    Yii{i} = Y; 

    start_colloc = start_colloc+M; 

end

end



