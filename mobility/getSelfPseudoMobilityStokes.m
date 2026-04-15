function [Uii,Yii,Lii,Lr] = getSelfPseudoMobilityStokes(P,q,rin_vec,rout_vec,rimage,pair_points,svd_opts)
%getSelfPseudoMobilityStokes computes two blocks per particle to build the MFS mobility single body pseudoinverse 
% based on the SVD factorisaiton of the single body matrix.          
%
% Syntax: [Uii,Yii,Lii] = getSelfPseudoMobilityStokes(P,q,rin_vec,rout_vec,rimage,pair_points)
%         [Uii,Yii,Lii] = getSelfPseudoMobilityStokes(P,q,rin_vec,rout_vec,rimage,pair_points,svd_opts)
%
% Input: 
% P         - Number of particles to compute pseudoinverse for 
% q         - Complex valued vector of center coordinates for the P particles
% rin_vec   - Complex valued vector of source points on everybody: body1,
%             body 2, ..., body P
% rout_vec  - Complex valued vector of collocation points on everybody:
%             body 1, body2,..., body P
% rimage    - vector containing complex valued locations for images sources
%             for body 1, body 2, ... , body P
% pair_points Matrix of size P x 2 containing the total number of image and
%             collocation locations for particle j on row j
%
% Output: 
% Uii         - Cell array of size P x 1, with each cell containing the matrix of
%             left singular vectors (only the ones we wish to keep, up to the numerical
%             rank set by a tolerance tol set inside this function
% Yii         - Cell array of size P x 1, with each cell containing the matrix formed
%             by VS^+, with S^⁺ a diagonal matrix containing
%             1/sigma for each sigma that is above tol and V the matrix
%             formed by right singular vectors.
% Lii         - Cell array of size P x 1, with each cell containing the projection matrices for a
%             single body onto the space of RBM
%
% Note: If the 1-body basis is used with 2-body precond, the 1-body basis
% from this function is determined with P = 1 and rimage, nimage empty,
%
% Anna Broms Feb 13, 2026

if nargin == 0
    test_projection;
    return
end
if nargin < 7 || isempty(svd_opts)
    svd_opts = struct();
end
    
mu = 1;
visualise = 0; % Keep singular-value plotting disabled by default

N = size(rin_vec,1)/P; %Proxy sources per body. Assumed the same for everyone.

% Choose truncation level for the SVD for the single body. Appropriate
% level may depend on if images are in play for the single body compression
% (this is only the case with plain 1-body preconditioning). 
if isempty(rimage) %compute just one factor, for a single particle
    tol = 1e-10;
    tol = 1e-14; %just remove the last singval.
    %tol = 1e-11;
    %tol = 1e-16;
else
    tol = 1e-10;
    tol = 1e-9;
    tol = 1e-14; 
    %tol = 1e-16;
    %tol = 1e-9;
    %tol = 1e-6;
end

%struct to store pseudoinverses - If plain 1-body with images, there will
%be one block per particle 
Uii = cell(P,1); 
Yii = cell(P,1);
Lii = cell(P,1);



if ~size(rimage,1)
    P = 1; %Do the compression for a single particle
    noim = 1;
else
    noim = 0;
end


%We assume equal circles so that discretisation is the same for all particles
if noim
    Kin = getKmat2D(rin_vec(1:N),0);
else
    Kin = getKmat2D(rin_vec(1:N)-q(1),0);
end

L = Kin*((Kin'*Kin)\Kin'); %This is x, then y for both rows and cols
Lii{1} = L;

start_im = 0;
start_colloc = 0;

% Get self-interaction matrix for each particle
for i = 1:P
    %Get image data
    if isempty(rimage)
        rim = [];
    else
        rim = rimage(start_im+1:start_im+pair_points(i,1));
       
        start_im = start_im+pair_points(i,1);
    end

    %Get collocation points
    M = pair_points(i,2);
    rout = rout_vec(start_colloc+1:start_colloc+M);

    %...and source points
    rin = [rin_vec(N*(i-1)+1:N*i); rim];

    %Get single layer
    Nio = stokSLPmat(rin,rout,mu);

    if noim
        Kout = getKmat2D(rout,0); %assume here that we use the base discretisation
    else
        Kin = getKmat2D(rin,q(i)); 
        Kout = getKmat2D(rout,q(i));
    end

    Lr = Kout*Kin'; %to close the system with unknown velocities. Remember: we use a recompleted formulation!

    Nii = Nio*(eye(size(L))-L)+Lr; %used for mobility
      

    %Compute SVD
    svd_opts_body = svd_opts;
    if logical(getOptField(svd_opts,'left_weight',false))
        svd_opts_body.row_weights = getPeriodicCurveWeights(rout);
    end
    [Y,Bi1]  = getPseudoFactors(Nii,tol,visualise,svd_opts_body);

    %Store factors
    Uii{i} = Bi1';
    Yii{i} = Y;

    start_colloc = start_colloc+M;

end

end

function test_projection

P = 1; 
q  = 0; 
Rp = 0.6;
N = 50;
t = linspace(0,2*pi,N+1);
t = t(1:end-1)';
rin = Rp*(cos(t)+1i*sin(t));

tout = linspace(0,2*pi,1.2*N+1);
tout = tout(1:end-1)';
rout = cos(tout)+1i*sin(tout); 
pair_points = [zeros(2,1) 1.2*N*ones(2,1)];

[Uii,Yii,Lii] = getSelfPseudoMobility(P,q,rin,rout,[],[],pair_points);

lambda = [rand(1)*ones(N,1); rand(1)*ones(N,1)]; %some random force density that corresponds to a net force on the particle
L = Lii{1};
disp('Projection supposed to give zero vector:')
norm((eye(size(L))-L)*lambda) 

end
