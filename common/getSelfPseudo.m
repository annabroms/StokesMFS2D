function [Uii,Yii] = getSelfPseudo(rin_vec,rout_vec,rimage,nimage,P,pairs,solve_xy,s)
%Get pseudoinverse matrices: rin: source points, rout: collocation points,
%rimage, image points ,nimage, "normal" (random) vectors needed for stresslets,
% np: number of particles, pairs contains data on the number of
%image and collocation points per particle
if nargin < 3
    rimage = [];
    nimage = [];
    P = 1; 
    pairs = [0,length(rout_vec)];
    s = 0; 
    solve_xy = 1; 
end

mu = 1; 
visualise = 0; 

N = size(rin_vec,1)/P;

if isempty(rimage) %compute just one factor, for a single particle
%     if ~wobbly
%         np = 1;
%     end
    tol = 1e-10;
    tol = 1e-12;
    tol = 1e-14; %best performance if this is small   
    %tol = 1e-15; %Not good as it leads to slow fourier decay on close neighbour.
    %should not choose eps here as it might be that we then don't remove the last singval 
else  
    tol = 1e-9;
    tol = 1e-14; 
    tol = 1e-9;
    tol = 1e-14;   
end

%struct to store pseudoinverses
Uii = cell(P,1); %One block for each particle and an additional block for the pair. TO be modified if for more pairs
Yii = cell(P,1); 

start_im = 0; 
start_colloc = 0; 
for i = 1:P
    %Get image data
    if isempty(rimage)
        rim = [];
        nim = [];
    else
        rim = rimage(start_im+1:start_im+pairs(i,1));
        nim = nimage(start_im+1:start_im+pairs(i,1));
        start_im = start_im+pairs(i,1); 
    end

    %Get collocation points
    M = pairs(i,2); 
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
    [Y,Bi1]  = getSVDFactors(Nii,tol,visualise);

    %Store factors
    Uii{i} = Bi1';
    Yii{i} = Y; 

    start_colloc = start_colloc+M; 

end

end



