function [Uii,Yii,Lii] = getSelfPseudoMobility(rin_vec,rout_vec,q,rimage,nimage,P,pairs,s,ver1)
%Get pseudoinverse matrices: rin: source points, rout: collocation points,
%rimage, P: number of particles, pairs contains data on the number of
%image and collocation points per particle
mu = 1; 
visualise = 1; 

N = size(rin_vec,1)/P;

if isempty(rimage) %compute just one factor, for a single particle
    tol = 1e-10;
    tol = 1e-14; %just remove the last singval. 
   % tol = 1e-11; 
else 
    tol = 1e-10;
    tol = 1e-9;
    tol = 1e-14; 
    %tol = 1e-9; 
    %tol = 1e-6;     
end

%struct to store pseudoinverses
Uii = cell(P,1); %One block for each particle and an additional block for the pair. TO be modified if for more pairs
Yii = cell(P,1); 
Lii = cell(P,1);



if ver1
    if ~size(rimage,1)
        P = 1; 
        noim = 1; 
    else 
        noim = 0;
    end
end

%We assume circles so that discretisation is the same for all particles
if noim
    Kin = getKmat2D(rin_vec(1:N),0);
else
    Kin = getKmat2D(rin_vec(1:N)-q(1),0);
end

L = Kin*((Kin'*Kin)\Kin'); %This is x y
Lii{1} = L; 

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

    %Get single layer
   % if solve_xy
     Nio = singleLayer(rin,rout,mu); 
%     else
%         Nio = singleLayerXY(rin,rout,mu);
%     end
    
    %Get image kernel
    if isempty(nim)
        Nimage = [];
    else
        Nimage = getImageKernels2D(rim,nim,rout,mu,s); %get the right singularities
    end
  
    if noim
        Kout = getKmat2D(rout,0);
    else
        Kout = getKmat2D(rout,q(i));
    end
  
    Lr = Kout*Kin';
    %Lr = zeros(size(Lr)); 

    if ver1

        %Nii = [Nio*(eye(size(L))-L) Nimage];
        Nii = [Nio*(eye(size(L))-L)+Lr Nimage];
        %Nii = [Nio*(eye(size(L))-L) Nimage];
    else
        
        Kin = getKmat2D([rin; rim;rim],q(i));
        L = Kin*((Kin'*Kin)\Kin');
        Lr = Kout*Kin';
        %Kim = getKmat2D(rim,q(i));
        Lii{i} = L; %This is x y
        %Need to be more careful here with the ordering of the unknowns... 
        %Nii = [Nio Nimage]*(eye(size(L))-L)+Lr;

        %x, then y
        Nii = [Nio(:,1:end/2) Nimage(:,1:pairs(i,1)) Nimage(:,2*pairs(i,1)+1:3*pairs(i,1)) ... 
            Nio(:,end/2+1:end) Nimage(:,pairs(i,1)+1:2*pairs(i,1)) Nimage(:,3*pairs(i,1)+1:end)];
        Nii = Nii*(eye(size(L))-L)+Lr;
    end

    %Compute SVD
    [Y,Bi1]  = getSVDFactors(Nii,tol,visualise);

    %Store factors
    Uii{i} = Bi1';
    Yii{i} = Y; 

    start_colloc = start_colloc+M; 

end



end



