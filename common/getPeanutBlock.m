function [DC,Y] = getPeanutBlock(rin_pair_c,rin_pair_f,rout_peanut,nimage,rimage,s,A,Lc)

mu = 1; 

%Considered for all pairs, also those not involving images.  
if size(rimage,1)
    Nimage = getImageKernels2D(rimage,nimage,rout_peanut,mu,s);
else
    Nimage =  [];
end
Nf = singleLayer(rin_pair_f,rout_peanut,mu);


if ~isempty(A)
    %Ntot = [Nf*A+Lr Nimage]; %what is Lr for the peanut??
    Ntot = [Nf*A Nimage];
    project = 1; 
else
    Ntot = [Nf Nimage]; %matrix in the rhs representing the fine grid
    project = 0; 
end

% Get peanut matrix

Npeanut = singleLayer(rin_pair_c,rout_peanut,mu); %matrix in lhs representing coarse grid evaluated at peanut boundary

if project
    Ac = getLpair(Lc);
    %Lrc = getLrPair(Bc,Bc,Kc); %this is not the correct Lrc. 
    %Npeanut = Npeanut*Ac+Lrc;
    Npeanut = Npeanut*Ac;
end

[UU,S,V] = svd(Npeanut); 
S = diag(S);

%tol = 1e-6;
tol = 1e-12;
tol = 1e-14;
tol = eps; 

%tol = 1e-8; 

%use relative tolerance 
ra = sum(S>max(S)*tol); 
S = S(1:ra);  
iS = 1./S; % rank
Y = V(:,1:ra)*diag(iS); 
DC = UU(:,1:ra)'*Ntot;


end