function [min_eig,c,no_image] = get_eigs(q,delta_pair)

rad = ones(size(q));

%% SET PARAMS
%GMRES params
maxit = 800; 
gmres_tol = 1e-6; 

gmres_tol = 1e-10; 

%Grid params
P = length(q); 
%Set coarse and fine grid. 
%Play with N_c, N_f, a (a_f). 
N_c = 80;  %100 better here? 
%N_c = 100; 
%N_c = 100; 
%N_c = 100; 
N_f = 150; 
N_f = 1000; 
%N_f = N_c; %debug
%N_c = 2; 
%N_c = 150; 

%N_c = 250; 
a_c = 1.2;
%a_c = 2; 
%a_c = 5;


%a_c = 4;
%a_c = 1.2; %results in horrible fourier decay if we 
           % evaluate on a particle close by in coarse colloc points 
%a = 2; 
a_f = 1.2; %upsampling factor for the fine grid
%a_f = 5; 
%a2 = 0.6; 
%a2 = a; 

tol_c = 1e-12; %I think this works reasonably
%tol_c = 1e-8; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 
%tol_c = 1e-11; 
% But smoother coarse 1-body basis to evaluate on neighbour


s = [0 0 1 1]; %set type of singularities at image points
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rg_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rg_f = max([1-sep_f,0.01]);  % and fine grid
%Rg_f = Rg_c; %debug


%Rg_c = 0.2; %debug! 
%Rg_c = 0.01; 

%To be sent to solver/matvec later
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.rad = rad; 
opt.s = s; 
opt.precomp = 0; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.


%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';

rbase_out_c = rad(1)*cos(tout_c)+1i*rad(1)*sin(tout_c);

%Construct image grid
basic = 1; %return only the basic outer grid.                               
[rout,~,~,~,~,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rad,Rg_c, a_c, N_c, 1,Rg_f,a_f,N_f,basic,delta_pair);

% %Now, very deliberately, remove the 1-3 pair!
% pairs = [1 2; 2 3]; 
% rimage_vec{1,3} = [];
% rimage_vec{3,1} = [];
% refine{1,3} = [];
% refine{3,1} = []; 


%% PREPARE PRECONDITIONING AND RHS
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rg_c*cos(tin)+Rg_c*1i*sin(tin);
rin = [];
for k = 1:P
    rin = [rin; q(k)+rbase_in_c];
end

rimage_in = []; 
nimage = [];
[U,Y,Lc] = getSelfPseudoMobility(rbase_in_c,rbase_out_c,q,rimage_in,nimage,1,[0,ceil(a_c*N_c)],s,1);

%K = getKmat2D(rbase_in_c,0);

Kf = getKmat2D(rbase_in_f,0);
Lf = Kf*((Kf'*Kf)\Kf'); %This is x y

%Get pair basis

[Upf,Ypf,~,~,nimage] = getPairBasis(q,a_f,N_f,rad,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,s,opt,0,Lf,Kf,1);


% Now, check pair basis up to the boundary. Is it nice and smooth?
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,Lc{1},Lf,N_c,N_f,a_c,a_f,rad)

%% SOLVE SYSTEM
% Build the matrix to check it out


x = zeros(2*size(rout,1),1);
for k = 1:2*size(rout)
    k
    x(:) = 0; 
    x(k) = 1; 
    uu = matvec_mob_pairprecond_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf);
    CC(:,k) = uu;
end
[V,D] = eig(CC);
D = diag(D); 
min_eig = min(D);
c = skeel(CC); 

figure()
imagesc(log10(abs(CC)))
colorbar

no_image =  isempty(rimage_vec{1,3});

debug = 0;
if debug


 
    figure()
    plot(real(D),imag(D),'+')
    xlabel('Re \lambda')
    ylabel('Im \lambda')

    [s,I] = mink(abs(D),3);
    Vsmall = V(:,I).*s';

    [UU,S,VV] = svd(CC);
    SS = diag(S);
    figure()
    semilogy(SS)
    semilogy(SS,'+')
    [s2,I2] = mink(SS,6);

    UUs = UU(:,I2);
    VVs = VV(:,I2); 

    Mc = round(a_c*N_c);
    t = linspace(0,2*pi,Mc+1);
    t = t(1:end-1)'; 

    %Visualise left and right singular vectors
    for i = 1:6
        for k = 1:P
            Vpx = VVs((k-1)*Mc+1:k*Mc,i);
            Vpy = VVs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
            Vp = abs(Vpx+1i*Vpy);
            figure(33)
            subplot(2,3,i)
            scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
            hold on
            colorbar
            view(0,90)
            sgtitle('Right sing vec','interpreter','latex')
            title(i)
            axis off
            axis square

            figure(34)
            subplot(2,3,i)
            quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
            hold on
            q2 = [q; q(1)];
            plot(real(q2),imag(q2),'k-')
            axis off
            sgtitle('Right sing vec','interpreter','latex')
            title(i)
            axis square


            Vpx = UUs((k-1)*Mc+1:k*Mc,i);
            Vpy = UUs((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
            Vp = abs(Vpx+1i*Vpy);
            figure(35)
            subplot(2,3,i)
            scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
            hold on
            colorbar
            view(0,90)
            sgtitle('Left sing vec','interpreter','latex')
            title(i)
            axis off
            axis square

            figure(36)
            subplot(2,3,i)
            quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
            hold on
            q2 = [q; q(1)];
            plot(real(q2),imag(q2),'k-')
            axis off
            sgtitle('Left sing vec','interpreter','latex')
            title(i)
            axis square
        end
    end
            
    %Visualise eigenvectors

    for i = 1:3
     
        for k = 1:P  
            Vpx = Vsmall((k-1)*Mc+1:k*Mc,i);
            Vpy = Vsmall((k-1)*Mc+Mc*P+1:k*Mc+Mc*P,i);
            Vp = abs(Vpx+1i*Vpy);
            figure(31)
            subplot(1,3,i)
            scatter3(real(q(k))+cos(t),imag(q(k))+sin(t),Vp,40,Vp,'filled');
            hold on
            colorbar
            view(0,90)
            sgtitle('Abs of eigvec','interpreter','latex')
            axis off
            axis square


            figure(32)
            subplot(1,3,i)
            quiver(real(q(k))+cos(t),imag(q(k))+sin(t),Vpx,Vpy)
            hold on
            q2 = [q; q(1)];
            plot(real(q2),imag(q2),'k-')
            sgtitle('Eigvec','interpreter','latex')
            axis off
            axis square



        end

    end
    

end