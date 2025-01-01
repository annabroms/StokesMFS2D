%test script for mobility and resistance)
clear; 
close all;
rng(3); 

%generate geometry
delta = 1; %minimum particle separation
delta = 0.01; 
P = 2000; % number of particles
P = 5000; 
Pvec = [10 30 50 80 100 200 500 2000];
reps = 10; %repeat experiment
Pvec = Pvec(1:6);
Pvec = 1000; %300; 
%Pvec = 2000;
reps = 1; 
for k = 1:length(Pvec)
    for i = 1:reps
        %% Set geometry
        P = Pvec(k); 
       
        %[q,B] = grow_cluster(P,delta,2); %harder to resolve...
       
        q = createDumbells(P,delta); %P even
% 
%         t = linspace(0,2*pi,200); 
%         figure()
%         for ii = 1:P
%             plot(real(q(ii))+cos(t),imag(q(ii))+sin(t));
%             hold on
%         end

         % q = 0; 
           %q = [0; 2+delta];

%          % Hexagonal lattice
%          x = 1+delta/2;
%          y = sqrt((2+delta)^2-(1+delta/2)^2);
%          q = [0; 2+delta; x+1i*y; x+2+delta+1i*y; 4+2*delta; x-1i*y; x+2+delta-1i*y]; 
% %         q = q(1:3);
        
        rads = ones(P,1);
        %% Set mobility data
        Fref = randn(P,2);
        Tref = randn(P,1);

        %% Solve mobility with 1-body precond
        images = 1; 
       % [UW,lambda_mob,it1,err1] = solve_2D_mob(q,Fref,Tref,rads,images);

        %% Solve mobility with pair corrections
        delta_pair = 0.35; %largest distance to be counted as a close pair
       % [UW2,lambda,it2,err2] = solve_mob_precond_images(q,Fref,Tref,rads,delta_pair);

        %% Solve with peanut compression
        N_peanut = 400; %number of collocation points on peanut. 
        %tic
        [UW3,lambda,it3,err3] = solve_mob_precond_peanut(q,Fref,Tref,rads,N_peanut,delta_pair);
       % toc
        
       %% Store stuff... 
       % iters(k,i) = it;
        %errvec(k,i) = err;
    end
end

%% If we want to determine 2-way error

U = [UW(1:3:end) UW(2:3:end)]; 
W = UW(3:3:end); 

%% Solve resistance only with 1-body precond
%Converges very slowly
images = 1; 
[F,T,err_res,it_res,ftest,lambda,rcheck_dom] = solve_2D_res(q,U,W,rads,images);

%% Solve resistance with pair-corrections
[F,T,err_im,it_im,ftest2] = solve_2D_precond_images(q,U,W,rads,delta_pair);
%larger residual than with peanut compression  
%% Solve resistance with peanut compression
[F,T,err_p,it_p,ftest_peanut,precond] = solve_precond_peanut(q,U,W,rads,400);

disp('Torque error')
max(abs(T-Tref'))/max(abs(Tref))

disp('Force error')
max(abs(Fref(:,1)+1i*Fref(:,2)-transpose(F)))/max(abs(Fref(:,1)+1i*Fref(:,2)))

