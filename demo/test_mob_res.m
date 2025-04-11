%test script for mobility and resistance with pair corrections 
clear; 
close all;
rng(8); 

%generate geometry
delta = 1; %minimum particle separation
delta = 0.01; 
P = 2000; % number of particles
P = 5000; 
Pvec = [10 30 50 80 100 200 500 2000];
reps = 10; %repeat experiment
Pvec = Pvec(1:6);
Pvec = 1000; %300; 
Pvec = 2; 
%Pvec = 101; 
%Pvec = 2000;
reps = 1; 

test_mob_first = 0; 

for k = 1:length(Pvec)
    for i = 1:reps
        %% Set geometry
        P = Pvec(k); 
       
       % [q,B] = grow_cluster(P,delta,2); %harder to resolve...
       
        %q = createDumbells(P,delta); %P even
        %q = (0:2+delta:(P-1)*(2+delta))';
% 
%         t = linspace(0,2*pi,200); 
%         figure()
%         for ii = 1:P
%             plot(real(q(ii))+cos(t),imag(q(ii))+sin(t));
%             hold on
%         end

         % q = 0; 
         q = [0; (2+delta)];
% 
% %          % Hexagonal lattice
         x = 1+delta/2;
         y = sqrt((2+delta)^2-(1+delta/2)^2);
         %q = [0; 2+delta; x+1i*y; x+2+delta+1i*y; 4+2*delta; x-1i*y; x+2+delta-1i*y]; 
         %q = q(1:3);

        % q = [q; 30+q1; -30+q1];
        % q = [q; 10+q; -10+q];

        % q = [0; 2+delta; -2-delta];
        % q = [0; 2+delta];
         %q = [0; 2+delta; (2+delta)*1i; 2+delta+(2+delta)*1i];
         %q = q(1:3); 
        % q = [0; 2+delta; (2+delta)*1i];
        %q = [q; 10; 12+delta; -10; -12-delta];



         
%         
        rads = ones(P,1);
        delta_pair = 3; %largest distance to be counted as a close pair

        if test_mob_first
            %% Set mobility data
            Fref = randn(P,2);
            Tref = randn(P,1);
            %Tref = zeros(size(Tref));
    
    
            %% Solve mobility with 1-body precond
            images = 1; 
            visualise = 1; 
            [UW,lambda_mob,it1,gmres_tol,err1] = solve_2D_mob(q,Fref,Tref,rads,images, visualise);
    
            %% Solve mobility with pair corrections
          
            [UW2,lambda,it2,gmres_tol,err2] = solve_mob_precond_images(q,Fref,Tref,rads,delta_pair,visualise);
    
            %% Solve with peanut compression
            N_peanut = 400; %number of collocation points on peanut. 
            %tic
            [UW3,lambda,it3,err3] = solve_mob_precond_peanut(q,Fref,Tref,rads,N_peanut,delta_pair,visualise);
           % toc
        else
            %% Set resistance data instead
            Uref = randn(P,2);
            Wref = randn(P,1);
            images = 1; 
            visualise = 0;
            [FT,lambda,it, gmres_tol, maxres] = solve_2D_res(q,Uref,Wref,rads,images,visualise);
            %% Solve resistance with pair-corrections
            [FT,lambda,it_im,tol_im,err_im] = solve_res_precond_images(q,Uref,Wref,rads,delta_pair,visualise);
            %larger residual than with peanut compression  
            %% Solve resistance with peanut compression
            [FT,lambda,it_p,tol_p,err_p] = solve_res_precond_peanut(q,Uref,Wref,rads,400,delta_pair,visualise);
        end
            
        
       %% Store stuff... 
       % iters(k,i) = it;
        %errvec(k,i) = err;
    end
end

if test_mob_first

    %% If we want to determine 2-way error
    
    U = [UW3(1:3:end) UW3(2:3:end)]; 
    W = UW3(3:3:end); 
    
    
    %U = rand(P,2);
    %W = rand(P,1); 
    %% Solve resistance only with 1-body precond
    %Converges very slowly
    %images = 1; 
    [F,T,err_res,it_res,ftest,lambda,rcheck_dom] = solve_2D_res(q,U,W,rads,images,visualise);
    
    %% Solve resistance with pair-corrections
   % [F,T,err_im,it_im,ftest2] = solve_res_precond_images(q,U,W,rads,delta_pair);
    %larger residual than with peanut compression  
    %% Solve resistance with peanut compression
    [F2,T2,err_p,it_p,ftest_peanut,precond] = solve_precond_peanut(q,U,W,rads,400,visualise);
    
    disp('Torque error')
    max(abs(T-Tref'))/max(abs(Tref))
    
    disp('Force error')
    max(abs(Fref(:,1)+1i*Fref(:,2)-transpose(F)))/max(abs(Fref(:,1)+1i*Fref(:,2)))
    
    disp('Torque error')
    max(abs(T2-Tref'))/max(abs(Tref))
    
    disp('Force error')
    max(abs(Fref(:,1)+1i*Fref(:,2)-transpose(F2)))/max(abs(Fref(:,1)+1i*Fref(:,2)))
else
    F = [FT(1:3:end),FT(2:3:end)];
    T = FT(3:3:end);
    [UW2,lambda,it2,gmres_mob,err2] = solve_mob_precond_images(q,F,T,rads,delta_pair,visualise);
    U2 = [UW2(1:3:end) UW2(2:3:end)];
    W2 = [UW2(3:3:end)];
    
    %% Solve with peanut compression
    N_peanut = 400; %number of collocation points on peanut. 
    %tic
    [UW3,lambda,it3,gmres_p,err3] = solve_mob_precond_peanut(q,F,T,rads,N_peanut,delta_pair,visualise);
    U3 = [UW3(1:3:end) UW3(2:3:end)];
    W3 = [UW3(3:3:end)];

    
    

%     disp('Ang error')
%     max(abs(T-T'))/max(abs(Tref))
%     
%     disp('Trans error')
%     max(abs(Fref(:,1)+1i*Fref(:,2)-transpose(F)))/max(abs(Fref(:,1)+1i*Fref(:,2)))
%     
%     disp('Ang error')
%     max(abs(T2-Tref'))/max(abs(Tref))
%     
%     disp('Trans error')
%     max(abs(Fref(:,1)+1i*Fref(:,2)-transpose(F2)))/max(abs(Fref(:,1)+1i*Fref(:,2)))
end



