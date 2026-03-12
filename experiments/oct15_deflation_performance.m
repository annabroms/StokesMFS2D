%Experiment with performance using long range precond
% This script is just to be able to time and find bottlenecks in code. 
close all; 
images = 0; %images not needed for well separated particles
delta = 0.1; 
rng(8); %reproducable

P = 6;
[q,B] = grow_cluster(P,delta,2);
% q = [0; (2+delta)];
% % 
% % %% Hexagonal lattice
%  x = 1+delta/2;
%  y = sqrt((2+delta)^2-(1+delta/2)^2);
%  q = [0; 2+delta; x+1i*y; x+2+delta+1i*y; 4+2*delta; x-1i*y; x+2+delta-1i*y]; 
q = hexagonal_lattice(delta, 2);
%q = [0:(2+delta):(P-1)*(2+delta)]';
 %q = q(1:3);
P = length(q); 
 
%q = 0; 
%q = [0; 2+delta]; %center coordinates
U = rand(P,2); %translational velocities 
W = rand(P,1); %angular velocities 
rads = ones(P,1);
visualise = 0; 
lr = 20; %lr = 3  %the number of modes will be lr-2. 
[FT_mu,lambda, it_res_mu, gmres_tol, err_res_mu] = solve_2D_res_lr(q,U,W,rads,images,lr,visualise);
[FT,lambda, it_res, gmres_tol, err_res] = solve_res_1B(q,U,W,rads,images,lr,visualise);
[FT2,lambda2, it_res2, gmres_tol, err_res2] = solve_res_1B(q,U,W,rads,images,0,visualise);
FT = rand(1,P*3); 
[UW,lambda_mob,it_mob,gmres_tol,err_mob] = solve_mob_1B(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rads,images, lr, visualise);
[UW2,lambda_mob,it_mob2,gmres_tol,err_mob2] = solve_mob_1B(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rads,images, 0, visualise);