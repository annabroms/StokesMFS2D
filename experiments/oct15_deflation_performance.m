%Experiment with performance using long range precond
close all; 
images = 0; %images not needed for well separated particles
delta = 1; 
rng(8); 

P = 200;
[q,B] = grow_cluster(P,delta,2);
%q = 0; 
%q = [0; 2+delta]; %center coordinates
U = rand(P,2); %translational velocities 
W = rand(P,1); %angular velocities 
rads = ones(P,1);
visualise = 0; 
lr = 50; %lr = 3  

[FT,lambda, it, gmres_tol, err] = solve_2D_res(q,U,W,rads,images,lr,visualise);