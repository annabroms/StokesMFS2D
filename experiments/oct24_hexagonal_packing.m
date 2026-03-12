close all; 
images = 0; %images not needed for well-separated particles
delta = 0.1; 
rng(8); %reproducable

% %% Hexagonal lattice
q = hexagonal_lattice(delta, 16); %the second parameter sets radius of the packing, which determines the number of particles
P = length(q); 
U = rand(P,2); %translational velocities 
W = rand(P,1); %angular velocities 
rad = ones(P,1);
visualise = 0; 
lr = 20; %the number of singular modes for the long range preconditioner will be lr-2. 
[FT,lambda, it_res, gmres_tol, err_res] = solve_res_1B(q,U,W,rad,images,lr,1);
[FT,lambda, it_res2, gmres_tol, err_res2] = solve_res_1B(q,U,W,rad,images,0,visualise);
FT = rand(1,P*3); 
[UW,lambda_mob,it_mob,gmres_tol,err_mob] = solve_mob_1B(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rad,images, lr, visualise);
[UW2,lambda_mob,it_mob2,gmres_tol,err_mob2] = solve_mob_1B(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rad,images, 0, visualise);