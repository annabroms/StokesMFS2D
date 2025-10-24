clear; 
close all;
delta = 0.1; 
P = 250; 
images = 0; 
q = [0:(2+delta):(P-1)*(2+delta)]'; % all particles along a line 

U = rand(P,2); %translational velocities 
W = rand(P,1); %angular velocities 
rads = ones(P,1);
visualise = 0; 
lr = 40;  %the number of modes will be lr-2. 
[FT,lambda, it_res, gmres_tol, err_res] = solve_2D_res(q,U,W,rads,images,lr,visualise);
[FT,lambda, it_res2, gmres_tol, err_res2] = solve_2D_res(q,U,W,rads,images,0,visualise);
FT = rand(1,P*3); 
[UW,lambda_mob,it_mob,gmres_tol,err_mob] = solve_2D_mob(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rads,images, lr, visualise);
[UW2,lambda_mob,it_mob2,gmres_tol,err_mob2] = solve_2D_mob(q,[FT(1:3:end)' FT(2:3:end)'],FT(3:3:end)',rads,images, 0, visualise);