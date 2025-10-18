function res = applyQmat(vel,rvec_in,rvec_out,Rinv,An,Am,Nx,Ny,Mx,Z,Y,opt)
%Applies the projection matrix needed for longrange preconditioning defined
%by %B = (R\eye(3*P))*AM'; Cmap = AN*B; Qmat = eye(round(2*P*Nc))-Cmap*G;
%The matrices AN and AM contain the projection matrices to the coarse space
%sorted as x then y. 



%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind
proj_vel = getVelocityField(rvec_in,rvec_out,vel(1:end/2),vel(end/2+1:end),[],[],[],[],[],[], []);

lambda = getCoarseSource(proj_vel,Rinv,Nx,Ny,Mx,Z,Y,opt.db,opt.P,opt.N_c,opt.a_c);
%lambda = An*Rinv*(Am'*proj_vel);

res = vel-lambda;


end