function res = applyPmat(vel,rvec_in,rvec_out,Rinv,Nx,Ny,Mx,Z,Y,opt)
%Applies the projection matrix needed for longrange preconditioning defined
%by the sequence of matrices %B = (R\eye(3*P))*AM'; Cmap = AN*B; P = eye(round(2*P*a*Nc))-G*Cmap;
%The matrices AN and AM contain the projection matrices to the coarse space
%sorted as x then y. 

%lambda = An*Rinv*(Am'*vel);

lambda = getCoarseSource(vel,Rinv,Nx,Ny,Mx,Z,Y,opt.db,opt.P,opt.N_c,opt.a_c);

%compute velocities using these source strengths

%remember that everything so far is implemented without images in mind
proj_vel = getVelocityField(rvec_in,rvec_out,lambda(1:end/2),lambda(end/2+1:end),[],[],[],[],[],[], []);

res = vel-proj_vel;


end