function res = lr_matvec_2D_mobility(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,L,Lr,pairs,s,Rinv,Nx,Ny,Mx,Zi,Yi,opt)
%LR_MATVEC_2D_MOBILITY applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem 

%Apply the preconditioner P, see applyPmat.m
project_proxy = 1; 
vel = matvec_2D_mobility(tau,rvec_in,rvec_out,rvec_out,rimage,nimage,q,Uii,Yii,L,pairs,s,1,project_proxy);


res = applyPmat_mob(vel,rvec_in,rvec_out,L{1},Lr,Rinv,Nx,Ny,Mx,Zi,Yi,opt);


end