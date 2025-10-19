function res = lr_matvec_2D_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s,Rinv,Nx,Ny,Mx,Zi,Yi,opt)
%LR_MATVECSTOKESMFS applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem 

%Apply the preconditioner P, see applyPmat.m

vel = matvec_2D_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s);

res = applyPmat(vel,rvec_in,rvec_out,Rinv,Nx,Ny,Mx,Zi,Yi,opt);


end