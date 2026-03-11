function res = lr_matvec_2D_Stokes_mu(tau,rvec_in,rvec_out,q,Uii,Yii,pairs,s,Sinv,Zi,Yi,opt)
%LR_MATVEC_2D_STOKES_MU applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem 

%Apply the preconditioner P, see applyPmat_mu.m

vel = matvec_res_Stokes(tau,rvec_in,rvec_out,[],[],q,Uii,Yii,pairs,s);

res = applyPmat_mu(vel,rvec_in,rvec_out,Sinv,q,Zi,Yi,Uii,Yii,pairs,opt);


end