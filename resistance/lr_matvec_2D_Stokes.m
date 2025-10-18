function res = lr_matvec_2D_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s,Rinv,Nx,Ny,Mx,Z,Y,opt)
%LR_MATVEC_2D_STOKES applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem using a grid of
%coarse sources in rvec_in enhanced with images stored in rimage.


vel = matvec_2D_Stokes(tau,rvec_in,rvec_out,rimage,nimage,q,Uii,Yii,pairs,s);

%Apply the preconditioner P defined as B = Rinv*AM';
    %B = (Rc\eye(3*P))*AM';%us coarse R instead
    %Cmap = AN*B;
   % P = eye(round(2*P*a*Nc))-G*Cmap;

res = applyPmat(vel,rvec_in,rvec_out,Rinv,Nx,Ny,Mx,Z,Y,opt);



end