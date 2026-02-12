function res = lr_matvec_2D_pairprecond(tau,rvec_in,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s,Sinv,Zi,Yi)
%LR_MATVEC_2D_PAIRPRECOND applies longrange preconditioning to the matrix 
% corresponding to one body preconditioning for the 2D Stokes resistance problem using a grid of
%coarse sources in rvec_in enhanced with images stored in rimage.


vel = matvec_2D_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);

%Apply the preconditioner P defined as B = Rinv*AM';
    %B = (Rc\eye(3*P))*AM';%us coarse R instead
    %Cmap = AN*B;
   % P = eye(round(2*P*a*Nc))-G*Cmap;

res = applyPmat(vel,rvec_in,rout,Sinv,[],[],[],Zi,Yi,opt);



end