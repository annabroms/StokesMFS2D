function res = lr_matvec_2D_peanut(tau,rvec_in,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all,YC_all,Cmap,Sinv,Zi,Yi)
%LR_MATVEC_2D_PEANUT applies longrange preconditioning to the matrix 
% corresponding to 2-body preconditioning with peanut compression for the 2D Stokes resistance problem using a grid of
%coarse sources in rvec_in enhanced with images stored in rimage.
debug = 0; 
vel = matvec_2D_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in,rbase_out_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,debug);

%vel = matvec_2D_peanut(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,pairs,Upf,Ypf,s);

%Apply the preconditioner P defined as B = Rinv*AM';
    %B = (Rc\eye(3*P))*AM';%us coarse R instead
    %Cmap = AN*B;
   % P = eye(round(2*P*a*Nc))-G*Cmap;
res = applyPmat_peanut(vel,rvec_in,rout,Sinv,q,Zi,Yi,rbase_in_c,...
    rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap);
%res = applyPmat_peanut(vel,rvec_in,rout,Sinv,[],[],[],Zi,Yi,opt);



end