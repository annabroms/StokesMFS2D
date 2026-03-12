function res = lr_matvec_mob_peanut(tau,rvec_in,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,q,UU,YY,L,pairs,UB_all,YB_all,UC_all,YC_all,Cmap,Lc_pair,Lf_pair,...
    Sinv,Zi,Yi,Vk)
%LR_MATVEC_MOB_PEANUT applies longrange preconditioning to the matrix 
% corresponding to 2-body preconditioning with peanut compression for the 2D Stokes mobility problem using a grid of
%coarse sources in rvec_in enhanced with images stored in rimage.

     %matvec_mob_pairprecond_peanut(x,rbase_in_c,rbase_in_f,rvec_in_c,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair)
vel = matvec_mob_pairprecond_peanut(tau,rbase_in_c,rbase_in_f,rvec_in,refine,rimage_vec,nimage,opt,rout,rout,q,UU,YY,L,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair);

%applyPmat_peanut_mob(vel,rvec_in,rvec_out,L,Sinv,q,Zi,Yi,rbase_in_c,...
  %  rbase_in_f,refine,rimage,nimage,opt,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair)
res = applyPmat_peanut_mob(vel,rvec_in,rout,L,Sinv,q,Zi,Yi,rbase_in_c,...
    rbase_in_f,refine,rimage_vec,nimage,opt,UU,YY,pairs,UB_all,YB_all,UC_all, YC_all,Cmap,Lc_pair,Lf_pair,Vk);




end