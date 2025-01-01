function res = getFMMVelocity(rvec_in,rcheck,stok_x,stok_y,rimage,nimage,pot_x,pot_y,stress_x,stress_y)
    [ufmm,vfmm] = stokesSLPfmm(stok_x,stok_y,real(rvec_in),imag(rvec_in),real(rcheck),imag(rcheck),...
            0,5);
    res = [ufmm; vfmm];

    if nargin > 4
        %Not computed with FMM... to be replaced
        u_stress = getStresslets(stress_x,stress_y,rimage,rcheck,real(nimage),imag(nimage));
        u_pot = getPotdip(pot_x, pot_y,rimage,rcheck);
    
        res = res+u_stress+u_pot; 
    end

end