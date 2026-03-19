function [tau_stokes_x, tau_stokes_nonpx,tau_self_x, tau_beta_x,tau_stokes_y,tau_stokes_nonpy,tau_self_y,tau_beta_y,u_corr] = transform_mob_peanut(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,rcheck_out,q,U,Y,Lc,pairs,Upf,Ypf,DC_all, YC_all,Cmap,Lc_pair,Lf_pair)
%Transform coarse \mu -> coarse \lambda using peanut compression, to be
%used in mobility matvec. 
                                                                              
P = length(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P; 
PM = length(rvec_out);
N_check = length(rcheck_out)/P;
PM2 = length(rcheck_out);

%s = opt.s; %[0 0 1 1]; %stresslet and potential dipoles at image points

%map densities back

%Preallocate
tau_stokes_x = zeros(N_c*P,1); 
tau_stokes_y = zeros(N_c*P,1);
tau_stokes_nonpx = zeros(N_c*P,1);
tau_stokes_nonpy = zeros(N_c*P,1);

tau_beta_x = zeros(N_f*P,1);
tau_beta_y = zeros(N_f*P,1);

rimage_k = cell(P,1);
nimage_k = cell(P,1);

u_corr = zeros(2*N_check*P,1); 
precomp = opt.precomp;  


for i = 1:P

    coarse_ind = (i-1)*N_c+1:i*N_c;  

    %Retrieve self evaluation blocks    
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %U{i}, Y{i}
    step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other?
    tau_mapped = Y{1}*step1; %this is the mapped density for this particle to throw in to the kernel

    tau_stokes_nonpx(coarse_ind) = tau_mapped(1:N_c);
    tau_stokes_nonpy(coarse_ind) = tau_mapped(N_c+1:end);

    %Project
    tau_i_x = tau_mapped(1:N_c);
    tau_i_y = tau_mapped(N_c+1:end);
    tau_mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y]; 

    tau_stokes_x(coarse_ind) = tau_mapped(1:N_c);
    tau_stokes_y(coarse_ind) = tau_mapped(N_c+1:end);

end


%For subtracting off self-contribution in matvec
tau_self_x = tau_stokes_x;
tau_self_y = tau_stokes_y; 

s_ind1_x = 1:opt.N_c; 
s_ind2_x = opt.N_c+1:2*opt.N_c;
s_ind1_y = 2*opt.N_c+1:3*opt.N_c;
s_ind2_y = 3*opt.N_c+1:4*opt.N_c;


for i = 1:P

    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %check if particle is in pair.
    if ~isempty(pairs)
        %neigh = find(pairs(:,1)==i);
        [I,J] = find(pairs==i);
        neigh = [I J];
    else
        neigh = [];
    end
    
    if ~isempty(neigh) %if isempty(neigh) has no neighbours 
        
        for k = 1:size(neigh,1)



            if neigh(k,2) == 2 %has already been counted in 
                break
            end
            %send in two particles with the fine grid 
            p2 = pairs(neigh(k),2); 

            rimage_k{i} = [rimage_k{i}; rimage_vec{i,p2}];
            rimage_k{p2} = [rimage_k{p2}; rimage_vec{p2,i}];
            nimage_k{i} = [nimage_k{i}; nimage{i,p2}];
            nimage_k{p2} = [nimage_k{p2}; nimage{p2,i}];

            tau_particle_x2 = tau((p2-1)*N_large+1:N_large*p2);
            tau_particle_y2 = tau(PM+(p2-1)*N_large+1:PM+N_large*p2);

            %U{i},Y{i}
            step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
            tau_mapped = Y{1}*step1;
            tau_i_x = tau_mapped(1:N_c);
            tau_i_y = tau_mapped(N_c+1:end);
            tau_mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y];

            %U{p2}, Y{p2}
            step1 = U{1}*[tau_particle_x2;tau_particle_y2]; 
            mapped = Y{1}*step1;
            tau_i_x = mapped(1:N_c);
            tau_i_y = mapped(N_c+1:end);
            mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y];

            %% First, determine beta, the fine sources
            if ~precomp                
                
                rout_fine_other = getFineOther(q,refine,i,p2,opt.a_f,opt.N_f); 
                %Nother = stokSLPmat(rbase_in_c+q(i),rout_fine_other,mu);
                %R2 = -Nother*tau_mapped; %read off on particle 2
    
                [u2,v2] = stokSLPdirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    tau_mapped(1:N_c),tau_mapped(N_c+1:2*N_c),N_c);
                R2  = -[u2; v2];
    
                %% Do a similar thing for the other order of the particles in the pair
                rout_fine_other = getFineOther(q,refine,p2,i,opt.a_f,opt.N_f); 
                %Nother2 = stokSLPmat(rbase_in_c+q(p2),rout_fine_other,mu);
                %R1 = -Nother2*mapped; %read off on particle 1
    
                %To be replaced with C implementation?
                [u1,v1] = stokSLPdirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    mapped(1:N_c),mapped(N_c+1:2*N_c),N_c);
                R1 = -[u1; v1];

            
                rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)]; 

            else              
                rhs = [tau_mapped(1:end/2); mapped(1:end/2); tau_mapped(end/2+1:end); mapped(end/2+1:end)];
            end
   
            pair_mapped = Upf{i,p2}*rhs; 
            tau_mapped_tot = Ypf{i,p2}*pair_mapped;

            %now we have beta. Lets do peanut compression from
            %here
            
            if opt.cmap
                tau_peanut_ntot = Cmap{i,p2}*rhs;
            else
                tau_peanut_temp = DC_all{i,p2}*tau_mapped_tot;
                tau_peanut_ntot = YC_all{i,p2}*tau_peanut_temp;
            end

            % Get projected equivalent sources
            tau_peanut_tot= Lc_pair*tau_peanut_ntot; 

            %Store
            tau_stokes_x((i-1)*N_c+1:N_c*i) = tau_stokes_x((i-1)*N_c+1:N_c*i)+...
                    tau_peanut_tot(s_ind1_x);
            tau_stokes_y((i-1)*N_c+1:N_c*i) = tau_stokes_y((i-1)*N_c+1:N_c*i)+...
                tau_peanut_tot(s_ind1_y);

            tau_stokes_x((p2-1)*N_c+1:N_c*p2) = tau_stokes_x((p2-1)*N_c+1:N_c*p2)+...
                tau_peanut_tot(s_ind2_x);
            tau_stokes_y((p2-1)*N_c+1:N_c*p2) = tau_stokes_y((p2-1)*N_c+1:N_c*p2)+...
                tau_peanut_tot(s_ind2_y);

            %Need to store to close system (impose bc)
            tau_beta_x((i-1)*N_f+1:N_f*i) = tau_beta_x((i-1)*N_f+1:N_f*i)+tau_mapped_tot(1:N_f);
            tau_beta_x((p2-1)*N_f+1:N_f*p2) = tau_beta_x((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(N_f+1:2*N_f);
            tau_beta_y((i-1)*N_f+1:N_f*i) = tau_beta_y((i-1)*N_f+1:N_f*i)+tau_mapped_tot(2*N_f+1:3*N_f);
            tau_beta_y((p2-1)*N_f+1:N_f*p2) = tau_beta_y((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(3*N_f+1:4*N_f);

            %% Evaluate flow field on pair and subtract this contribution. 
            % Replace with fine grid on the pair
            rout_pair = [rcheck_out((i-1)*N_check+1:i*N_check,:); rcheck_out((p2-1)*N_check+1:p2*N_check,:)];
                       
            rimage = [rimage_vec{i,p2}; rimage_vec{p2,i}];
            nimage_pair = [nimage{i,p2}; nimage{p2,i}]; 
            rim = length(rimage);
            
            %get contribution from image points 
            u_stress = getStresslets(tau_mapped_tot(4*N_f+1:4*N_f+rim),...
                tau_mapped_tot(4*N_f+1+rim:4*N_f+2*rim),rimage,...
                rout_pair,real(nimage_pair),imag(nimage_pair));

            u_pot = getPotdip(tau_mapped_tot(4*N_f+1+2*rim:4*N_f+3*rim),...
                tau_mapped_tot(4*N_f+1+3*rim:4*N_f+4*rim),rimage,rout_pair);

            %... and from fine grid of Stokeslets 
            rin_pair = [rbase_in_f+q(i); rbase_in_f+q(p2)];

            tau_mapped_proj = Lf_pair*tau_mapped_tot(1:4*N_f); 


            [u1,v1] = stokSLPdirect(real(rin_pair),imag(rin_pair),...
                real(rout_pair),imag(rout_pair),...
                tau_mapped_proj(1:2*N_f),tau_mapped_proj(2*N_f+1:4*N_f),2*N_f);
            u_stok = [u1; v1];
            u_pair = u_stress+u_pot+u_stok;

            %Determine contribution to be subtracted

            rin_pair = [rbase_in_c+q(i); rbase_in_c+q(p2)];               

            [u1,v1] = stokSLPdirect(real(rin_pair),imag(rin_pair),...
                real(rout_pair),imag(rout_pair),...
                tau_peanut_tot(1:2*N_c),tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
            u_peanut_corr = [u1; v1];
 
            pair_ind = [(i-1)*N_check+1:i*N_check (p2-1)*N_check+1:p2*N_check ...
                (i-1)*N_check+PM2+1:i*N_check+PM2 ...
                (p2-1)*N_check+PM2+1:p2*N_check+PM2]'; 
            u_corr(pair_ind) = u_corr(pair_ind)+u_pair-u_peanut_corr;


            
        end
        


    end
       
end

end
