function [rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx, tau_stokes_nonpy,tau_stress_x,tau_stress_y,tau_stress_all_x,tau_stress_all_y,tau_stress_all_px,tau_stress_all_py,...
    tau_pot_x,tau_pot_y,tau_pot_all_x,tau_pot_all_y,tau_pot_all_px,tau_pot_all_py,u_corr] = getMobPairTransformation(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,nimage,opt,rvec_out,q,U,Y,Lc,Lf,pairs,Upf,Ypf)
%GETPAIRTRANSFORMATION maps data at coarse collocation nodes back to coarse and fine source
%strengths, in preparation for the mobility matvec. This is very similar to
% getPairTransformation.m for the resistance problem, with the difference
% that probjections are done. 

%Very long function call...

P = length(q);
N_coarse = opt.N_c;
N_large = length(rvec_out)/P; 
PM = length(rvec_out); 

two_parts = 0; 
precomp = opt.precomp; 
use_matrix_free_projection = true; % set false to use the original K-based projector


%map densities back
%Bookkeeping stuff

%should be preallocated...
tau_stokes_x = []; 
tau_stokes_y = []; 
tau_stokes_nonpx = [];
tau_stokes_nonpy = [];

rvec_in = []; 
rimage_in = []; 
nimage_in = []; 



coarse_ind = cell(P,1); 
tau_stokes_fine_x = cell(P,1); 
tau_stokes_fine_y = cell(P,1); 

tau_stokes_fine_nonpx = cell(P,1); 
tau_stokes_fine_nonpy = cell(P,1); 

tau_stress_x = cell(P); 
tau_stress_y = cell(P); 

tau_pot_x = cell(P); 
tau_pot_y = cell(P); 

rimage_k = cell(P,1);
nimage_k = cell(P,1);

u_corr = zeros(2*PM,1); 
%Store local fine grid correction


for i = 1:P
    %Retrieve self evaluation blocks

    %check if particle is in any pair.
    if ~isempty(pairs)
        %neigh = find(pairs(:,1)==i);
        [I,J] = find(pairs==i);
        neigh = [I J];
    else
        neigh = [];
    end
    
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %if isempty(neigh) %has no neighbours 
    step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other?
    tau_mapped = Y{1}*step1; %this is the mapped density for this particle to throw in to the kernel
    
    %check residual for the self-interaction only
    %NN = stokSLPmat(rbase_in_c+q(i),rvec_out((i-1)*N_large+1:i*N_large,:),mu);
    %disp('Self-interaction error')
    %norm(NN*tau_mapped-[tau_particle_x;tau_particle_y],inf)
    %norm((NN*tau_mapped-[tau_particle_x;tau_particle_y])./[tau_particle_x;tau_particle_y],inf)
    %norm(tau_mapped,inf) %large magnitude?

    tau_stokes_nonpx = [tau_stokes_nonpx; tau_mapped(1:N_coarse)];
    tau_stokes_nonpy = [tau_stokes_nonpy; tau_mapped(N_coarse+1:end)];

    start_ind = size(tau_stokes_x,1);   
    
    tau_i_x = tau_mapped(1:N_coarse);
    tau_i_y = tau_mapped(N_coarse+1:end);
    tau_mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y]; 
    

    %debug mode to comment out here
    tau_stokes_x = [tau_stokes_x; tau_mapped(1:N_coarse)];
    tau_stokes_y = [tau_stokes_y; tau_mapped(N_coarse+1:end)];

    %store indices to later subtract self-interaction blocks in the matvec
    coarse_ind{i} = start_ind+1:start_ind+N_coarse; 
    rvec_in = [rvec_in; q(i)+rbase_in_c]; %Keep track of coarse source points
    
    if ~isempty(neigh)
        
        for k = 1:size(neigh,1)



            if neigh(k,2) == 2 %has already been counted in 
                break
            end
            %send in two particles with the fine grid 
            p2 = pairs(neigh(k),2); %get neighbour

            rimage_k{i} = [rimage_k{i}; rimage_vec{i,p2}];
            rimage_k{p2} = [rimage_k{p2}; rimage_vec{p2,i}];
            nimage_k{i} = [nimage_k{i}; nimage{i,p2}];
            nimage_k{p2} = [nimage_k{p2}; nimage{p2,i}];

            tau_particle_x2 = tau((p2-1)*N_large+1:N_large*p2);
            tau_particle_y2 = tau(PM+(p2-1)*N_large+1:PM+N_large*p2);

            %Keep track of local ordering of source vector for the pair
            im_nr = length(rimage_vec{i,p2});
            s_ind1_x = 1:opt.N_f; 
            s_ind2_x = opt.N_f+1:2*opt.N_f;
            s_ind1_y = 2*opt.N_f+1:3*opt.N_f;
            s_ind2_y = 3*opt.N_f+1:4*opt.N_f;

            %The ones below are different for every pair (depend on the
            %number of image points)

            t_ind_x = 4*opt.N_f+1:4*opt.N_f+2*im_nr;
            t_ind_y = 4*opt.N_f+2*im_nr+1:4*opt.N_f+4*im_nr;

            p_ind_x = 4*opt.N_f+4*im_nr+1:4*opt.N_f+6*im_nr;
            p_ind_y = 4*opt.N_f+6*im_nr+1:4*opt.N_f+8*im_nr;
 
            if two_parts 
                %Here, we build \chi 1,2 and \chi 2,1 separately. NOT
                %needed to do it like this...

                %Already computed above: 

                % %Here we use that the onebody basis is same for everybody. Otherwise, replace U{1} -> U{i}:
%                 step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
%                 tau_mapped = Y{1}*step1; 
%                 
%                 %We construct data corresponding to 1-body basis on first
%                 %particle. 

%               tau_i_x = tau_mapped(1:N_coarse);
%               tau_i_y = tau_mapped(N_coarse+1:end);
%               tau_mapped = [tau_i_x; tau_i_y]-L*[tau_i_x; tau_i_y]; 

    
    
                %Read off coarse grid contribution on other particle fine
                %grid
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,i,p2);            
                Nother = stokSLPmat(rbase_in_c+q(i),rout_fine_other,mu);
                
    
    
                R2 = -Nother*tau_mapped; %rhs on particle 2
                block = R2(1:end/2);
                A2 = [zeros(size(block)); block; zeros(size(block)); R2(end/2+1:end)]; 
                %store as x x y y 
    
                %Evaluate pseudo-inverse corresponding to fine grid
                %backward-stably.
                pair_mapped = Upf{i,p2}*A2; 
                tau_pair = Ypf{i,p2}*pair_mapped; 
    
                %% Do the same thing for the other order of the particles in the pair
                %Should be U{p2} and Y{p2} here, but now we assume every
                %particle has the same basic (coarse) discretisation. 
                step1 = U{1}*[tau_particle_x2; tau_particle_y2]; 
                mapped = Y{1}*step1;
    
                
                tau_i_x = mapped(1:N_coarse);
                tau_i_y = mapped(N_coarse+1:end);
                mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y]; 
                
    
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,p2,i);            
                Nother = stokSLPmat(rbase_in_c+q(p2),rout_fine_other,mu);
               
                R1 = -Nother*mapped;
                block = R1(1:end/2); %x-contribution           
                A1 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
                pair_mapped = Upf{i,p2}*A1;
                tau_pair2 = Ypf{i,p2}*pair_mapped;
    
                %% Project
    
                tau_mapped_tot_xi = tau_pair(s_ind1_x) + tau_pair2(s_ind1_x);
                tau_mapped_tot_yi = tau_pair(s_ind1_y) + tau_pair2(s_ind1_y);
    
                tau_mapped_tot_xp2 = tau_pair(s_ind2_x) + tau_pair2(s_ind2_x);
                tau_mapped_tot_yp2 = tau_pair(s_ind2_y) + tau_pair2(s_ind2_y);

                %Store
                tau_stress_x{i,p2} = tau_pair(t_ind_x) + tau_pair2(t_ind_x);
                tau_stress_y{i,p2} = tau_pair(t_ind_y) + tau_pair2(t_ind_y);

                tau_pot_x{i,p2} = tau_pair(p_ind_x) + tau_pair2(p_ind_x);
                tau_pot_y{i,p2} = tau_pair(p_ind_y) + tau_pair2(p_ind_y);
    
                
            else
                %Build \chi 1,2 and \chi 2,1 jointly
                %step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
                %tau_mapped = Y{1}*step1;    
                %Already determined!
 

                %Same thing with different orders for the particles in the
                %pair
                step1 = U{1}*[tau_particle_x2; tau_particle_y2]; 
                mapped = Y{1}*step1;
                %project
                tau_i_x = mapped(1:N_coarse);
                tau_i_y = mapped(N_coarse+1:end);
                mapped = [tau_i_x; tau_i_y]-Lc*[tau_i_x; tau_i_y];


                if ~precomp
                    %Read off coarse grid contribution on other particle fine
                    %grid
                    rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,i,p2);            
%                   Nother2 = stokSLPmat(rbase_in_c+q(i),rout_fine_other2,mu);
%                   R2 = -Nother2*tau_mapped; 
                
                    [udirect,vdirect] = StokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
                        real(rout_fine_other2),imag(rout_fine_other2),tau_mapped(1:end/2),...
                        tau_mapped(end/2+1:end),opt.N_c);
                    R2 = -[udirect; vdirect];
        
                    rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,opt.rad,refine,q,p2,i);            
    %                 Nother1 = stokSLPmat(rbase_in_c+q(p2),rout_fine_other1,mu);               
    %                 R1 = -Nother1*mapped;
    
                    [udirect,vdirect] = StokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
                        real(rout_fine_other1),imag(rout_fine_other1),mapped(1:end/2),...
                        mapped(end/2+1:end),opt.N_c);
                    R1 = -[udirect; vdirect];
    
                    coarse_to_fine_tot = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)]; 
                else
                    coarse_to_fine_tot = [tau_mapped(1:end/2); mapped(1:end/2); 
                        tau_mapped(end/2+1:end); mapped(end/2+1:end)];
                end
                %store as x x y y 
    
                %Take pseudoinverse of the fine representation to determine
                %beta for BOTH \chi 1,2 and \chi 2,1.
                pair_mapped = Upf{i,p2}*coarse_to_fine_tot; 
                beta_tot = Ypf{i,p2}*pair_mapped; 


                %% Project
    
                tau_mapped_tot_xi = beta_tot(s_ind1_x);
                tau_mapped_tot_yi = beta_tot(s_ind1_y);
    
                tau_mapped_tot_xp2 = beta_tot(s_ind2_x);
                tau_mapped_tot_yp2 = beta_tot(s_ind2_y);

                %Store
                tau_stress_x{i,p2} = beta_tot(t_ind_x);
                tau_stress_y{i,p2} = beta_tot(t_ind_y);

                tau_pot_x{i,p2} = beta_tot(p_ind_x);
                tau_pot_y{i,p2} = beta_tot(p_ind_y);
    
            end

            tau_fine_i = [tau_mapped_tot_xi; tau_mapped_tot_yi];
            tau_fine_p2 = [tau_mapped_tot_xp2; tau_mapped_tot_yp2];

            if use_matrix_free_projection
                proj_ok_i = true;
                proj_ok_p2 = true;
                try
                    tau_mapped_proj_i = projectOutRigid2D(tau_fine_i, rbase_in_f+q(i), q(i));
                    tau_mapped_proj_p2 = projectOutRigid2D(tau_fine_p2, rbase_in_f+q(p2), q(p2));
                    proj_ok_i = isvector(tau_mapped_proj_i) && ...
                        numel(tau_mapped_proj_i) == numel(tau_fine_i) && ...
                        all(isfinite(tau_mapped_proj_i));
                    proj_ok_p2 = isvector(tau_mapped_proj_p2) && ...
                        numel(tau_mapped_proj_p2) == numel(tau_fine_p2) && ...
                        all(isfinite(tau_mapped_proj_p2));
                catch
                    proj_ok_i = false;
                    proj_ok_p2 = false;
                end
            else
                proj_ok_i = false;
                proj_ok_p2 = false;
            end

            if ~(proj_ok_i && proj_ok_p2)
                tau_mapped_proj_i = tau_fine_i-Lf*tau_fine_i;
                tau_mapped_proj_p2 = tau_fine_p2-Lf*tau_fine_p2;
            end

            %% Evaluate flow field on the pair itself
            % Should be computed with the fine grid on the pair
            rout_pair = [rvec_out((i-1)*N_large+1:i*N_large,:); rvec_out((p2-1)*N_large+1:p2*N_large,:)];
                   
            rimage = [rimage_vec{i,p2}; rimage_vec{p2,i}];
            nimage_pair = [nimage{i,p2}; nimage{p2,i}];
        
        

            if two_parts 
                % Build matrices explicitly -> slower! 
                N_pair = stokSLPmat([rbase_in_f+q(i); rbase_in_f+q(p2)],rout_pair,mu);
                if size(rimage,1)
                    N_image = getImageKernels2D(rimage,nimage_pair,rout_pair,mu,s);
                else
                    N_image = [];
                end            

                u_pair = [N_pair N_image]*[tau_mapped+ tau_mapped2];
            else
                %Avoid constructing matrices explicitly. 
                rim = length(rimage);
                
                %get contribution from image points 
                

                u_stress = getStresslets(beta_tot(4*opt.N_f+1:4*opt.N_f+rim),...
                    beta_tot(4*opt.N_f+1+rim:4*opt.N_f+2*rim),rimage,...
                    rout_pair,real(nimage_pair),imag(nimage_pair));

                u_pot = getPotdip(beta_tot(4*opt.N_f+1+2*rim:4*opt.N_f+3*rim),...
                    beta_tot(4*opt.N_f+1+3*rim:4*opt.N_f+4*rim),rimage,rout_pair);

                %... and from fine grid of Stokeslets 
                rin_pair = [rbase_in_f+q(i); rbase_in_f+q(p2)];

%                 [u1,v1] = StokesletDirect(real(rin_pair),imag(rin_pair),...
%                     real(rout_pair),imag(rout_pair),...
%                     beta_tot(1:2*opt.N_f),beta_tot(2*opt.N_f+1:4*opt.N_f),2*opt.N_f);
                [u1,v1] = stokSLPdirect(real(rin_pair),imag(rin_pair),...
                    real(rout_pair),imag(rout_pair),...
                    [tau_mapped_proj_i(1:end/2); tau_mapped_proj_p2(1:end/2)],...
                    [tau_mapped_proj_i(end/2+1:end); tau_mapped_proj_p2(end/2+1:end)],2*opt.N_f);
                u_stok = [u1; v1];
                u_pair = u_stress+u_pot+u_stok;

                
            end
            ind1x = (i-1)*N_large+1:i*N_large;
            ind2x =  (p2-1)*N_large+1:p2*N_large;
            ind1y = (i-1)*N_large+PM+1:i*N_large+PM;
            ind2y = (p2-1)*N_large+PM+1:p2*N_large+PM;
            pair_ind = [ind1x ind2x ind1y ind2y]'; 

            u_corr(pair_ind) = u_corr(pair_ind)+u_pair; 
           

            %also subtract self-interaction on other.
            N2 = stokSLPmat(rbase_in_c+q(i),rvec_out((p2-1)*N_large+1:p2*N_large,:),1);
            u2 = N2*tau_mapped;

            N1 = stokSLPmat(rbase_in_c+q(p2),rvec_out((i-1)*N_large+1:i*N_large,:),1);
            u1 = N1*mapped; 

            u_corr(ind1x) = u_corr(ind1x)+u1(1:end/2);
            u_corr(ind2x) = u_corr(ind2x)+u2(1:end/2);
            u_corr(ind1y) = u_corr(ind1y)+u1(end/2+1:end);
            u_corr(ind2y) = u_corr(ind2y)+u2(end/2+1:end);

            %% Store source strengths at source points for the different
             % source types corresponding to the fine grid.          
            if isempty(tau_stokes_fine_x{i})
  
                tau_stokes_fine_nonpx{i} = tau_mapped_tot_xi; %to be replaced by the projected version
                tau_stokes_fine_nonpy{i} = tau_mapped_tot_yi; 

                tau_stokes_fine_x{i} = tau_mapped_proj_i(1:end/2);
                tau_stokes_fine_y{i} = tau_mapped_proj_i(end/2+1:end);
            else
                tau_stokes_fine_nonpx{i} = tau_stokes_fine_nonpx{i} + tau_mapped_tot_xi;
                tau_stokes_fine_nonpy{i} = tau_stokes_fine_nonpy{i} + tau_mapped_tot_yi;

                tau_stokes_fine_x{i} = tau_stokes_fine_x{i}+tau_mapped_proj_i(1:end/2);
                tau_stokes_fine_y{i} = tau_stokes_fine_y{i}+tau_mapped_proj_i(end/2+1:end);
            end

            if isempty(tau_stokes_fine_x{p2})
%                 tau_stokes_fine_x{i} = tau_mapped(1:end/4);
%                 tau_stokes_fine_y{i} = tau_mapped(end/2+1:3*end/4);   
                tau_stokes_fine_nonpx{p2} = tau_mapped_tot_xp2;
                tau_stokes_fine_nonpy{p2} = tau_mapped_tot_yp2; 

                tau_stokes_fine_x{p2} = tau_mapped_proj_p2(1:end/2);
                tau_stokes_fine_y{p2} = tau_mapped_proj_p2(end/2+1:end);
            else
                tau_stokes_fine_nonpx{p2} = tau_stokes_fine_nonpx{p2} + tau_mapped_tot_xp2;
                tau_stokes_fine_nonpy{p2} = tau_stokes_fine_nonpy{p2} + tau_mapped_tot_yp2;

                tau_stokes_fine_x{p2} = tau_stokes_fine_x{p2}+tau_mapped_proj_p2(1:end/2);
                tau_stokes_fine_y{p2} = tau_stokes_fine_y{p2}+tau_mapped_proj_p2(end/2+1:end);
            end 



    
        end
        
    end
       
end

has_neigh = sort(unique(pairs(:)));

tau_stress_all_x = [];
tau_stress_all_y = [];

tau_pot_all_x = [];
tau_pot_all_y = [];

tau_stress_all_px = [];
tau_stress_all_py = [];

tau_pot_all_px = [];
tau_pot_all_py = [];

%Collect all source points and source data

for i = 1:length(has_neigh)
    k = has_neigh(i);    
    
    tau_stokes_x = [tau_stokes_x; tau_stokes_fine_x{k}];
    tau_stokes_y = [tau_stokes_y; tau_stokes_fine_y{k}];

    tau_stokes_nonpx = [tau_stokes_nonpx; tau_stokes_fine_nonpx{k}];
    tau_stokes_nonpy = [tau_stokes_nonpy; tau_stokes_fine_nonpy{k}];
    
    %Keep track of the source points corresponding to the fine source data
    rvec_in = [rvec_in; q(k)+rbase_in_f];

    %How to find all neigbours of this guy
    [I,J] = find(pairs==k);
    neigh = [I J];

    for l = 1:size(neigh,1)

        if neigh(l,2) == 2 %has already been counted in 
            break
        end
        p2 = pairs(neigh(l),2);

        tau_stress_all_x = [tau_stress_all_x; tau_stress_x{k,p2}];
        tau_stress_all_y = [tau_stress_all_y; tau_stress_y{k,p2}];
    
        tau_pot_all_x = [tau_pot_all_x; tau_pot_x{k,p2}];
        tau_pot_all_y = [tau_pot_all_y; tau_pot_y{k,p2}];

        if opt.proj_all
            pot_kx = tau_pot_x{k,p2};
            pot_ky = tau_pot_y{k,p2};
            pot_1 = [pot_kx(1:end/2); pot_ky(1:end/2)];
            pot_2 = [pot_kx(end/2+1:end); pot_ky(end/2+1:end)];

            stress_kx = tau_stress_x{k,p2};
            stress_ky = tau_stress_y{k,p2};
            stress_1 = [stress_kx(1:end/2); stress_ky(1:end/2)];
            stress_2 = [stress_kx(end/2+1:end); stress_ky(end/2+1:end)];

            rim1 = rimage_vec{k,p2};
            if use_matrix_free_projection
                ok_pot1 = true;
                ok_stress1 = true;
                try
                    pot_1_new = projectOutRigid2D(pot_1, rim1, q(k));
                    stress_1_new = projectOutRigid2D(stress_1, rim1, q(k));
                    ok_pot1 = isvector(pot_1_new) && ...
                        numel(pot_1_new) == numel(pot_1) && ...
                        all(isfinite(pot_1_new));
                    ok_stress1 = isvector(stress_1_new) && ...
                        numel(stress_1_new) == numel(stress_1) && ...
                        all(isfinite(stress_1_new));
                    if ok_pot1
                        pot_1 = pot_1_new;
                    end
                    if ok_stress1
                        stress_1 = stress_1_new;
                    end
                catch
                    ok_pot1 = false;
                    ok_stress1 = false;
                end
            else
                ok_pot1 = false;
                ok_stress1 = false;
            end

            if ~(ok_pot1 && ok_stress1)
                Kim = getKmat2D(rim1,q(k)); 
                Lim = Kim*((Kim'*Kim)\Kim');
                pot_1 = pot_1-Lim*pot_1;
                stress_1 = stress_1-Lim*stress_1;
            end

            rim2 = rimage_vec{p2,k};
            if use_matrix_free_projection
                ok_pot2 = true;
                ok_stress2 = true;
                try
                    pot_2_new = projectOutRigid2D(pot_2, rim2, q(p2));
                    stress_2_new = projectOutRigid2D(stress_2, rim2, q(p2));
                    ok_pot2 = isvector(pot_2_new) && ...
                        numel(pot_2_new) == numel(pot_2) && ...
                        all(isfinite(pot_2_new));
                    ok_stress2 = isvector(stress_2_new) && ...
                        numel(stress_2_new) == numel(stress_2) && ...
                        all(isfinite(stress_2_new));
                    if ok_pot2
                        pot_2 = pot_2_new;
                    end
                    if ok_stress2
                        stress_2 = stress_2_new;
                    end
                catch
                    ok_pot2 = false;
                    ok_stress2 = false;
                end
            else
                ok_pot2 = false;
                ok_stress2 = false;
            end

            if ~(ok_pot2 && ok_stress2)
                Kim = getKmat2D(rim2,q(p2)); 
                Lim = Kim*((Kim'*Kim)\Kim');
                pot_2 = pot_2-Lim*pot_2;
                stress_2 = stress_2-Lim*stress_2;
            end

            tau_stress_all_px = [tau_stress_all_px; stress_1(1:end/2); stress_2(1:end/2)];
            tau_stress_all_py = [tau_stress_all_py; stress_1(end/2+1:end); stress_2(end/2+1:end)];

            tau_pot_all_px = [tau_pot_all_px; pot_1(1:end/2); pot_2(1:end/2)];
            tau_pot_all_py = [tau_pot_all_py; pot_1(end/2+1:end); pot_2(end/2+1:end)];



        else
            tau_stress_all_px = tau_stress_all_x;
            tau_stress_all_py = tau_stress_all_y; 

            tau_pot_all_px = tau_pot_all_x;
            tau_pot_all_py = tau_pot_all_y; 
        end
        
        %Add all image points from particle k and its neighbour
        rimage_in = [rimage_in;  rimage_vec{k,p2}; rimage_vec{p2,k}];
        nimage_in = [nimage_in;  nimage{k,p2}; nimage{p2,k}];
    end
    
end

end
