function [rvec_in,coarse_ind,tau_stokes_x,tau_stokes_y,u_corr] = getPairTransformationStokes(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,opt,rvec_out,q,U,Y,pairs,Upf,Ypf)
%GETPAIRTRANSFORMATIONSTOKES maps data at coarse collocation nodes back to coarse and fine source
%strengths, in preparation for the resistance matvec. Uses stokeslet
%sources only. 

%Very long function call...


P = opt.P; 
PM = length(rvec_out);
N_large = PM/P; 
N_coarse = opt.N_c;

two_parts = 0; %Determine \chi 1,2 and \chi 2,1 separately? 
%Makes more sense computationally to determine both at once. 

%Bookkeeping stuff

tau_stokes_x = []; %should be preallocated...
tau_stokes_y = []; 

rvec_in = []; 

coarse_ind = cell(P,1); 
tau_stokes_fine_x = cell(P,1); 
tau_stokes_fine_y = cell(P,1); 

tau_stokes_extra_x = cell(P,1); 
tau_stokes_extra_y = cell(P,1); 

rimage_k = cell(P,1);

u_corr = zeros(2*PM,1); 
%Store local fine grid correction

precomp = opt.precomp; 

 
for i = 1:P
    
    %Get data on coarse collocation points for this particle
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %Retrieve self evaluation blocks
    step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
    tau_mapped = Y{1}*step1; %this is the mapped density (data at source points)
    % for this particle to throw in to the Stokeslet kernel,corresponding to the coarse grid.
    
    %check residual for the self-interaction only
    %NN = singleLayer(rbase_in_c+q(i),rvec_out((i-1)*N_large+1:i*N_large,:),mu);
    %disp('Self-interaction error')
    %norm(NN*tau_mapped-[tau_particle_x;tau_particle_y],inf)
    %norm((NN*tau_mapped-[tau_particle_x;tau_particle_y])./[tau_particle_x;tau_particle_y],inf)
    %norm(tau_mapped,inf)
    start_ind = size(tau_stokes_x,1);   

    tau_stokes_x = [tau_stokes_x; tau_mapped(1:N_coarse)];
    tau_stokes_y = [tau_stokes_y; tau_mapped(N_coarse+1:end)];

    %store indices to later subtract self-interaction blocks in the matvec
    coarse_ind{i} = start_ind+1:start_ind+N_coarse; 
    rvec_in = [rvec_in; q(i)+rbase_in_c]; %Keep track of coarse source points


    %check if particle is in any pair.
    if ~isempty(pairs)
        %neigh = find(pairs(:,1)==i);
        [I,J] = find(pairs==i);
        neigh = [I J];
    else
        neigh = [];
    end
    
    if ~isempty(neigh)
        
        %loop over all close neigbours of this particle
        for k = 1:size(neigh,1)

            if neigh(k,2) == 2 %has already been counted in 
                break
            end
            
            p2 = pairs(neigh(k),2); %get neighbour index

            %image points and the random vectors used as partial data for the stresslets
            rimage_k{i} = [rimage_k{i}; rimage_vec{i,p2}];
            rimage_k{p2} = [rimage_k{p2}; rimage_vec{p2,i}];

            %coarse collocation point data on the other particle
            tau_particle_x2 = tau((p2-1)*N_large+1:N_large*p2);
            tau_particle_y2 = tau(PM+(p2-1)*N_large+1:PM+N_large*p2);

            %Keep track of local ordering of source vector for the pair
            %TODO: modify here!
            im_nr = length(rimage_vec{i,p2});
            % s_ind1_x = 1:opt.N_f+im_nr; 
            % s_ind2_x = opt.N_f+im_nr+1:2*opt.N_f+2*im_nr;
            % s_ind1_y = 2*opt.N_f+2*im_nr+1:3*opt.N_f+3*im_nr;
            % s_ind2_y = 3*opt.N_f+3*im_nr+1:4*opt.N_f+4*im_nr;

            %nodes on uniform curve
            s_ind1_x = 1:opt.N_f; 
            s_ind2_x = opt.N_f+im_nr+1:2*opt.N_f+im_nr;
            s_ind1_y = 2*opt.N_f+2*im_nr+1:3*opt.N_f+2*im_nr;
            s_ind2_y = 3*opt.N_f+3*im_nr+1:4*opt.N_f+3*im_nr;

            %enhancing nodes
            e_ind1_x = opt.N_f+1:opt.N_f+im_nr; 
            e_ind2_x = 2*opt.N_f+im_nr+1:2*opt.N_f+2*im_nr;
            e_ind1_y = 3*opt.N_f+2*im_nr+1:3*opt.N_f+3*im_nr;
            e_ind2_y = 4*opt.N_f+3*im_nr+1:4*opt.N_f+4*im_nr;


            if two_parts
                %Here, we build \chi 1,2 and \chi 2,1 separately. NOT
                %needed to do it like this...

                % Use that the one body basis is same for everybody. Otherwise, replace U{1} -> U{i}:
                step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
                tau_mapped = Y{1}*step1;     
                %Already determined?
    
                %Read off coarse grid contribution on other particle fine
                %grid
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,i,p2);            
                Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,mu);
    
                R2 = -Nother*tau_mapped; %rhs on particle 2
                block = R2(1:end/2);
                a2 = [zeros(size(block)); block; zeros(size(block)); R2(end/2+1:end)]; 
                %store as x x y y 
    
                %Evaluate pseudo-inverse corresponding to fine grid
                %backward-stably.
                pair_mapped = Upf{i,p2}*a2; 
                tau_mapped = Ypf{i,p2}*pair_mapped; 
    
                %% Do the same thing for the other order of the particles in the pair
                %Should be U{p2} and Y{p2} here, but now we assume every
                %particle has the same basic (coarse) discretisation. 
                step1 = U{1}*[tau_particle_x2; tau_particle_y2]; 
                mapped = Y{1}*step1;
    
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,p2,i);            
                Nother = singleLayer(rbase_in_c+q(p2),rout_fine_other,mu);
               
                R1 = -Nother*mapped;
                block = R1(1:end/2); %x-contribution
                a1 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];           
                pair_mapped = Upf{i,p2}*a1;
                tau_mapped2 = Ypf{i,p2}*pair_mapped;
                
                %Store source strengths at source points for the different
                %source types corresponding to the fine grid. 

                warn('warning: old!')
                if isempty(tau_stokes_fine_x{i})
      
                    tau_stokes_fine_x{i} = tau_mapped2(s_ind1_x) + tau_mapped(s_ind1_x);
                    tau_stokes_fine_y{i} = tau_mapped2(s_ind1_y) + tau_mapped(s_ind1_y); 
                else
                    tau_stokes_fine_x{i} = tau_stokes_fine_x{i} + tau_mapped(s_ind1_x)+tau_mapped2(s_ind1_x);
                    tau_stokes_fine_y{i} = tau_stokes_fine_y{i} + tau_mapped(s_ind1_y) + tau_mapped2(s_ind1_y);
                end
    
                if isempty(tau_stokes_fine_x{p2})  
                    tau_stokes_fine_x{p2} = tau_mapped2(s_ind2_x) + tau_mapped(s_ind2_x);
                    tau_stokes_fine_y{p2} = tau_mapped2(s_ind2_y) + tau_mapped(s_ind2_y); 
                else
                    tau_stokes_fine_x{p2} = tau_stokes_fine_x{p2} + tau_mapped(s_ind2_x) + tau_mapped2(s_ind2_x);
                    tau_stokes_fine_y{p2} = tau_stokes_fine_y{p2} + tau_mapped(s_ind2_y) + tau_mapped2(s_ind2_y);
                end 
    

                
            
            else
                %Build \chi 1,2 and \chi 2,1 jointly
                %step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
                %tau_mapped = Y{1}*step1;    
                %Already determined!

                %Same thing with different orders for the particles in the
                %pair
                step1 = U{1}*[tau_particle_x2; tau_particle_y2]; 
                mapped = Y{1}*step1;

                if ~precomp    
                    %Read off coarse grid contribution on other particle fine grid
                    rout_fine_other2 = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,i,p2);            
                    %Nother2 = singleLayer(rbase_in_c+q(i),rout_fine_other2,mu);   
                    %R2 = -Nother2*tau_mapped; 
    
                    [u2,v2] = StokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
                    real(rout_fine_other2),imag(rout_fine_other2),...
                    tau_mapped(1:N_coarse),tau_mapped(N_coarse+1:2*N_coarse),N_coarse);
                    R2  = -[u2; v2];
    
                    rout_fine_other1 = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,p2,i);            
    %                 Nother1 = singleLayer(rbase_in_c+q(p2),rout_fine_other1,mu);               
    %                 R1 = -Nother1*mapped;
    
                    [u1,v1] = StokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
                        real(rout_fine_other1),imag(rout_fine_other1),...
                        mapped(1:N_coarse),mapped(N_coarse+1:2*N_coarse),N_coarse);
                    R1  = -[u1; v1];
    
                    coarse_to_fine_tot = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)]; 
                    %store as x x y y 
                else
                    coarse_to_fine_tot = [tau_mapped(1:end/2); mapped(1:end/2); 
                        tau_mapped(end/2+1:end); mapped(end/2+1:end)];
                end

    
                %Take pseudoinverse of the fine representation to determine
                %beta for BOTH \chi 1,2 and \chi 2,1.
                pair_mapped = Upf{i,p2}*coarse_to_fine_tot; 
                beta_tot = Ypf{i,p2}*pair_mapped; 

  
                %% Evaluate flow field on the pair itself
                % Should be computed with the fine grid on the pair
                rout_pair = [rvec_out((i-1)*N_large+1:i*N_large,:); rvec_out((p2-1)*N_large+1:p2*N_large,:)];

                if two_parts 
                    % Build matrices explicitly -> slower! 
                    N_pair = singleLayer([rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}],rout_pair,mu);
                    u_pair = N_pair*(tau_mapped+ tau_mapped2);
                else
                    %Avoid constructing matrices explicitly. 
  
                    rin_pair = [rbase_in_f+q(i); rimage_vec{i,p2}; rbase_in_f+q(p2); rimage_vec{p2,i}];
                    st_all = length(rin_pair);    
                    [u1,v1] = stokesletDirect(real(rin_pair),imag(rin_pair),...
                        real(rout_pair),imag(rout_pair),...
                        beta_tot(1:st_all),beta_tot(st_all+1:2*st_all),st_all);
                    u_pair = [u1; v1];
    
                    
                end
                ind1x = (i-1)*N_large+1:i*N_large;
                ind2x =  (p2-1)*N_large+1:p2*N_large;
                ind1y = (i-1)*N_large+PM+1:i*N_large+PM;
                ind2y = (p2-1)*N_large+PM+1:p2*N_large+PM;
                pair_ind = [ind1x ind2x ind1y ind2y]'; 

                u_corr(pair_ind) = u_corr(pair_ind)+u_pair; 
               

                %also subtract self-interaction on other.
                N2 = singleLayer(rbase_in_c+q(i),rvec_out((p2-1)*N_large+1:p2*N_large,:),1);
                u2 = N2*tau_mapped;

                N1 = singleLayer(rbase_in_c+q(p2),rvec_out((i-1)*N_large+1:i*N_large,:),1);
                u1 = N1*mapped; 

                u_corr(ind1x) = u_corr(ind1x)+u1(1:end/2);
                u_corr(ind2x) = u_corr(ind2x)+u2(1:end/2);
                u_corr(ind1y) = u_corr(ind1y)+u1(end/2+1:end);
                u_corr(ind2y) = u_corr(ind2y)+u2(end/2+1:end);
                
               
                %tau_stokes_fine is not any more the same for every pair!

                %Store source strengths for later evaluation.
                if isempty(tau_stokes_fine_x{i})
      
                    tau_stokes_fine_x{i} = beta_tot(s_ind1_x);
                    tau_stokes_fine_y{i} = beta_tot(s_ind1_y);
                else
                    tau_stokes_fine_x{i} = tau_stokes_fine_x{i} + beta_tot(s_ind1_x);
                    tau_stokes_fine_y{i} = tau_stokes_fine_y{i} + beta_tot(s_ind1_y);
                end

    
                if isempty(tau_stokes_fine_x{p2})   
                    tau_stokes_fine_x{p2} = beta_tot(s_ind2_x);
                    tau_stokes_fine_y{p2} = beta_tot(s_ind2_y); 
                else
                    tau_stokes_fine_x{p2} = tau_stokes_fine_x{p2} + beta_tot(s_ind2_x);
                    tau_stokes_fine_y{p2} = tau_stokes_fine_y{p2} + beta_tot(s_ind2_y);
                end 

                tau_stokes_extra_x{i} = [tau_stokes_extra_x{i}; beta_tot(e_ind1_x)];
                tau_stokes_extra_y{i} = [tau_stokes_extra_y{i}; beta_tot(e_ind1_y)];
                tau_stokes_extra_x{p2} = [tau_stokes_extra_x{p2}; beta_tot(e_ind2_x)];
                tau_stokes_extra_y{p2} = [tau_stokes_extra_y{p2}; beta_tot(e_ind2_y)];
    
            end

        end
       
    end

end


has_neigh = sort(unique(pairs(:)));

%Collect all source points and source data

for i = 1:length(has_neigh)
    k = has_neigh(i);    
    
    %are tau_stokes_x and tau_stokes_fine_x at the same set of nodes or
    %different?

    tau_stokes_x = [tau_stokes_x; tau_stokes_fine_x{k}; tau_stokes_extra_x{k}];
    tau_stokes_y = [tau_stokes_y; tau_stokes_fine_y{k}; tau_stokes_extra_y{k}];

    %Keep track of the source points corresponding to the fine source data
    rvec_in = [rvec_in; q(k)+rbase_in_f; rimage_k{k}];  

    %How to find all neigbours of this guy
    [I,J] = find(pairs==k);
    neigh = [I J];

    for l = 1:size(neigh,1)

        if neigh(l,2) == 2 %has already been counted in 
            break
        end
        p2 = pairs(neigh(l),2);

        
        
        %Add all extra points from particle k. Not sure these will be
        %needed.
       % rvec_in = [rvec_in;  rimage_vec{k,p2}; rimage_vec{p2,k}];

    end
    
end

end
