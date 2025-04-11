function [tau_stokes_x, tau_self_x, tau_beta_x,tau_cf_x,tau_stokes_y,tau_self_y,tau_beta_y,tau_cf_y,u_corr] = transform_peanut(tau,rbase_in_c,rbase_in_f,rbase_out_f,refine,rimage_vec,nimage,opt,rvec_out,rcheck_out,q,U,Y,pairs,Upf,Ypf,DC_all, YC_all,Cmap,exterior,debug)
%Transform coarse \mu -> coarse \lambda using peanut compression, to be
%used in resistance matvec

P = length(q);
N_c = opt.N_c;
N_f = opt.N_f;
N_large = length(rvec_out)/P; 
PM = length(rvec_out);
N_check = length(rcheck_out)/P;
PM2 = length(rcheck_out);


%s = opt.s; %[0 0 1 1]; %stresslet and potential dipoles at image points.
%Currently hardcoded.

%Preallocate
tau_stokes_x = zeros(N_c*P,1); 
tau_stokes_y = zeros(N_c*P,1); 

tau_beta_x = zeros(N_f*P,1);
tau_beta_y = zeros(N_f*P,1);
tau_cf_x = zeros(N_c*P,1);
tau_cf_y = zeros(N_c*P,1);


rimage_k = cell(P,1);
nimage_k = cell(P,1);

u_corr = zeros(2*N_check*P,1); 
%Store local fine grid correction (peanut compression not used locally on the pair itself)

two_parts = 0; 
precomp = opt.precomp; 
 

for i = 1:P
    coarse_ind = (i-1)*N_c+1:i*N_c; 

    %Retrieve self evaluation blocks    
    tau_particle_x = tau((i-1)*N_large+1:N_large*i);
    tau_particle_y = tau(PM+(i-1)*N_large+1:PM+N_large*i);

    %U{i}, Y{i}
    step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other?
    tau_mapped = Y{1}*step1; %this is the mapped density for this particle to throw in to the kernel
%         tau_stokes_x((i-1)*N_small+1:i*N_small) = tau_mapped(1:N_small);
%         tau_stokes_y((i-1)*N_small+1:i*N_small) = tau_mapped(N_small+1:2*N_small);
    
    % %check residual for the self-interaction only
    % NN = singleLayer(rbase_in_c+q(i),rvec_out((i-1)*N_large+1:i*N_large,:),mu);
    % disp('Self-interaction error')
    % norm(NN*tau_mapped-[tau_particle_x;tau_particle_y],inf)
    % norm((NN*tau_mapped-[tau_particle_x;tau_particle_y])./[tau_particle_x;tau_particle_y],inf)
    % norm(tau_mapped,inf)  

    tau_stokes_x(coarse_ind) = tau_mapped(1:N_c);
    tau_stokes_y(coarse_ind) = tau_mapped(N_c+1:end);
    

end


%For subtracting off self-contribution in matvec
tau_self_x = tau_stokes_x;
tau_self_y = tau_stokes_y; 

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
 
            %% First, determine beta, the fine sources
            %U{i},Y{i}
            step1 = U{1}*[tau_particle_x;tau_particle_y]; %here I assume x and y follow each other
            tau_mapped = Y{1}*step1;

           % Do a similar thing for the other order of the particles in the pair
            %U{p2}, Y{p2}
            step1 = U{1}*[tau_particle_x2;tau_particle_y2]; 
            mapped = Y{1}*step1;



            if precomp || two_parts
                %Read off coarse contribution in fine grid of other
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,i,p2); 
                %Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,mu);
                %R2 = -Nother*tau_mapped; %read off on particle 2
    
                [u2,v2] = StokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    tau_mapped(1:N_c),tau_mapped(N_c+1:2*N_c),N_c);
                R2  = -[u2; v2];
    
                rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,p2,i); 
                %Nother2 = singleLayer(rbase_in_c+q(p2),rout_fine_other,mu);
                %R1 = -Nother2*mapped; %read off on particle 1
    
                [u1,v1] = StokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
                    real(rout_fine_other),imag(rout_fine_other),...
                    mapped(1:N_c),mapped(N_c+1:2*N_c),N_c);
                R1 = -[u1; v1];
            end

            if two_parts
                %determine contributions from \chi_1 and \chi_2 separately
                block = R2(1:end/2);   
                rhs1 = [zeros(size(block)); R2(1:end/2); zeros(size(block)); R2(end/2+1:end)]; 
         
                pair_mapped = Upf{i,p2}*rhs1;
                tau_mapped = Ypf{i,p2}*pair_mapped; 
            
                %now we have beta. Lets do peanut compression from
                %here

                tau_peanut_temp = DC_all{i,p2}*tau_mapped;
                tau_peanut1 = YC_all{i,p2}*tau_peanut_temp;

                rhs2 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))]; 
                pair_mapped = Upf{i,p2}*rhs2;
                tau_mapped2 = Ypf{i,p2}*pair_mapped; 

                tau_peanut_temp = DC_all{i,p2}*tau_mapped2;
                tau_peanut2 = YC_all{i,p2}*tau_peanut_temp;
            else
                %... or jointly (faster)
                if precomp
                    rhs = [tau_mapped(1:end/2); mapped(1:end/2); 
                        tau_mapped(end/2+1:end); mapped(end/2+1:end)];
                else
                    rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
                end
       
                pair_mapped = Upf{i,p2}*rhs; 
                tau_mapped_tot = Ypf{i,p2}*pair_mapped;
    
                %now we have beta. Lets do peanut compression from
                %here

                if opt.cmap
                    tau_peanut_tot = Cmap{i,p2}*rhs;
                    
                else 
                    tau_peanut_temp = DC_all{i,p2}*tau_mapped_tot;
                    tau_peanut_tot = YC_all{i,p2}*tau_peanut_temp;
                end

                



                %store to later be used to determine force and torque
                %Need to store to close system (impose bc)
                tau_beta_x((i-1)*N_f+1:N_f*i) = tau_beta_x((i-1)*N_f+1:N_f*i)+tau_mapped_tot(1:N_f);
                tau_beta_x((p2-1)*N_f+1:N_f*p2) = tau_beta_x((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(N_f+1:2*N_f);
                tau_beta_y((i-1)*N_f+1:N_f*i) = tau_beta_y((i-1)*N_f+1:N_f*i)+tau_mapped_tot(2*N_f+1:3*N_f);
                tau_beta_y((p2-1)*N_f+1:N_f*p2) = tau_beta_y((p2-1)*N_f+1:N_f*p2)+tau_mapped_tot(3*N_f+1:4*N_f);

                %Note! Only to be done for the postprocessing
%                 tau_peanut_stokes = Cmap_F{i,p2}*rhs;
%                 tau_cf_x((i-1)*N_c+1:N_c*i) = tau_cf_x((i-1)*N_c+1:N_c*i)+tau_peanut_stokes(1:N_c);
%                 tau_cf_x((p2-1)*N_c+1:N_c*p2) = tau_cf_x((p2-1)*N_c+1:N_c*p2)+tau_peanut_stokes(N_c+1:2*N_c);
%                 tau_cf_y((i-1)*N_c+1:N_c*i) = tau_cf_y((i-1)*N_c+1:N_c*i)+tau_peanut_stokes(2*N_c+1:3*N_c);
%                 tau_cf_y((p2-1)*N_c+1:N_c*p2) = tau_cf_y((p2-1)*N_c+1:N_c*p2)+tau_peanut_stokes(3*N_c+1:4*N_c);
              
            end


            debug = 0; 
            %% Check residual in the peanut compression (if two_parts)
            if debug

                %For debugging, evaluate tau_peanut.
                %Need the peanut colloc points
                rout_peanut = create_peanut(q(i),q(p2),1060,0);
                N_peanut = singleLayer([rbase_in_c+q(i); rbase_in_c+q(p2)],rout_peanut,mu);
                u_peanut = N_peanut*tau_peanut1;
    
                %evaluate rhs
                rimage = [rimage_vec{i,p2}; rimage_vec{p2,i}];
                nimage_pair = [nimage{i,p2}; nimage{p2,i}];
                N_image = getImageKernels2D(rimage,nimage_pair,rout_peanut,mu,s);
                N_peanut_f = singleLayer([rbase_in_f+q(i); rbase_in_f+q(p2)],rout_peanut,mu);
                u_rhs = [N_peanut_f N_image]*tau_mapped;
                u_peanut2 = N_peanut*tau_peanut2;
                u_rhs2 = [N_peanut_f N_image]*tau_mapped2;


                % figure(50);
                % clf;
                % subplot(1,3,1)                
                % semilogy(abs(u_rhs+u_rhs2-u_peanut+u_peanut2)./abs(u_rhs+u_rhs2))
                % hold on
                % %semilogy(abs(u_rhs2-u_peanut2)./abs(u_rhs2+u_rhs))
                % 
                % subplot(1,3,2)                
                % semilogy(abs(u_rhs-u_peanut))
                % hold on
                % semilogy(abs(u_rhs2-u_peanut2))
                % 
                % subplot(1,3,3)
                % semilogy(abs(u_rhs))
                % hold on
                % semilogy(abs(u_rhs2))

                norm(u_rhs2-u_peanut2)
                norm(u_rhs-u_peanut)

                %magnitude of u_peanut.
                % figure(51)
                % clf;
                % subplot(1,2,1)
                % magn = abs(u_peanut(1:end/2)+1i*u_peanut(end/2+1:end));
                % scatter3(real(rout_peanut),imag(rout_peanut),magn,10,magn);
                % view(0,90)
                % subplot(1,2,2)
                % magn = abs(u_peanut2(1:end/2)+1i*u_peanut2(end/2+1:end));
                % scatter3(real(rout_peanut),imag(rout_peanut),magn,10,magn);
                % view(0,90)

            end

            %% Evaluate flow field on pair and subtract this contribution. 
            % Should be computed with the fine grid on the pair
            rout_pair = [rcheck_out((i-1)*N_check+1:i*N_check,:); rcheck_out((p2-1)*N_check+1:p2*N_check,:)];
                       
            rimage = [rimage_vec{i,p2}; rimage_vec{p2,i}];
            nimage_pair = [nimage{i,p2}; nimage{p2,i}];
            
            

            if two_parts 
                % Build matrices explicitly -> slower! 
                N_pair_c = singleLayer([rbase_in_c+q(i); rbase_in_c+q(p2)],rout_pair,mu);
                N_pair = singleLayer([rbase_in_f+q(i); rbase_in_f+q(p2)],rout_pair,mu);
                if size(rimage,1)
                    N_image = getImageKernels2D(rimage,nimage_pair,rout_pair,mu,s);
                else
                    N_image = [];
                end            
                u_peanut_corr = N_pair_c*(tau_peanut1+ tau_peanut2);
                u_pair = [N_pair N_image]*[tau_mapped+ tau_mapped2];
            else
                %Avoid constructing matrices explicitly. 
                rim = length(rimage);
                
                %get contribution from image points 
                u_stress = getStresslets(tau_mapped_tot(4*N_f+1:4*N_f+rim),...
                    tau_mapped_tot(4*N_f+1+rim:4*N_f+2*rim),rimage,...
                    rout_pair,real(nimage_pair),imag(nimage_pair));

                u_pot = getPotdip(tau_mapped_tot(4*N_f+1+2*rim:4*N_f+3*rim),...
                    tau_mapped_tot(4*N_f+1+3*rim:4*N_f+4*rim),rimage,rout_pair);

                %... and from fine grid of Stokeslets 
                rin_pair = [rbase_in_f+q(i); rbase_in_f+q(p2)];

                [u1,v1] = StokesletDirect(real(rin_pair),imag(rin_pair),...
                    real(rout_pair),imag(rout_pair),...
                    tau_mapped_tot(1:2*N_f),tau_mapped_tot(2*N_f+1:4*N_f),2*N_f);
                u_stok = [u1; v1];
                u_pair = u_stress+u_pot+u_stok;

                %Determine contribution to be subtracted

                rin_pair = [rbase_in_c+q(i); rbase_in_c+q(p2)];

                [u1,v1] = StokesletDirect(real(rin_pair),imag(rin_pair),...
                    real(rout_pair),imag(rout_pair),...
                    tau_peanut_tot(1:2*N_c),tau_peanut_tot(2*N_c+1:4*N_c),2*N_c);
                u_peanut_corr = [u1; v1];

                
            end

            pair_ind = [(i-1)*N_check+1:i*N_check (p2-1)*N_check+1:p2*N_check ...
                (i-1)*N_check+PM2+1:i*N_check+PM2 ...
                (p2-1)*N_check+PM2+1:p2*N_check+PM2]'; 
            u_corr(pair_ind) = u_corr(pair_ind)+u_pair-u_peanut_corr; 



            use_coarse = 0; 
            if use_coarse %if one uses instead the coarse data on the pair
%            %     if (sum(ismember(pairs_far,[i, p2],'rows')))
%                     %U{i}, Y{i}
                    step1 = U{1}*[tau_particle_x;tau_particle_y]; 
                    mapped = Y{1}*step1;
                    Nother = singleLayer(rbase_in_c+q(i),rout_pair(end/2+1:end,:),mu);
                    R1 = -Nother*mapped;
                    block = R1(1:end/2);
                    A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)]; 
                    
                    %For the other particle, p2
                    step1 = U{1}*[tau_particle_x2;tau_particle_y2]; 
                    mapped = Y{1}*step1;
                    Nother = singleLayer(rbase_in_c+q(p2),rout_pair(1:end/2,:),mu);
                    R1 = -Nother*mapped;
                    block = R1(1:end/2);
                    A1 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
    
    
                    u_pair2 = A1+A2;
                    norm(u_pair2-u_pair) %
            end


            s_ind1_x = 1:opt.N_c; 
            s_ind2_x = opt.N_c+1:2*opt.N_c;
            s_ind1_y = 2*opt.N_c+1:3*opt.N_c;
            s_ind2_y = 3*opt.N_c+1:4*opt.N_c;

            if two_parts

                tau_stokes_x((i-1)*N_c+1:N_c*i) = tau_stokes_x((i-1)*N_c+1:N_c*i)+...
                    tau_peanut1(s_ind1_x) + tau_peanut2(s_ind1_x);
                tau_stokes_y((i-1)*N_c+1:N_c*i) = tau_stokes_y((i-1)*N_c+1:N_c*i)+...
                    tau_peanut1(s_ind1_y) + tau_peanut2(s_ind1_y);
    
                tau_stokes_x((p2-1)*N_c+1:N_c*p2) = tau_stokes_x((p2-1)*N_c+1:N_c*p2)+...
                    tau_peanut1(s_ind2_x) + tau_peanut2(s_ind2_x);
                tau_stokes_y((p2-1)*N_c+1:N_c*p2) = tau_stokes_y((p2-1)*N_c+1:N_c*p2)+...
                    tau_peanut1(s_ind2_y) + tau_peanut2(s_ind2_y);
            else

                tau_stokes_x((i-1)*N_c+1:N_c*i) = tau_stokes_x((i-1)*N_c+1:N_c*i)+...
                    tau_peanut_tot(s_ind1_x);
                tau_stokes_y((i-1)*N_c+1:N_c*i) = tau_stokes_y((i-1)*N_c+1:N_c*i)+...
                    tau_peanut_tot(s_ind1_y);
    
                tau_stokes_x((p2-1)*N_c+1:N_c*p2) = tau_stokes_x((p2-1)*N_c+1:N_c*p2)+...
                    tau_peanut_tot(s_ind2_x);
                tau_stokes_y((p2-1)*N_c+1:N_c*p2) = tau_stokes_y((p2-1)*N_c+1:N_c*p2)+...
                    tau_peanut_tot(s_ind2_y);

            end


            
        end
        


    end
       
end

end

