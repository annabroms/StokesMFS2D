function res = matvec_mob_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rvec_out,rcheck,q,U,Y,Lc,Lf,pairs,Upf,Ypf)

[rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y, ...
    tau_stokes_nonpx, tau_stokes_nonpy,tau_stress_x,tau_stress_y,tau_stress_all_x,tau_stress_all_y,tau_stress_all_px,tau_stress_all_py,...
    tau_pot_x,tau_pot_y,tau_pot_all_x,tau_pot_all_y,tau_pot_all_px,tau_pot_all_py,u_corr]= getMobPairTransformation(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,nimage,opt,rvec_out,q,U,Y,Lc,Lf,pairs,Upf,Ypf); 
%[rvec_in,rimage_in,nimage_in,coarse_ind,tau_stokes_x,tau_stokes_y,tau_stokes_nonpx, tau_stokes_nonpy,tau_stress_x,tau_stress_y,tau_stress_all_x,tau_stress_all_y,...
 %   tau_pot_x,tau_pot_y,tau_pot_all_x,tau_pot_all_y,u_corr] 
 
P = length(q); 
PM = length(rvec_out);
N_large = PM/P;
mu = 1; 
N_c = opt.N_c;
N_f = opt.N_f;




%% Get flow field from all source types

rot = []; % no rotlets at image points so far

res = getVelocityField(rvec_in,rcheck,tau_stokes_x,tau_stokes_y,rimage_in,nimage_in,rot,...
    tau_stress_all_px,tau_stress_all_py,tau_pot_all_px,tau_pot_all_py);


two_corr = 0;

if isequal(rcheck,rvec_out)

    B = getKmat2D(rvec_out(1:N_large)-q(1),0);
    K = getKmat2D(rbase_in_f,0);
    
    %This part is already taken care of... 
    % for k= 1:P
    %     bcvec = B*K'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
    %     res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
    %     res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);
    % end
    
    %Have not taken care of this part. 
    
    if ~two_corr

        has_neigh = sort(unique(pairs(:)));
        for i = 1:length(has_neigh)
            k = has_neigh(i); 
            bcvec = B*K'*[tau_stokes_nonpx((k-1)*N_f+1+P*N_c:k*N_f+P*N_c); 
                tau_stokes_nonpy((k-1)*N_f+1+P*N_c:k*N_f+P*N_c)];
        
            res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/2);
            res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:end);


            if opt.proj_all %project also from images. To be done for all images belonging to this particle
                for j = 1:P %make sure not to add contrib twice! 
                    rim1 = rimage_vec{k,j};
                    if ~isempty(rim1)
                        rim2 = rimage_vec{j,k};
                        Kim1 = getKmat2D(rim1,q(k)); 
                        Kim2 = getKmat2D(rim2,q(j));

                        Lr = getLrPair(B,B,Kim1,Kim2);

                        stress = [tau_stress_x{k,j};tau_stress_y{k,j}]; % This will give me stresslets for both particles in the pair
                        if ~isempty(stress)
                            bcvec = Lr*stress;
    
                            pot = [tau_pot_x{k,j};tau_pot_y{k,j}];
                            bcvec = bcvec+Lr*pot;
                            res((k-1)*N_large+1:k*N_large) = res((k-1)*N_large+1:k*N_large) + bcvec(1:end/4);
                            res((j-1)*N_large+1:j*N_large) = res((j-1)*N_large+1:j*N_large) + bcvec(end/4+1:end/2);
                            res((k-1)*N_large+PM+1:k*N_large+PM) = res((k-1)*N_large+PM+1:k*N_large+PM) + bcvec(end/2+1:3*end/4);
                            res((j-1)*N_large+PM+1:j*N_large+PM) = res((j-1)*N_large+PM+1:j*N_large+PM) + bcvec(3*end/4+1:end);
                        end
                    end

                end

            end

        end


    end


    %Need to subract off the part computed twice, only from pairs not involved
    %in a pair. 
    rout = rvec_out(1:N_large)-q(1);
    rin = rbase_in_c;
    Nii = singleLayer(rin,rout,mu);
    %Nii_fine = singleLayer(rbase_in_f,rout,mu);

    for i = 1:P
        % Get sources on this particle from single layer
       % if ~isempty(coarse_ind{i})
            tau_xy = [tau_stokes_x(coarse_ind{i}); tau_stokes_y(coarse_ind{i})];
            
            %Not needed to redo for all! 
            %rout = rvec_out((i-1)*N_large+1:i*N_large,:);
            %rin = rbase_in_c+q(i);
        
            %Nii = singleLayer(rin,rout,mu);
            uii = Nii*tau_xy;
    
            %subract contribution in x
            res((i-1)*N_large+1:i*N_large) = res((i-1)*N_large+1:i*N_large)-uii(1:end/2)+tau((i-1)*N_large+1:i*N_large);
    
            %subract contribution in y
            res((i-1)*N_large+1+PM:i*N_large+PM) = res((i-1)*N_large+1+PM:i*N_large+PM)-...
                uii(end/2+1:end)+tau((i-1)*N_large+PM+1:i*N_large+PM);
     %   end
    
    %       if ~isempty(pairs)
    %           %neigh = find(pairs(:,1)==i);
    %           [I,J] = find(pairs==i);
    %           neigh = [I J];
    %       else
    %           neigh = [];
    %       end
    
        
    end

    if two_corr
        res = res-u_corr;
    end 

end


end
