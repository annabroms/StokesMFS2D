function u = eval_Stokes_solution(lambda,rin,rcheck,rimage,nimage,s,N_c,P)
%EVAL_STOKES_SOLUTION(lambda,rin,rchech,rimage,nimage)
if length(s)<6
    s(6) = 0; s(7) = 0; %no Stokes doulets, no Potential doublets
elseif length(s)<7
    s(7) = 0; 
end

    if nargin < 8
        P = 1; 
    end
    stok_x = lambda(1:P*N_c); 
    stok_y = lambda(P*N_c+1:2*P*N_c);
    im_nbr = length(rimage);
    ind_start = 2*P*N_c;
    if im_nbr
        if s(1)
            ind_stok = P*2*N_c+1:2*P*N_c+2*im_nbr;
            ind_start = 2*P*N_c+2*im_nbr;
        else
            ind_stok = [];
        end
        if s(2)
            ind_rot = ind_start+1:ind_start+im_nbr;
            ind_start = ind_start+im_nbr;
        else
            ind_rot = [];
        end
        if s(3)
            ind_stress = ind_start+1:ind_start+2*im_nbr;
            ind_start = ind_start+2*im_nbr;            
        else
            ind_stress = [];
        end
        
        if s(5)
            ind_stress = [ind_stress ind_start+1:ind_start+4*im_nbr];
            ind_start = ind_start+4*im_nbr;           
        end

        % if ~s(3) && ~s(5)
        %     ind_stress = [];
        % end

        if s(4)
            ind_pot = ind_start+1:ind_start+2*im_nbr;
            ind_start = ind_start+2*im_nbr;
        else
            ind_pot = [];
        end

        if s(7)
            ind_sd = ind_start+1:ind_start+4*im_nbr;
        else
            ind_sd = [];
        end


        
    
        %ind_stress = P*2*N_c+1:2*P*N_c+2*im_nbr;
        %ind_pot = ind_stress(end)+1:ind_stress(end)+2*im_nbr;
        stok_x = [stok_x; lambda(ind_stok(1:end/2))];
        stok_y = [stok_y; lambda(ind_stok(end/2+1:end))];





        if s(1)
            rin = [rin; rimage];
        end
        rot = lambda(ind_rot); 
        pot_x = lambda(ind_pot(1:end/2));
        pot_y = lambda(ind_pot(end/2+1:end));
                    
        if s(5)
            stress_x = [lambda(ind_stress(1:end/4)); lambda(ind_stress(end/2+1:3*end/4))];
            stress_y = [lambda(ind_stress(end/4+1:end/2)); lambda(ind_stress(3*end/4+1:end))];
            %nimage = [nimage; zeros(size(nimage))];
        else
            stress_x = lambda(ind_stress(1:end/2));
            stress_y = lambda(ind_stress(end/2+1:end));
           
        end

       % sd_x = [lambda(ind_sd(1:end/4)); lambda(ind_sd(end/2+1:3*end/4))];
        %sd_y = [lambda(ind_sd(end/4+1:end/2)); lambda(ind_sd(3*end/4+1:end))];
        sd = lambda(ind_sd);
    else
        pot_x = [];
        pot_y  = [];
        stress_x = [];
        stress_y = [];
        rot = []; 
    end

    % figure(3)
    % clf;
    % scatter(real(rin),imag(rin),20,log10(abs(stok_x+1i*stok_y)),'filled')
    % axis equal
    % colorbar
    
    %% Get and compare velocity field

    


    if s(5)
        u = getVelocityField(rin, rcheck, stok_x, stok_y, rimage, nimage, ...
                               rot,[], [],pot_x, pot_y);
        u_stress = getStresslets(stress_x,stress_y,[rimage;rimage],rcheck,real(nimage),imag(nimage));
        u = u+u_stress;
    else

        u = getVelocityField(rin, rcheck, stok_x, stok_y, rimage, nimage, ...
                               rot,stress_x, stress_y,pot_x, pot_y);
    end

    if s(7) && size(rimage,1)
        u_sd = getSourceDoublet(sd,rimage,rcheck);
        u = u+u_sd;
    end

    %if constant term
    u(1:end/2) = u(1:end/2)+lambda(end-1); 
    u(end/2+1:end) = u(end/2+1:end)+lambda(end); 




end