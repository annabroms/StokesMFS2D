function [Lr_S,Lr_R,Lr_T,Lr_D,Sim,Rim,Tim,Dim] = getProjectedImages(q,rimage,nimage,rout,s)
%project also the image sources
warning('Description needed')

Mtot = length(rout);
pair_points(:,1) = length(rimage)/2*ones(2,1); 
pair_points(:,2) = Mtot/2*ones(2,1); 

M = 0; 
for k = 1:2

    M_k = pair_points(k,2);
    
    rout_k = rout(M+1:M+M_k);
    Kout = getKmat2D(rout_k,q(k));

    %Loop over source particles
    N = 0; 

    for i = 1:2
        ind = N+1:N+pair_points(i,1); 
        rim_i = rimage(ind);
        Kin_im = getKmat2D(rim_i,q(i));
        L_im = Kin_im*((Kin_im'*Kin_im)\Kin_im');
        
        if s(1)
            Npk_im = singleLayer(rim_i,rout_k,1);
            A_im = Npk_im-Npk_im*L_im;

            if i == k
                Lr_im = Kout*Kin_im'; 
                Lr_S(M+1:M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(1:pair_points(k,2),1:end/2);
                ind_y = sum(pair_points(:,1))+N+1:sum(pair_points(:,1))+N+pair_points(k,1);
                Lr_S(M+1:M+pair_points(k,2),ind_y) = Lr_im(1:pair_points(k,2),end/2+1:end);
                Lr_S(Mtot+M+1:Mtot+M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(pair_points(k,2)+1:end,1:end/2);
                Lr_S(Mtot+M+1:Mtot+M+pair_points(k,2),ind_y) = Lr_im(pair_points(k,2)+1:end,end/2+1:end);
            end

            Sim(M+1:M+M_k,ind) = A_im(1:end/2,1:end/2);
            Sim(M+1:M+M_k,sum(pair_points(:,1))+ind) = A_im(1:end/2,end/2+1:end);
            Sim(M+1+Mtot:M+Mtot+M_k,ind) = A_im(end/2+1:end,1:end/2);
            Sim(M+1+Mtot:M+Mtot+M_k,sum(pair_points(:,1))+ind) = A_im(end/2+1:end,end/2+1:end);

        else
            Sim = [];
            Lr_S = [];
        end

        if s(2)
           % Kin_im = [zeros(pair_points(i,1),2) ones(pair_points(i,1),1)];
            Kin_R = ones(pair_points(i,1),1);
            L_imR = Kin_R*((Kin_R'*Kin_R)\Kin_R');
            R = getImageKernels2D(rim_k,[],rout_k,1,[0 1 0 0]);
            A_im  = R-R*L_imR;
            Rim(M+1:M+M_k,ind) = A_im(1:end/2,:);
            Rim(M+1+Mtot:M+M_k+Mtot,ind) = A_im(end/2+1:end,:);

            % Need to close the system
            if i == k
                Kin_R = [zeros(pair_points(i,1),2) ones(pair_points(i,1),1)];
                Lr_Rk = Kout*Kin_R';
                Lr_R(M+1:M+M_k,ind) = Lr_Rk(1:end/2,:);
                Lr_R(M+1+Mtot:M+M_k+Mtot,ind) = Lr_Rk(end/2+1:end,:);
            end
      
        else
            Rim = [];
            Lr_R = [];
        end

        % If we want to project also for the stresslets image sources
        if s(3) 
            nim_i = nimage(ind); 
            Tpk_im = getImageKernels2D(rim_i,nim_i,rout_k,1,[0 0 1 0 0 0 0]); 
            A_im = Tpk_im-Tpk_im*L_im;

            if s(1)
                Lr_T = Lr_S;
            else
                if i == k
                    Lr_im = Kout*Kin_im'; 
                    Lr_T(M+1:M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(1:pair_points(k,2),1:end/2);
                    ind_y = sum(pair_points(:,1))+N+1:sum(pair_points(:,1))+N+pair_points(k,1);
                    Lr_T(M+1:M+pair_points(k,2),ind_y) = Lr_im(1:pair_points(k,2),end/2+1:end);
                    Lr_T(Mtot+M+1:Mtot+M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(pair_points(k,2)+1:end,1:end/2);
                    Lr_T(Mtot+M+1:Mtot+M+pair_points(k,2),ind_y) = Lr_im(pair_points(k,2)+1:end,end/2+1:end);
                end
            end

            Tim(M+1:M+M_k,ind) = A_im(1:end/2,1:end/2);
            Tim(M+1:M+M_k,sum(pair_points(:,1))+ind) = A_im(1:end/2,end/2+1:end);
            Tim(M+1+Mtot:M+Mtot+M_k,ind) = A_im(end/2+1:end,1:end/2);
            Tim(M+1+Mtot:M+Mtot+M_k,sum(pair_points(:,1))+ind) = A_im(end/2+1:end,end/2+1:end);

        else
            Tim = [];
            Lr_T = [];
        end

        %... project also for the potential dipoles?
        if s(4) 
            Dpk_im = getImageKernels2D(rim_i,nim_i,rout_k,1,[0 0 0 1 0 0 0]); 
            A_im = Dpk_im-Dpk_im*L_im;

            if s(1)
                Lr_D = Lr_S;
            elseif s(3)
                Lr_D = Lr_T;
            else
                if i == k
                    Lr_im = Kout*Kin_im'; 
                    Lr_D(M+1:M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(1:pair_points(k,2),1:end/2);
                    ind_y = sum(pair_points(:,1))+N+1:sum(pair_points(:,1))+N+pair_points(k,1);
                    Lr_D(M+1:M+pair_points(k,2),ind_y) = Lr_im(1:pair_points(k,2),end/2+1:end);
                    Lr_D(Mtot+M+1:Mtot+M+pair_points(k,2),N+1:N+pair_points(k,1)) = Lr_im(pair_points(k,2)+1:end,1:end/2);
                    Lr_D(Mtot+M+1:Mtot+M+pair_points(k,2),ind_y) = Lr_im(pair_points(k,2)+1:end,end/2+1:end);
                end
            end

            Dim(M+1:M+M_k,ind) = A_im(1:end/2,1:end/2);
            Dim(M+1:M+M_k,sum(pair_points(:,1))+ind) = A_im(1:end/2,end/2+1:end);
            Dim(M+1+Mtot:M+Mtot+M_k,ind) = A_im(end/2+1:end,1:end/2);
            Dim(M+1+Mtot:M+Mtot+M_k,sum(pair_points(:,1))+ind) = A_im(end/2+1:end,end/2+1:end);

        else
            Dim = [];
            Lr_D = [];
        end



        N = N+pair_points(i,1);


    
    end

    M = M + M_k;



end


end