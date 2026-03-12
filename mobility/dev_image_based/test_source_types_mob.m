clear;
close all;
%Tests out combinations of fundamental solutions at image points for a
%Stokes 2D mobility problem

delta = 0.001; % particle - particle distance

% Set parameters
% N_c = 60; 
% tol = 1e-12;
% sep = (1/N_c)*log(1/tol);
% Rp_c = max([1-sep,0.01]); %radius of proxy surface
% a_c = 1.2; % determines number of uniform collocation points, M_c =
% a_c*N_c.
opt = get2Dparams();
opt.M_image = 35; % sets extra collocation points for close to touching region
N_c = opt.N_c;

mu = 1; %viscosity
image = 1; 
opt.image = image; 
opt.pc = 0; % no pair corr, we solve densely

%Two test geometries: 2 or 3 particles
P = 3;
if P == 2
    q = [0; 2+delta];
else
    x = 1+delta/2;
    y = sqrt((2+delta)^2-(1+delta/2)^2);
    q = [0; 2+delta; x+1i*y; x+2+delta+1i*y; 4+2*delta; x-1i*y; x+2+delta-1i*y]; 
    q = q(1:3);
   % q(2) = q(2)+0.01; %perturb the position of one of the particles
end

rad = ones(P,1); %set radii

%Set truncation level in SVD. Low truncation important for high accuracy
tol = 1e-12; 
%tol = eps; 
tol = 1e-17; %don't truncate!

% Represent all combinations of fundamental solutions at image locations
smat = dec2bin(0:15) - '0';

project_all = 1; %do the projection for all source types 

debug = 0; %plot residual and coefficient magnitude
visualise = 0; %visualise singular values?

sind = [3 7 11 13]+1; 

sind = 4;

rhs_num = 10; % Do test for many rhs
for i = 1%:size(smat,1) %loop over the combinations.
    i
    s = smat(sind(i),:); 
    %s = smat(i,:); 
    opt.s = s; 

    %% Get grid
    [rout, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base] = get2DImageGrid(q, rad, opt);

    if debug
        figure()
        plot(real(rout),imag(rout),'.')
        axis equal
        hold on
        plot(real(rimage),imag(rimage),'.')
    end

    %% Prepare system matrix

    %do projection for all Stokeslet and rotlet sources. Project also other
    %source types? Benefits?
    Ntot = buildMobilitySystem(q,rin,rout,rimage,nimage,N_c,P,pair_points,s,project_all);

    % Factorize to be able to get pseudo inverse
    [YY,UU] = getPseudoFactors(Ntot,tol,visualise);
    rng(6); %generate same rhs for each combination of source types
    for j = 1:rhs_num
        F = -1+2*rand(P,2); %random force
        T = -1+2*rand(P,1); %random torque
        %Also make sure to include these forces in the ones tested
        if (j == rhs_num) && (P == 3)
            F = [0 0; 1 0; 0 1]; T = [1; 0; 0]; %for P = 3 
        elseif (j == rhs_num) && (P == 2)
            F = [1 0; 0 0]; T = [0; 0]; %for P = 2
        end
        %% create rhs vector
        [fout,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T);    
        %% Solve system
        lambda = YY*(UU'*fout);
    
        if debug
            figure()
            semilogy(abs(lambda))
        end
        %% Extract source types
        stok_x = lambda(1:P*N_c); 
        stok_y = lambda(P*N_c+1:2*P*N_c);
        im_nbr = length(rimage);
        ind_start = 2*P*N_c;

        if s(1) && im_nbr
            rin_tot = [rin; rimage];
        else
            rin_tot = [rin];
        end
    
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
                ind_stress= [];
            end
            if s(4)
                ind_pot = ind_start+1:ind_start+2*im_nbr;
            else
                ind_pot = [];
            end
    
            stok_x = [stok_x; lambda(ind_stok(1:end/2))];
            stok_y = [stok_y; lambda(ind_stok(end/2+1:end))];
    
    
            rot = lambda(ind_rot); 
            pot_x = lambda(ind_pot(1:end/2));
            pot_y = lambda(ind_pot(end/2+1:end));
            stress_x = lambda(ind_stress(1:end/2));
            stress_y = lambda(ind_stress(end/2+1:end));
        else
            pot_x = [];
            pot_y  = [];
            stress_x = [];
            stress_y = [];
            rot = []; 
        end
    
        %% Generate flow due to RBM
    
        % Get test targets
        rcheck = [];
        N = 0;
        ntest = 4000; 
        t = linspace(0,2*pi,ntest+1);
        t = t(1:end-1)'; 
        Ktest = getKmat2D(cos(t)+1i*(sin(t)),0);
    
        u_ref = zeros(ntest*P*2,1); 

        % Projection for the proxy sources 
        Kin = getKmat2D(rin(1:N_c),0);
        L = Kin*((Kin'*Kin)\Kin');

        for k = 1:P
    
            %Loop over particles to project out the sources
            ind = (k-1)*N_c+1:k*N_c;
            lambda_k = [stok_x(ind); stok_y(ind)];
            lambda_k_tot = (eye(size(L))-L)*lambda_k;%+[lambda_x(ind);lambda_y(ind)];
            stok_x_tot(ind) = lambda_k_tot(1:end/2);
            stok_y_tot(ind) = lambda_k_tot(end/2+1:end);
    
            stok_x_all(ind) = lambda_k_tot(1:end/2)+lambda_x(ind);
            stok_y_all(ind) = lambda_k_tot(end/2+1:end)+lambda_y(ind);
    
            lambda_x_all(ind) = lambda_x(ind);
            lambda_y_all(ind) = lambda_y(ind);


            Uvec = -Kin'*lambda_k;
    
            if s(1)    
            
                rim_k = rimage(N+1:N+pair_points(k,1));
                lambda_imk = [stok_x(P*N_c+N+1:P*N_c+N+pair_points(k,1));
                    stok_y(P*N_c+N+1:P*N_c+N+pair_points(k,1))];
    
                Kim = getKmat2D(rim_k,q(k));
                Uvec = Uvec -Kim'*lambda_imk;
    
                
                Lim = Kim*((Kim'*Kim)\Kim');
                ind = P*N_c+N+1:P*N_c+N+pair_points(k,1);
                stok_im = (eye(size(Lim))-Lim)*lambda_imk;
    
                stok_x_tot(ind) = stok_im(1:end/2);
                stok_y_tot(ind) = stok_im(end/2+1:end); 
    
                %Only used for debug
                stok_x_all(ind) = stok_im(1:end/2); 
                stok_y_all(ind) = stok_im(end/2+1:end);

            end

           
            if s(2)
                rot_k = rot(N+1:N+pair_points(k,1));
                Kim = ones(pair_points(k,1),1); 
                Lim =  Kim*((Kim'*Kim)\Kim');
                rot_p = rot_k-Lim*rot_k; % + rot_vec(N+1:N+pair_points(k,1)); 
                rot_tot(N+1:N+pair_points(k,1)) = rot_p;
    
                rot_proj(N+1:N+pair_points(k,1)) = rot_p;
    
       
                Uvec = Uvec-[0;0;sum(rot_k)];
            else
                rot_tot = [];
                rot_proj = [];
            end

            if s(3) && project_all
                ind = N+1:N+pair_points(k,1);
                rim_k = rimage(ind);
                Kim = getKmat2D(rim_k,q(k));
                Lim = Kim*((Kim'*Kim)\Kim');
                stress_imk = [stress_x(ind); stress_y(ind)];
                stress_proj = stress_imk-Lim*stress_imk;
               
                stress_x_tot(ind) = stress_proj(1:end/2);
                stress_y_tot(ind) = stress_proj(end/2+1:end); 
                
                Uvec = Uvec-Kim'*stress_imk;
            elseif s(3)
                %Here, the stresslets don't contribute to velocities
                ind = N+1:N+pair_points(k,1);               
                stress_x_tot(ind) = stress_x(ind);
                stress_y_tot(ind) = stress_y(ind); 
            else
                stress_x_tot = [];
                stress_y_tot = [];
            end


            if s(4) && project_all
                ind = N+1:N+pair_points(k,1);
                rim_k = rimage(N+1:N+pair_points(k,1));
                Kim = getKmat2D(rim_k,q(k));
                Lim = Kim*((Kim'*Kim)\Kim');
                pot_imk = [pot_x(ind); pot_y(ind)];
                pot_proj = pot_imk-Lim*pot_imk;
               
                pot_x_tot(ind) = pot_proj(1:end/2);
                pot_y_tot(ind) = pot_proj(end/2+1:end); 
    
                Uvec = Uvec -Kim'*pot_imk;
            elseif s(4)
                ind = N+1:N+pair_points(k,1);
                pot_x_tot(ind) = pot_x(ind);
                pot_y_tot(ind) = pot_y(ind); 
            else
                pot_x_tot = [];
                pot_y_tot = [];
            end

            %% Now, we have the velocities for the paricle and can determine the flow at the boundary
            Up = Uvec(1:2);
            Wp = Uvec(3); 
            Uall((k-1)*3+1:k*3) = [Up; Wp]; 
            u = Ktest*[Up; Wp];
            u_ref((k-1)*ntest+1:k*ntest) = u(1:end/2);
            u_ref((k-1)*ntest+P*ntest+1:k*ntest+P*ntest) = u(end/2+1:end);  
            rcheck = [rcheck; cos(t)+1i*(sin(t))+q(k) ];
    
            N = N+pair_points(k,1); 
        
        end
    

        %% Get flow from representation
        if s(1)
            lambda_x_all = [lambda_x; zeros(sum(pair_points(:,1)),1)];
            lambda_y_all = [lambda_y; zeros(sum(pair_points(:,1)),1)];
        end
    
        u = getVelocityField(rin_tot, rcheck, stok_x_all', stok_y_all', rimage, nimage, ...
                                   rot_proj',stress_x_tot', stress_y_tot',pot_x_tot, pot_y_tot);
    
    
%         u2 = getVelocityField(rin_tot, rcheck, stok_x_tot, stok_y_tot, rimage, nimage, rot_proj',stress_x, stress_y,pot_x, pot_y);
%         [comp_flow,lambda_x,lambda_y] = getRecompletionFlow(rin,rcheck,q,F,T);
%         u2 = u2-comp_flow;  %Correct! (same as u)
    


        if debug
            figure()
            plot(u)
            hold on
            %plot(u2)
           % hold on
            plot(u_ref)
            title('Velocities')
            
            figure()
            plot(u-u_ref);
            title('Velocity error')
        end
    
    
        %% Get error
        diff_vec = u-u_ref;
        max_ref = max(abs(u_ref(1:end/2)+1i*u_ref(end/2+1:end)));
        diff_vec = (diff_vec(1:end/2)+1i*diff_vec(end/2+1:end))./max_ref; 
        
        if debug
            figure()
          
            scatter(real(rcheck),imag(rcheck),5,log10(abs(diff_vec)));
            colorbar
        end
        
        err(i,j) = max(abs(diff_vec));
        lambda_norm(i,j) = norm(lambda); 
    end
    coll_points(i) = length(rout)/P; 
end
%%

col = lines(16); % generate 16 different colors.

figure()
subplot(1,2,1)
%semilogy(Mvec,err)
for k = 1:16
    semilogy((k-1)*ones(1,rhs_num),err(k,:),'.','Color',col(k,:));
    hold on
end
ylabel('Relative residual','interpreter','latex')
xlabel('Source types at image points')
sgtitle('Mobility','interpreter','latex')
max(err)
grid on
ylim([1e-10,1e0])



subplot(1,2,2)
%semilogy(Mvec,lambda_norm)
for k = 1:16
    semilogy((k-1)*ones(1,rhs_num),lambda_norm(k,:),'.','Color',col(k,:));
    hold on
end
ylabel('$\|\vec \lambda\|$','interpreter','latex')
xlabel('Source types at image points','interpreter','latex')
grid on
ylim([1e1, 1e16])
max(lambda_norm)

%subplot(1,3,3)
%plot(Mvec,coll_points)
%semilogy(1:4,coll_points)

function Ntot = buildMobilitySystem(q,rin,rout,rim,nim,Nc,P,pair_points,s,project_all)
%BUILDMOBILITYMATRIX builds dense system matrix for a 2D Stokes mobility problem. 
%
% Nc is the number of coarse source points on the proxy boundary of each particle.
% pair_points contains in its second column and row j information on the number of
% collocation points for particle j. s is a vector that sets the
% singularity types along image lines.

%Image enhancement?
if size(rim)
    image = 1;
else
    image = 0;
end

%Create projection for the particle(s) for which a force is provided (the
%moving bodys)
rin_k = rin(1:Nc);
Kin = getKmat2D(rin_k,q(1));
L = Kin*((Kin'*Kin)\Kin'); %Projection matrix, same for everyone as it depends only on the proxy sources

M = 0;
Mtot = sum(pair_points(:,2));

%Loop over all target particles
Lr_all = zeros(2*sum(pair_points(:,2)),2*Nc*(P));
for k = 1:P

    M_k = pair_points(k,2);
    
    rout_k = rout(M+1:M+M_k);
    Kout = getKmat2D(rout_k,q(k));
    Lr = Kout*Kin'; 

    Lr_all(M+1:M+pair_points(k,2),Nc*(k-1)+1:Nc*k) = Lr(1:pair_points(k,2),1:end/2);
    Lr_all(M+1:M+pair_points(k,2),Nc*P+Nc*(k-1)+1:Nc*P+Nc*k) = Lr(1:pair_points(k,2),end/2+1:end);

    Lr_all(Mtot+M+1:Mtot+M+pair_points(k,2),Nc*(k-1)+1:Nc*k) = Lr(pair_points(k,2)+1:end,1:end/2);
    Lr_all(Mtot+M+1:Mtot+M+pair_points(k,2),Nc*P+Nc*(k-1)+1:Nc*P+Nc*k) = Lr(pair_points(k,2)+1:end,end/2+1:end);

    %Loop over source particles
    N = 0; 
    for i = 1:P
        
        Npk = stokSLPmat(rin((i-1)*Nc+1:i*Nc),rout_k,1);
        A = Npk*(eye(size(L))-L);

        Sproxy(M+1:M+M_k,(i-1)*Nc+1:i*Nc) = A(1:end/2,1:end/2);
        Sproxy(M+1:M+M_k,P*Nc+(i-1)*Nc+1:P*Nc+i*Nc) = A(1:end/2,end/2+1:end);
        Sproxy(M+1+Mtot:M+Mtot+M_k,(i-1)*Nc+1:i*Nc) = A(end/2+1:end,1:end/2);
        Sproxy(M+1+Mtot:M+Mtot+M_k,P*Nc+(i-1)*Nc+1:P*Nc+i*Nc) = A(end/2+1:end,end/2+1:end);
    
        if size(rim,1)
            % Project also for image sources
            if s(1)
                ind = N+1:N+pair_points(i,1);
                rim_i = rim(ind);
                Kin_im = getKmat2D(rim_i,q(i));
                L_im = Kin_im*((Kin_im'*Kin_im)\Kin_im');
                Npk_im = stokSLPmat(rim_i,rout_k,1);
                A_im = Npk_im*(eye(size(L_im))-L_im);

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
                ind = N+1:N+pair_points(i,1);
                rim_k = rim(ind);
               % Kin_im = [zeros(pair_points(i,1),2) ones(pair_points(i,1),1)];
                Kin_im = ones(pair_points(i,1),1);
                L_im = Kin_im*((Kin_im'*Kin_im)\Kin_im');
                R = getImageKernels2D(rim_k,[],rout_k,1,[0 1 0 0]);
                A_im  = R*(eye(size(L_im))-L_im);
                Rim(M+1:M+M_k,ind) = A_im(1:end/2,:);
                Rim(M+1+Mtot:M+M_k+Mtot,ind) = A_im(end/2+1:end,:);

                % Need to close the system
                if i == k
                  %  Lr_im = Kout(:,1)*Kin_im';
                    Kin_im = [zeros(pair_points(i,1),2) ones(pair_points(i,1),1)];
                    %Kout_R = [zeros(2*pair_points(i,2),2) ones(2*pair_points(i,2),1)];
                    %Kout_R = [ones(2*pair_points(i,2),2) ones(2*pair_points(i,2),1)];
                    Lr_im = Kout*Kin_im';
                    %Lr_im = ones(2*pair_points(k,2),pair_points(k,1));
                    Lr_R(M+1:M+M_k,ind) = Lr_im(1:end/2,:);
                    Lr_R(M+1+Mtot:M+M_k+Mtot,ind) = Lr_im(end/2+1:end,:);

                end

                
            else
                Rim = [];
                Lr_R = [];
            end

            % If we want to project also for the stresslets image sources
            if s(3) && project_all
                ind = N+1:N+pair_points(i,1);
                rim_i = rim(ind);
                nim_i = nim(ind); 
                Kin_im = getKmat2D(rim_i,q(i));
                L_im = Kin_im*((Kin_im'*Kin_im)\Kin_im');
                Tpk_im = getImageKernels2D(rim_i,nim_i,rout_k,1,[0 0 1 0]); 
                A_im = Tpk_im*(eye(size(L_im))-L_im);

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
            if s(4) && project_all
                ind = N+1:N+pair_points(i,1);
                rim_i = rim(ind);
                nim_i = nim(ind); 
                Kin_im = getKmat2D(rim_i,q(i));
                L_im = Kin_im*((Kin_im'*Kin_im)\Kin_im');
                Dpk_im = getImageKernels2D(rim_i,nim_i,rout_k,1,[0 0 0 1]); 
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

        else
            Sim = [];
            Rim = [];
            Tim = [];
            Dim  = [];

        end

        N = N+pair_points(i,1);


    
    end

    M = M + M_k;



end

N = Sproxy+Lr_all;


%%


if image && ~project_all  

    Tim = getImageKernels2D(rim,nim,rout,1,[0 0 s(3) 0]); 
    Dim = getImageKernels2D(rim,nim,rout,1,[0 0 0 s(4)]); 
    Lr_D = zeros(size(Dim));
    Lr_T = zeros(size(Tim));
end
Ntot = [N Sim+Lr_S Rim+Lr_R Tim+Lr_T Dim+Lr_D];

end


