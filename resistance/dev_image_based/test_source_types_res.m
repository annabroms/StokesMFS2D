clear;
close all;

% This script tests the combination of fundamental solutions to be used for
% the fine grid for a close to touching pair of particles in the Stokes
% resistance problem. 

delta = 0.001; % particle-particle distance
D = 2+delta;
S1 = D^2;
xa1 = (S1/D/2-sqrt((S1^2)/D^2/4-1))/D;
acc = xa1*D;

%get parameters
opt = get2Dparams();
opt.N_c = 150;
%opt.N_c = 60;

%set separation between true boundary and proxy surface
tol = 1e-12;
sep = (1/opt.N_c)*log(1/tol); %separation between surfaces based on rule of thumb in Stein & Barnett QFS paper (2021) 
opt.Rp_c = max([1-sep,0.01]); %radius of proxy surface

%opt.Rp_c = 0.8151;
%opt.Rp_c = 0.75; % proxy radius
%opt.Rp_c = 0.53;
N_c = opt.N_c;  % should be 60 for the test

image = 1; %use images!
normalize = 1; % Do column scaling
left_precond = 0; % Do row weighting

mu = 1; % viscosity
P = 2; % number of particles

% Prepare for analysing singularities
if P == 2
    q = [0; 2+delta];
    N = 1e4;

    %beta = log(10/sqrt(delta));
    %beta = log(1/delta);
    %N = 50*round(beta);
    %N = max([N,ceil(1.2*im_types*sample_nbr)]);
   % [t, ~] = gradedptr(N, beta);
    y = [logspace(-12,0,N)'; -logspace(-12,0,N)'];
    y = sort(y);  
   % y = t'-pi;
    rfield = 1+delta/3+1i*y;
    tau = 1e-4;
    rfield2 = rfield*exp(1i*tau);
    rfield3 = rfield*exp(-1i*tau);
   % rfield = [rfield; rfield2; rfield3];
   rfield = rfield3; 
else
    % three particles in a triangle
    x = 1+delta/2;
    y = sqrt((2+delta)^2-(1+delta/2)^2);
    q = [0; 2+delta; x+1i*y; x+2+delta+1i*y; 4+2*delta; x-1i*y; x+2+delta-1i*y]; 
    q = q(1:3);

    mid1 = (2+delta-(1+delta/2+1i*y))/2;
    mid2 = (1+delta/2+1i*y)/2;



    N = 5e3;
    y = 5*[logspace(-12,-1,N)'; -logspace(-12,-1,N)'];

  
    y = sort(y);  
   % y = t'-pi;
    rfield = 1i*y;

    rfield2 = rfield*exp((pi/6+pi/2)*1i)+mid1'+1+delta/2;

    rfield3 = rfield*exp((pi/6+pi/2+2*pi/3)*1i)+mid2;%+1+delta/2;

    %rfield = [1+delta/2+rfield;rfield2; rfield3];


   % q(2) = q(2)+0.01; %perturb the position of one of the particles
end

rads = ones(P,1);  % particle radii
%Mvec = [5:5:100];  % sets different number of collocation points in close to touching region
tol = 1e-14; %low truncation here is important! 

%Binary vector to represent combination fo fundamental solutions.
smat = dec2bin(0:15) - '0'; % generates matrix for all possible combination of fundamental solutions [S R Tr D T]
smat = [smat zeros(16,1)];
smat = dec2bin(0:31) - '0';


debug = 0; 
sind = [3 7 11 13 15]+1;
%sind = 11; 
sind = 4;
rhs_num = 1;

alpha_vec = linspace(0.01,pi/6,10);
alpha_vec = linspace(0.01,pi/3,10);
opt.n_clusters = 40; 

smat = [1 0 1 1 0 0 0];

for i = 1:size(smat,1)
%for i = length(alpha_vec) % Need to modify get2DImageGrid first to loop
%over alpha

    i
    
    %M_image = 35;
%    opt.alpha = alpha_vec(i); 
    %opt.alpha = pi/6;
    %opt.alpha = 0.01; 
    %opt.M_image = 30;
    
   % s = smat(sind(i),:); 
    s = smat(i,:); 
    if s(:,3) && s(:,5)
        err(i,j) = inf;
        lambda_norm(i,j) = inf; 

        continue;
    end
    s = smat(i,:); 
   % %s = [0 0 0 1 1];
    %s = [0 0 1 1 0];
     
    visualise = 0; 
    opt.s =s;
    
    %% Get grid
    
    % The number of sources will have an impact on the number of
    % collocation points. 
    [rout, weights, rin, rimage, nimage, pair_points, pairs, rimage_pairs, refine, rin_base] = get2DImageGrid( ...
        q, rads,opt);
    
    if debug
        figure(5)
        clf;
        plot(real(rout),imag(rout),'.')
        axis equal
        hold on
        plot(real(rin),imag(rin),'.')
        plot(real(rimage),imag(rimage),'.')
    end
    
    
    %% Prepare system matrix
    Nio = stokSLPmat(rin,rout,mu);
    if size(rimage,1)
        Nimage = getImageKernels2D(rimage,nimage,rout,mu,s);
    else
        Nimage = [];
    end
    
    Ntot = [Nio Nimage];
    
    if normalize
        col_norms = vecnorm(Ntot);
        D = diag(1 ./ col_norms);
        
        Ntot = Ntot*D;
        
    end
    
    w = repmat(sqrt(weights),2,1);
    if left_precond
        
        Ntot = diag(w)*Ntot;
    end
    
    
    
    % Get pseudo inverse
    tic
    [YY,UU] = getPseudoFactors(Ntot,tol,visualise);
    toc    
    
    rng(6); 
    p = 0.05; %determines velocity magnitudes
    p = 1; 
    %Repeat experiment for multiple rhs

    if s(5)
        nimage = [ones(size(nimage)); 1i*ones(size(nimage))];
    end

    for j = 1:rhs_num % Repeat for different rhs
    
        Up = -p+2*p*rand(P,2);
        Wp = -p+2*p*rand(P,1); 
        if (j == rhs_num) && (P == 3)
            Up = [0 0; 1 0; 0 1]; Wp = [1; 0; 0]; %for P = 3
        elseif (j == rhs_num) && (P == 2)
            Up = [0 0; 1 0]; Wp = [1; 0]; %for P = 2
        end
        %% create rhs vector
        start = 1; 
        ntest = 4000; 
        %ntest = 1000; %testing to see the difference
        %ntest = N_c*a_c; 
        t = linspace(0,2*pi,ntest+1);
        t = t(1:end-1)'; 
        u_ref = zeros(P*ntest*2,1); 
        Ktest = getKmat2D(cos(t)+1i*(sin(t)),0);
        rcheck = [];
        fout = zeros(sum(pair_points(:,2)),1); 
        for k = 1:P
           stop = sum(pair_points(1:k,2));
           rhs_f = @(x) [Up(k,1)-Wp(k)*(imag(x)-imag(q(k))); Up(k,2)+Wp(k)*(real(x-q(k)))];
           rhs_k = rhs_f(rout(start:stop));
           fout(start:stop) = rhs_k(1:end/2);
           fout(sum(pair_points(:,2))+start:sum(pair_points(:,2))+stop) = rhs_k(end/2+1:end);
        
           start = stop+1;
           u = Ktest*[Up(k,:)'; Wp(k)];
           u_ref((k-1)*ntest+1:k*ntest) = u(1:end/2);
           u_ref((k-1)*ntest+P*ntest+1:k*ntest+P*ntest) = u(end/2+1:end);  
           rcheck = [rcheck; cos(t)+1i*(sin(t))+q(k) ];
        
        end
        
        
        %% Solve system
        if left_precond
            fout = repmat(sqrt(weights),2,1).*fout;
        end
        lambda = YY*(UU'*fout);

    
        if normalize
            lambda = D*lambda;
        end
        
        %% Evaluate resulting flow
        if debug
            figure()
            semilogy(abs(lambda))
        end


        u = eval_Stokes_solution(lambda,rin,rcheck,rimage,nimage,s,N_c,P);

        if debug
    
            %Identify singularity locations using AAA
    
            u_field = getVelocityField(rin, rfield, stok_x, stok_y, rimage, nimage, ...
                                       rot,stress_x, stress_y,pot_x, pot_y);
        
            %Try to fit a rational approximation to ftest an locate the poles 
            comp_f = u_field(end/2+1:end);
            %comp_f = u_field(1:end/2);
          
          %  [r, pol, res, zer, zj, fj, wj, errvec, wt] = aaa(comp_f,rfield,'degree',14);
            [r, pol, res, zer, zj, fj, wj, errvec, wt] = aaa(comp_f,rfield,'tol',1e-6);
        
            
            figure()
            subplot(1,2,1)
            plot(real(rfield),real(r(rfield)),'+')
            hold on
            plot(real(rfield),comp_f,'o')
        
            subplot(1,2,2)
            semilogy(imag(rfield),abs(real(r(rfield))-comp_f),'+-')
        
        
            figure(5)
            %plot(real(pol),imag(pol),'k.','MarkerSize',10);
            scatter(real(pol),imag(pol),30,log10(abs(res)),'filled');
            colorbar
            plot(real(rfield),imag(rfield),'k.')
            title('Stokes: AAA residue magnitude log10')
        
        
            figure()
            subplot(1,2,1)
            scatter(real(rfield),imag(rfield),5,real(r(rfield)));
            colorbar
        
            subplot(1,2,2)
            scatter(real(rfield),imag(rfield),5,imag(r(rfield)));
            colorbar
        end
    
     
        if debug
            figure()
            plot(u)
            hold on
            plot(u_ref)
            
            figure(6)
            plot(u-u_ref);
            title('Surface residual')
        end
        
        diff_vec = u-u_ref;
        max_ref = max(abs(u_ref(1:end/2)+1i*u_ref(end/2+1:end)));
        diff_vec = (diff_vec(1:end/2)+1i*diff_vec(end/2+1:end))./max_ref; 
        
        if debug
            figure()
          
            scatter(real(rcheck),imag(rcheck),5,log10(abs(diff_vec)));
            colorbar
        end
        
        err(i,j) = max(abs(diff_vec));
        lambda_norm(i,j) = norm(lambda,inf); 
    end
    coll_points(i) = length(rout)/P; 
end
%%
%fprintf('Particles delta = %1.2e apart, surface residual %1.2e, max lambda %1.2e\n AAA error %1.2e, degree %u, pol to acc dist: %1.2e\n',...
 %   delta,err(end),max(lambda_norm),errvec(end),length(pol),abs(min(pol-acc)))


col = lines(32); 

figure()
subplot(1,2,1)
%semilogy(Mvec,err)
for k = 1:32
    semilogy((k-1)*ones(1,rhs_num),err(k,:),'.','Color',col(k,:));
    hold on
end
xlim([0,30])
ylabel('Relative residual','interpreter','latex')
xlabel('Source types at image points')
sgtitle('Resistance','interpreter','latex')
max(err)
grid on

subplot(1,2,2)
%semilogy(Mvec,lambda_norm)
for k = 1:32
    semilogy((k-1)*ones(1,rhs_num),lambda_norm(k,:),'.','Color',col(k,:));
    hold on
end
xlim([0,30])
ylabel('$\|\vec \lambda\|$','interpreter','latex')
xlabel('Source types at image points','interpreter','latex')
grid on
max(lambda_norm)

%subplot(1,3,3)
%plot(Mvec,coll_points)
%semilogy(1:4,coll_points)


%% Visualise against alpha

figure()
subplot(1,2,1)
%semilogy(Mvec,err)
semilogy(alpha_vec,err,'r.')
xticks([pi/6,pi/3])
xticklabels({'\pi/6','\pi/3'});
% for k = 1:16
%     semilogy(alpha_vec,err(k,:),'b.');
%     hold on
% end
ylabel('Relative residual','interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
sgtitle('Resistance','interpreter','latex')
max(err)
grid on

subplot(1,2,2)
semilogy(alpha_vec,lambda_norm,'r.')
%semilogy(Mvec,lambda_norm)
% for k = 1:16
%     semilogy(alpha_vec,lambda_norm(k,:),'b.');
%     hold on
% end
ylabel('$\|\vec \lambda\|$','interpreter','latex')
xlabel('$\alpha$','interpreter','latex')
xticks([pi/6,pi/3])
xticklabels({'\pi/6','\pi/3'});

 
grid on
max(lambda_norm)

figure(3)
clf;
scatter(real(rin),imag(rin),10,log10(abs(stok_x+1i*stok_y)),'filled')
axis equal
colorbar
