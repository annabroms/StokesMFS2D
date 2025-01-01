function [UW,lambda,it,err] = solve_mob_precond_images(q,F,T,rads,delta_pair,rcheck_dom,precond)
%Intended to be the same thing as solve_2D_precond, but based on images for
%the fine grid. Tests an MFS generalisation of the Cheng-Greengard idea. 
%Applies FMM for the evaluation of Stokeslets but direct summation for
%other source types. 


%% SET PARAMS
%GMRES params
maxit = 800; 
gmres_tol = 1e-6; 

gmres_tol = 1e-10; 

%Grid params
P = length(q); 
%Set coarse and fine grid. 
%Play with N_c, N_f, a (a_f). 
N_c = 100;  %100 better here? 
%N_c = 100; 
%N_c = 100; 
%N_c = 100; 
N_f = 100; 
%N_f = N_c; %debug
%N_c = 2; 
%N_c = 150; 

%N_c = 250; 
a_c = 1.2;
%a_c = 2; 
%a_c = 5;


%a_c = 4;
%a_c = 1.2; %results in horrible fourier decay if we 
           % evaluate on a particle close by in coarse colloc points 
%a = 2; 
a_f = 1.2; %upsampling factor for the fine grid
%a_f = 5; 
%a2 = 0.6; 
%a2 = a; 

tol_c = 1e-12; %I think this works reasonably
%tol_c = 1e-8; %Curve moves closer to the surface -> smaller coeff 
%tol_c = 1e-16; %Curve moves further from surface -> larger coeff. 
%tol_c = 1e-11; 
% But smoother coarse 1-body basis to evaluate on neighbour


s = [0 0 1 1]; %set type of singularities at image points
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

sep_c = (1/N_c)*log(1/tol_c);
sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rg_c = max([1-sep_c,0.01]); %radius of proxy surface for coarse grid
Rg_f = max([1-sep_f,0.01]);  % and fine grid
%Rg_f = Rg_c; %debug


%Rg_c = 0.2; %debug! 
%Rg_c = 0.01; 

%To be sent to solver/matvec later
opt.a_c = a_c; 
opt.a_f = a_f; 
opt.N_c = N_c;
opt.N_f = N_f; 
opt.rads = rads; 
opt.s = s; 
opt.precomp = 1; %faster if evaluation of one body basis on fine grid is compted only once. 
% %Less storage required.
 
%accumulation point, given Rg and delta. Closed formula from fixed point of reflection formula
accstop = (1-Rg_c)^2/Rg_c;  

if nargin < 5
    delta_pair = accstop; %We want to use the pair correction for all gaps smaller than delta_pair. (or accstop).
end

%% CREATE GRID
%Outer basic grid
tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
tout_c = tout_c(1:end-1)';

rbase_out_c = rads(1)*cos(tout_c)+1i*rads(1)*sin(tout_c);

%Construct image grid
basic = 1; %return only the basic outer grid.                               
[rout,~,~,~,~,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rads,Rg_c, a_c, N_c, 1,Rg_f,a_f,N_f,basic,delta_pair);


%% PREPARE PRECONDITIONING AND RHS
tin = linspace(0,2*pi,N_c+1);
tin = tin(1:end-1)';
rbase_in_c = Rg_c*cos(tin)+Rg_c*1i*sin(tin);
rin = [];
for k = 1:P
    rin = [rin; q(k)+rbase_in_c];
end

%get evaluation of lambda0
[urhs,comp] = getRecompletionFlow(rin,rout,q,F,T); 
rimage_in = []; 
nimage = [];
[U,Y,Lc] = getSelfPseudoMobility(rbase_in_c,rbase_out_c,q,rimage_in,nimage,1,[0,ceil(a_c*N_c)],s,1);

%K = getKmat2D(rbase_in_c,0);

Kf = getKmat2D(rbase_in_f,0);
Lf = Kf*((Kf'*Kf)\Kf'); %This is x y

%Get pair basis

[Upf,Ypf,~,~,nimage] = getPairBasis(q,a_f,N_f,rads,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,s,opt,0,Lf,Kf,1);



% Now, check pair basis up to the boundary. Is it nice and smooth?
%warning('Deactivate opt.precomp');
%viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,Lc{1},Lf,N_c,N_f,a_c,a_f,rads)
%% Construct chech boundaries
%rcheck_dom = 100+100i;

% Create new grid points, for which the accuracy of the solution is
% to be evaluated. 
rcheck_b = [];
n_bound = 803;
t = linspace(0,2*pi,n_bound)';
for k = 1:P
    rcheck_b = [rcheck_b; q(k)+rads(k)*(cos(t)+1i*sin(t))];
end

%% SOLVE SYSTEM
% Build the matrix to check it out

debug = 0;
if debug
    x = zeros(size(urhs));
    for k = 1:size(urhs,1)
        k
        x(:) = 0; 
        x(k) = 1; 
        uu = matvec_mob_pairprecond_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf);
        CC(:,k) = uu;
    end
    figure()
    imagesc(log10(abs(CC)))
    colorbar
    skeel(CC)
end



[tau,it,resvec,real_res] = helsing_gmres_mv(@(x) matvec_mob_pairprecond_images(x,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf),urhs,2*length(rout),maxit,gmres_tol,1,rout);
debug = 1; 

%Modify to build with krylov preconditioning
%[tau, e2, precond] = precond_gmres(@(x) matvec_2D_pairprecond_images_noacc(x,rbase_in_c,refine,rimage_vec,nimage,opt,rvec_out,q,U,Y,pairs,Upf,Ypf,s), fout, zeros(2*size(rvec_out,1),1), 2*size(rvec_out,1), gmres_tol, precond,debug);
%it = length(e2); 

if debug
      figure()
%     semilogy(e2);
      semilogy(resvec);
      title('GMRES convergence mobility, pair corrections', 'Interpreter','latex')
      u2 = matvec_mob_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rout,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf);
end 

%% POSTPROCESS
[~,~,~,~,tau_stokes_x,tau_stokes_y,tau_stokes_nonpx,tau_stokes_nonpy,tau_stress_all_x,tau_stress_all_y,tau_pot_all_x,tau_pot_all_y] = getMobPairTransformation(tau,rbase_in_c,rbase_in_f,refine,...
    rimage_vec,nimage,opt,rout,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf);

tau_image = [tau_stress_all_x; tau_stress_all_y; tau_pot_all_x; tau_pot_all_y];
% tau_proxy = [tau_stokes_x; tau_stokes_y]; %This is the fine density
tau_proxy = [tau_stokes_nonpx; tau_stokes_nonpy];
%And evaluate in new points rcheck_dom and rcheck_b

lambda = [tau_proxy; tau_image]; %not just for other things than visualisation

%% Get rigid body motion. 

%First due to all coarse sources
Kc = getKmat2D(rbase_in_c,0);
UW= zeros(3*P,1); 
for k= 1:P
    UW((k-1)*3+1:3*k) = -Kc'*[tau_stokes_nonpx((k-1)*N_c+1:k*N_c); tau_stokes_nonpy((k-1)*N_c+1:k*N_c)];
end

%Then, due to all fine sources.
has_neigh = sort(unique(pairs(:)));
for i = 1:length(has_neigh)
    k = has_neigh(i); 
    UW((k-1)*3+1:3*k) = UW((k-1)*3+1:3*k)-Kf'*[tau_stokes_nonpx((k-1)*N_f+1+P*N_c:k*N_f+P*N_c); 
        tau_stokes_nonpy((k-1)*N_f+1+P*N_c:k*N_f+P*N_c)];
end


%% CHECK RESIDUAL AT SURFACE
%Compute velocity at surface
B = getKmat2D(rcheck_b(1:n_bound)-q(1),0); %same for all particles
u_lhs = zeros(2*P*n_bound,1);
for k = 1:P  
    res = B*UW(3*(k-1)+1:3*k);
    u_lhs((k-1)*n_bound+1:k*n_bound) = res(1:end/2);
    u_lhs(P*n_bound+(k-1)*n_bound+1:P*n_bound+k*n_bound) = res(end/2+1:end); 
end

%Using representation
u_rhs = matvec_mob_pairprecond_images(tau,rbase_in_c,rbase_in_f,refine,rimage_vec,nimage,opt,rout,rcheck_b,q,U,Y,Lc{1},Lf,pairs,Upf,Ypf);

S_0 = getRecompletionFlow(rin,rcheck_b,q,F,T); 
u_rhs = u_rhs-S_0;  %Note! Sign here due to how we have defined the completion flow. 
                    %This is accordinng to the representation of the flow

%debug
% figure()
% subplot(1,2,1)
% plot(u_rhs);
% hold on
% plot(u_lhs);
% subplot(1,2,2)
% semilogy(abs(u_rhs-u_lhs));


disp('Surface residual')
diff_vec = u_rhs-u_lhs;
max_abs = max(abs(u_rhs(1:end/2)+1i*u_rhs(end/2+1:end)));
err_vec =  abs(diff_vec(1:end/2)+1i*diff_vec(end/2+1:end))/max_abs;
err = max(err_vec)


%Some visualisation stuff... 
visualise = 1; 
if visualise
    figure(9)
    subplot(2,2,1)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(diff_vec(1:end/2))),30,log10(abs(diff_vec(1:end/2))),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('error in x velocity')


    subplot(2,2,2)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(diff_vec(end/2+1:end))),30,log10(abs(diff_vec(end/2+1:end))),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('error in y velocity')

    % Visualise the actual velocity in the rhs and lhs
    subplot(2,2,3)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(u_rhs(1:end/2))),30,log10(abs(u_rhs(1:end/2))),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('x velocity rhs')

    subplot(2,2,4)
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(u_lhs(1:end/2))),30,log10(abs(u_lhs(1:end/2))),'filled')
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('x velocity lhs')


    sgtitle('Error on boundary mob pair corr','interpreter','latex')

    figure()
    scatter3(real(rcheck_b),imag(rcheck_b),log10(abs(err_vec)),30,log10(abs(err_vec)),'filled');
    colorbar
    axis equal
    view(0,90)
    grid off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('Relative residual mob pair corr','interpreter','latex')


    %% Visualise source strengths
    figure()
    semilogy(abs(lambda))
    title('Source strengths mobility with pair corr')
    
end

end


%Old stuff from resistance


function doPairBasisTest(Upf,Ypf,i,p2,q,U,Y,rbase_in_c,rbase_out_f,rpair_fine)
%Test pair basis with a smooth coarse density
warning('Check consistency with the matvec! NOT in use')
        mu = 1; 

        %generate density
        t = linspace(0,2*pi,ceil(1.2*size(rbase_in_c,1))+1);
        t = t(1:end-1)';
        tau_1 = sum([exp(1i*t),exp(2*1i*t),exp(3*1i*t),exp(4*1i*t),exp(4*1i*t)]*rand(5),2);
        tau_2 = sum([exp(1i*t),exp(2*1i*t),exp(3*1i*t),exp(4*1i*t),exp(4*1i*t)]*rand(5),2);

        tau_p1_x = real(tau_1);
        tau_p1_y = imag(tau_1); 
        tau_p2_x = real(tau_2);
        tau_p2_y = imag(tau_2);

        figure()
        plot(tau_p1_x);
        hold on
        plot(tau_p1_y);
        plot(tau_p2_x);
        plot(tau_p2_y);

        %Get A1:
        step1 = U{i}'*[tau_p1_x;tau_p1_y]; %here I assume x and y follow each other?
        tau_mapped = Y{i}*step1;
        Nother = singleLayer(rbase_in_c+q(i),rbase_out_f+q(p2),mu);
        R1 = -Nother*tau_mapped; %not correct. Should read off on particle 2
        block = R1(1:end/2);
    
        A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)]; 
        %A1 = [tau_particle_x; zeros(size(tau_particle_x)); tau_particle_y; zeros(size(tau_particle_x))];
        pair_mapped = Upf{i,p2}'*A2;

        %pair_mapped2 = Upf{i,p2}'*A2;
        tau_mapped = Ypf{i,p2}*pair_mapped;

        %Get A2
        step1 = U{p2}'*[tau_p2_x;tau_p2_y]; %here I assume x and y follow each other?
        mapped = Y{p2}*step1;
        Nother = singleLayer(rbase_in_c+q(p2),rbase_out_f+q(i),mu);
        R1 = -Nother*mapped;
        block = R1(1:end/2);
        A1 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
        %A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)];            
        %A2 = [R1(1:end/2); zeros(size(block)); R1(end/2+1:end); zeros(size(block))];
        pair_mapped = Upf{i,p2}'*A1;
        tau_mapped2 = Ypf{i,p2}*pair_mapped;

        tau_tot = tau_mapped+tau_mapped2; 

        %Now apply. Want to compute what this is exterior to the two
        %particles in the pair
        x = linspace(-2,5);
        y = linspace(-2,5);
        [X,Y] = meshgrid(x,y); 
        rcheck = X(:)+1i*Y(:);
        ind1 = find(abs(rcheck-q(i))<1);
        ind2 = find(abs(rcheck-q(p2))<1);
        ind_keep_1 = setdiff(1:size(rcheck,1),ind1);
        ind_keep_2 = setdiff(ind_keep_1,ind2);
        %rcheck = rcheck(ind_keep_2);

        Npair = singleLayer(rpair_fine,rcheck,mu);
        ucheck = Npair*tau_tot;
      %  ucheck(ind1) =  nan;
      %  ucheck(ind2) = nan; 
        ucheck = ucheck(1:end/2)+1i*ucheck(end/2+1:end);
        ucheck(ind1) = nan;
        ucheck(ind2) = nan;


        figure()
        surfir(real(rcheck),imag(rcheck),abs(ucheck));
        colorbar
        axis equal
        hold on
        plot(real(q),imag(q),'k+')
        view(0,90);

        figure()
        surfir(real(rcheck),imag(rcheck),log10(abs(ucheck)));
        colorbar
        axis equal
        hold on
        plot(real(q),imag(q),'k+')
        view(0,90);

        %Also, compute residuals from the solve steps above.
        %ucheck is the velocity computed with the pair basis?
        Npair = singleLayer(rpair_fine,[rbase_out_f+q(i); rbase_out_f+q(p2)],mu);
        ulhs = Npair*tau_mapped;
        urhs = A2;

        ulhs2 = Npair*tau_mapped2;
        urhs2 = A1;

        norm(ulhs-A2,inf) %seems pretty accurate... 
        norm(ulhs2-A1,inf)

%         figure(10)
%         clf;
%         plot((A1(1:end/4)-ulhs2(1:end/4))./A1(1:end/4))
%         hold on
%         plot((A1(end/2+1:3*end/4)-ulhs2(end/2+1:3*end/4))./A1(end/2+1:3*end/4))
%         plot((A2(end/4+1:end/2)-ulhs(end/4+1:end/2))./A2(end/4+1:end/2))
%         plot((A2(3*end/4+1:end)-ulhs(3*end/4+1:end))./A2(3*end/4+1:end))

        figure(10)
        clf;
        plot((A1(1:end/4)-ulhs2(1:end/4)))
        hold on
        plot((A1(end/2+1:3*end/4)-ulhs2(end/2+1:3*end/4)))
        plot((A2(end/4+1:end/2)-ulhs(end/4+1:end/2)))
        plot((A2(3*end/4+1:end)-ulhs(3*end/4+1:end)))


end



