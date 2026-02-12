function viewPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,L,Lf,Nc,Nf,a_c,a_f,rads)
%Take U,Y as input, or recompute for an even coarser grid?

%% Evaluation points
%create grid to check basis in, at the boundaries of the particles in the
%pair
%parameterisation to visualise in at boundary
tself = linspace(0,2*pi,a_c*Nc+1); 
tself = tself(1:end-1)'; 
teval = linspace(0,2*pi,50*a_c*Nc+1); 
teval = teval(1:end-1)'; 
rcheck_b = [q(1)+cos(teval)+1i*sin(teval); q(2)+cos(teval)+1i*sin(teval)];

% We also want to evaluate at a third particle, centered at x,y, touching
% the other two, or, alternatively, close by. Evaluation points in rcheck3
% and rcheck4. 
delta = abs(q(1)-q(2))-2; 
alpha = acos((2+delta)/4);
x = 1+delta/2;
y = 2*sin(alpha);
rcheck3 = x+cos(teval)+y*1i+1i*sin(teval);
%rcheck4 = x+cos(tfine)+(y+delta)*1i+1i*sin(tfine); offset y for third
%particle in triangle
rcheck4 = x+cos(teval)+(y+delta)*1i+1i*sin(teval); 

%Also look in a bunch of points exterior to the particles
N = 600;
t = linspace(0,2*pi,N);
t = t(1:end-1)';
rcheck = [q(1)+cos(t)+1i*sin(t)];% q(2)+cos(t)+1i*sin(t)];
N = 200; 
x = linspace(-2,4,N);
y = linspace(-2,2,N);
x = linspace(0,2,N);
y = linspace(-1,1,N);
[XX,YY] = meshgrid(x,y);
r = XX(:)+1i*YY(:);
ind1 = find(abs(r-q(1))<1);
ind2 = find(abs(r-q(2))<1);
%ind2 = []; 
indok = setdiff(1:N*N,ind1);
indok = setdiff(indok,ind2); 
rcheck = [rcheck; r(indok)]; 

%later use surfir for visualisation 

% figure(9)
% clf; 
% plot(real(q(1))+cos(tfine),imag(q(1))+sin(tfine))
% hold on
% plot(x+cos(tfine),y+sin(tfine))
% plot(real(q(2))+cos(tfine),imag(q(2))+sin(tfine))
% plot(real(rcheck3),imag(rcheck3),'b*')

%% Tools needed to determine the pair-basis

%Initialize data on boundary
mu = zeros(2*a_c*Nc,1);
i = 1; %particle 1
p2 = 2; %particle 2

%relevant image indices for stokeslets
im_nr = length(rimage_vec{i,p2});
s_ind1_x = 1:Nf; 
s_ind2_x = Nf+1:2*Nf;
s_ind1_y = 2*Nf+1:3*Nf;
s_ind2_y = 3*Nf+1:4*Nf;

%The ones below are different depending on delta and therefore we
%store them in a different way. This is the bookkeeping done in the
%matvec. Would be better not to copy this (and use same as in matvec)... 

%stresslet indices
t_ind_x = 4*Nf+1:4*Nf+2*im_nr;
t_ind_y = 4*Nf+2*im_nr+1:4*Nf+4*im_nr;

%potential dipole indices
p_ind_x = 4*Nf+4*im_nr+1:4*Nf+6*im_nr;
p_ind_y = 4*Nf+6*im_nr+1:4*Nf+8*im_nr;

%proxy source locations
rvec_in = [rbase_in_f+q(1); rbase_in_f+q(2)]; 

%vectors needed for the "stresslets" with random "normal" vector
nvec = [nimage{i,p2};nimage{p2,i}]; 

%Keep track of what coarse point on boundary that is nonzero
t = linspace(0,2*pi,a_c*Nc+1); 
t = t(1:end-1)';

%for visualisation, this is the parametrisation used for the coarse
%collocation points.
tcoarse = linspace(0,2*pi,a_c*Nc+1);
tcoarse = tcoarse(1:end-1)'; 


%these matrices do not change if density changes
[rout_fine_other,tother] = getFineOther(a_f,Nf,rads,refine,q,i,p2);  
rout_fine_self = getFineOther(a_f,Nf,rads,refine,q,p2,i); 

%Needed to determine pair-basis
Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,1);

%Only used to check properties of the pair-basis
Nother_c = singleLayer(rbase_in_c+q(i),q(p2)+cos(teval)+1i*sin(teval),1);

tfiner = linspace(0,2*pi,1000+1);
tfiner = tfiner(1:end-1)';
Nother_finer = singleLayer(rbase_in_c+q(i),q(p2)+cos(tfiner)+1i*sin(tfiner),1);

%coarse self-interaction
Nself = singleLayer(rbase_in_c,cos(tcoarse)+1i*sin(tcoarse),1);
Nself = singleLayer(rbase_in_c,cos(teval)+1i*sin(teval),1);
Nself3 = singleLayer(rbase_in_c,rcheck3,1);
Nself4 = singleLayer(rbase_in_c,rcheck4,1);

Ncheck = singleLayer(rbase_in_c,rcheck,1); 

%% Check basis
% loop over x, then y, one nonzero coordinate in mu at the time
for k = 1 :2*a_c*Nc 
    try
        ind = mod(k,a_c*Nc+1);
        if floor(k/(a_c*Nc+1))
            ind = ind+1;
        end
        tstar = t(ind); %for visualisation in what point on coarse boundary data is non-zero
    catch
        disp('index error')
    end
    mu(:) = 0; %set coarse data
    mu(k) = 1;
    step1 = U{1}*mu; %Map to data at proxy points
    tau_mapped = Y{1}*step1;

    %If mobility, need to project
    if ~isempty(L)
        tau_mapped_non = tau_mapped(1:2*Nc);
        tau_i_x = tau_mapped(1:Nc);
        tau_i_y = tau_mapped(Nc+1:end);
        tau_mapped = [tau_i_x; tau_i_y]-L*[tau_i_x; tau_i_y];
    end

    %% Look first at the 1-body basis...
    % ...on the body itself. 
    uself = Nself*tau_mapped; 

    figure(6);
    clf;
%     plot(tself,uself(1:end/2))
%     hold on
%     plot(tself,uself(end/2+1:end))
    plot(teval,uself(1:end/2),'LineWidth',2)
    hold on
    plot(teval,uself(end/2+1:end),'LineWidth',2)
    plot(tcoarse,zeros(a_c*Nc,1),'k-*')
    xlabel('t')
    legend('x','y','interpreter','latex')
    title('1-body basis on self')
    xlim([0,2*pi])

    figure()
    f1 = fft(uself(1:end/2))
    
    %A smooth function, but it has to be resolved on the neighbour...

    uother = Nother*tau_mapped; 
    uother_f = Nother_finer*tau_mapped; 
    figure(12);
    clf;
    plot(tother,uother(1:end/2),'.')
    hold on
    plot(tother,uother(end/2+1:end),'.')
    plot(tfiner,uother_f(1:end/2),'k')
    hold on
    plot(tfiner,uother_f(end/2+1:end),'k')
    plot(tcoarse,zeros(a_c*Nc,1),'k-*')
    xlabel('t')
    legend('x','y','interpreter','latex')
    title('1-body basis evaluated on fine grid of neighbour')
    axis tight

    %... Check in the exterior of the particle
    ufield = Ncheck*tau_mapped;

    figure(13)
    subplot(1,3,1)
    surfir(real(rcheck),imag(rcheck),log10(abs(ufield(1:end/2))))
    %surfir(real(rcheck),imag(rcheck),ufield(1:end/2));
    colorbar
    axis tight
    view(0,90)
    subplot(1,3,2)
    surfir(real(rcheck),imag(rcheck),log10(abs(ufield(end/2+1:end))))
    %surfir(real(rcheck),imag(rcheck),ufield(end/2+1:end))
    removePatches([real(q(1)) imag(q(1))],0,1,1)
    view(0,90)
    colorbar
    axis tight

    figure(14)
    subplot(1,2,1)
    %surfir(real(rcheck),imag(rcheck),log10(abs(ufield(1:end/2))))
    surfir(real(rcheck),imag(rcheck),ufield(1:end/2));
    colorbar
    axis tight
    view(0,90)
    subplot(1,2,2)
    %surfir(real(rcheck),imag(rcheck),log10(abs(ufield(end/2+1:end))))
    surfir(real(rcheck),imag(rcheck),ufield(end/2+1:end))
    removePatches([real(q(1)) imag(q(1))],0,1,1)
    view(0,90)
    colorbar
    axis tight


    %% Compute pair basis. 
    %Determine contribution from the 1-body basis on the other particle.
    

    %First evaluate coarse basis on fine grid on the neighbour
    %Set bc for the pair-problem
    if ~isempty(L)
        B = getKmat2D(rout_fine_self,q(i)); %
        K = getKmat2D(rbase_in_c,q(i));
        ucomp = B*K'*tau_mapped_non;
    end


    R1 = -Nother*tau_mapped;     
    block = R1(1:end/2);

    %Add contribution from completion flow for mob?
    A2 = [zeros(size(block)); block; zeros(size(block)); R1(end/2+1:end)]; 
    %store as x x y y 

    %Evaluate pseudo-inverse from fine sources and collocation points.
    pair_mapped = Upf{i,p2}*A2; 
    tau_mapped_pair = Ypf{i,p2}*pair_mapped; 

    if ~isempty(L)
        %Compute (I-L)\beta
            tau_mapped_tot_xi = tau_mapped_pair(s_ind1_x);
            tau_mapped_tot_yi = tau_mapped_pair(s_ind1_y);

            tau_mapped_pair_noni = [tau_mapped_tot_xi; tau_mapped_tot_yi];

            tau_mapped_tot_xp2 = tau_mapped_pair(s_ind2_x);
            tau_mapped_tot_yp2 = tau_mapped_pair(s_ind2_y);

            tau_mapped_pair_nonp2 = [tau_mapped_tot_xp2; tau_mapped_tot_yp2];

            tau_mapped_i = [tau_mapped_tot_xi;tau_mapped_tot_yi]-Lf*[tau_mapped_tot_xi;tau_mapped_tot_yi];
            tau_mapped_p2 = [tau_mapped_tot_xp2;tau_mapped_tot_yp2]-Lf*[tau_mapped_tot_xp2;tau_mapped_tot_yp2];
            tau_mapped_pair(1:4*Nf) = [tau_mapped_i(1:end/2); tau_mapped_p2(1:end/2); tau_mapped_i(end/2+1:end); tau_mapped_p2(end/2+1:end)];
    
    end

    debug = 0; 

%     %bild matrix for pair to check if close to identity
%     if debug && ~isempty(L)
% 
% 
%         % Debug state for mobility
%         a2 = zeros(size(A2));
%     
%         warning('Should be evaluated in fine outer grid')
% 
%         for ii = 1:length(a2)
%             ii
%         
%             a2(:) = 0; 
%             a2(ii) = 1; 
%         
%             %Evaluate pseudo-inverse from fine sources and collocation points.
%             pair_mapped = Upf{i,p2}*a2; 
%             tau_mapped_pair = Ypf{i,p2}*pair_mapped; 
% 
%             
%         
%             if ~isempty(L)
%                 %Compute (I-L)\beta
%                     tau_mapped_tot_xi = tau_mapped_pair(s_ind1_x);
%                     tau_mapped_tot_yi = tau_mapped_pair(s_ind1_y);
%         
%                     tau_mapped_pair_noni = [tau_mapped_tot_xi; tau_mapped_tot_yi];
%                     tau_mapped_tot_xp2 = tau_mapped_pair(s_ind2_x);
%                     tau_mapped_tot_yp2 = tau_mapped_pair(s_ind2_y);
%         
%                     tau_mapped_i = [tau_mapped_tot_xi;tau_mapped_tot_yi]-Lf*[tau_mapped_tot_xi;tau_mapped_tot_yi];
%                     tau_mapped_p2 = [tau_mapped_tot_xp2;tau_mapped_tot_yp2]-Lf*[tau_mapped_tot_xp2;tau_mapped_tot_yp2];
%                     tau_mapped_pair(1:4*Nf) = [tau_mapped_i(1:end/2); tau_mapped_p2(1:end/2); tau_mapped_i(end/2+1:end); tau_mapped_p2(end/2+1:end)];
%                     
% 
%                     
%             end
% 
%             [ufmm,vfmm] = stokesSLPfmm([tau_mapped_pair(s_ind1_x); tau_mapped_pair(s_ind2_x)],...
%                 [tau_mapped_pair(s_ind1_y); tau_mapped_pair(s_ind2_y)],real(rvec_in),imag(rvec_in),real(rcheck_b),imag(rcheck_b),...
%                     0,5);
%             u_stress = getStresslets(tau_mapped_pair(t_ind_x),tau_mapped_pair(t_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],...
%                 rcheck_b,real(nvec),imag(nvec));
%             u_pot = getPotdip(tau_mapped_pair(p_ind_x),tau_mapped_pair(p_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],rcheck_b);
% 
%     
%             u_pair = [ufmm; vfmm]+u_stress+u_pot;
%             res(:,ii) = u_pair;
%         end
%     
%         
%     end


    
    %% Evaluate flow field only on the boundary of the two particles
    %Not really reasonable to use fmm here and in subsequent calls... 
    [ufmm,vfmm] = stokesSLPfmm([tau_mapped_pair(s_ind1_x); tau_mapped_pair(s_ind2_x)],...
        [tau_mapped_pair(s_ind1_y); tau_mapped_pair(s_ind2_y)],real(rvec_in),imag(rvec_in),real(rcheck_b),imag(rcheck_b),...
         0,5);
    u_stress = getStresslets(tau_mapped_pair(t_ind_x),tau_mapped_pair(t_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],...
        rcheck_b,real(nvec),imag(nvec));
    u_pot = getPotdip(tau_mapped_pair(p_ind_x),tau_mapped_pair(p_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],rcheck_b);

    ub = [ufmm; vfmm]+ u_stress+u_pot; 

    u1x = ub(1:end/4);
    u2x = ub(end/4+1:end/2);
    u1y = ub(end/2+1:3*end/4);
    u2y = ub(3*end/4+1:end); 
    
    if ~isempty(L)
        %add Lr*beta 
        B = getKmat2D(rcheck_b(1:end/2),q(i)); %
        K = getKmat2D(rbase_in_f,q(i));
        ucomp1 = B*K'*tau_mapped_pair_noni; 
        ucomp2 = B*K'*tau_mapped_pair_nonp2;

        u1x = u1x+ucomp1(1:end/2);
        u1y = u1y+ucomp1(end/2+1:end); 

        u2x = u2x+ucomp2(1:end/2);
        u2y = u2y+ucomp2(end/2+1:end); 


    end

    %evaluate 1-body basis in same points on body 2, I expect the sum with \chi 1,2
    % to be close to zero 

    u2_selfother = Nother_c*tau_mapped;  

    figure(7)
    clf;
    subplot(1,3,1)
%     plot(tfine,u1x);
%     hold on
%     plot(tfine,u1y);
    semilogy(teval,abs(u1x));
    hold on
    semilogy(teval,abs(u1y));
    plot(tstar,0,'k*','MarkerSize',10)
    title('On particle 1')
    xlabel('t')
    subplot(1,3,2)
    plot(teval,u2x);
    hold on
    plot(teval,u2y); 
    legend('x','y','interpreter','latex')
    title('On particle 2')
    xlabel('t')

    subplot(1,3,3)
    plot(real(rcheck_b),imag(rcheck_b));
    hold on
    plot(real(q(1))+cos(tstar),imag(q(1))+sin(tstar),'k*','MarkerSize',10)   
    sgtitle('Evaluation of pair-basis \chi 1,2')
   
    %Include also the contribution from the 1-body basis
    figure(9)
    clf;
    subplot(1,2,1)
    plot(teval,u2x+u2_selfother(1:end/2));
    xlabel('t')
    subplot(1,2,2)
    plot(teval,u2y+u2_selfother(end/2+1:end));
    sgtitle('\chi 1,2 at Particle 2 with added 1-body basis, close to zero?')
    %pause(2);

    %% Same thing, but evaluate on a third close particle
    %Remember that all stokeslets dont live on the same curve now! 
    
    [ufmm,vfmm] = stokesSLPfmm([tau_mapped_pair(s_ind1_x); tau_mapped_pair(s_ind2_x)],...
        [tau_mapped_pair(s_ind1_y); tau_mapped_pair(s_ind2_y)],real(rvec_in),imag(rvec_in),real(rcheck3),imag(rcheck3),...
         0,5);
    u_stress = getStresslets(tau_mapped_pair(t_ind_x),tau_mapped_pair(t_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],...
        rcheck3,real(nvec),imag(nvec));
    u_pot = getPotdip(tau_mapped_pair(p_ind_x),tau_mapped_pair(p_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],rcheck3);

    uthird = [ufmm; vfmm]+ u_stress+u_pot; 

    [ufmm,vfmm] = stokesSLPfmm([tau_mapped_pair(s_ind1_x); tau_mapped_pair(s_ind2_x)],...
        [tau_mapped_pair(s_ind1_y); tau_mapped_pair(s_ind2_y)],real(rvec_in),imag(rvec_in),real(rcheck4),imag(rcheck4),...
         0,5);
    u_stress = getStresslets(tau_mapped_pair(t_ind_x),tau_mapped_pair(t_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],...
        rcheck4,real(nvec),imag(nvec));
    u_pot = getPotdip(tau_mapped_pair(p_ind_x),tau_mapped_pair(p_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],rcheck4);

    uthird2 = [ufmm; vfmm]+ u_stress+u_pot; 

    figure(10);
    clf;
    subplot(1,2,1)
    x1 = uthird(1:end/2);
    y1 = uthird(end/2+1:end);
    x2 = uthird2(1:end/2);
    y2 = uthird2(end/2+1:end);
    plot(teval,x1);
    hold on
    plot(teval,y1);
    plot(teval,x2,'LineWidth',2);
    plot(teval,y2,'LineWidth',2);
    legend('touching particle, x coord','touching particle, y coord','slight offset, x coord','slight offset, y coord',...
        'interpreter','latex')
    axis tight
    
    subplot(1,2,2)
    f1 = fft(x1);
    g1 = fft(y1);
    f2 = fft(x2);
    g2 = fft(y2);
    semilogy(abs(f1(1:end/2)),'*');
    hold on
    semilogy(abs(g1(1:end/2)),'*');
    semilogy(abs(f2(1:end/2)),'*');
    semilogy(abs(g2(1:end/2)),'*');
    grid on
    title('well-resolved in coarse colloc points?')
    legend('touching particle, x coord','touching particle, y coord','slight offset, x coord','slight offset, y coord',...
        'interpreter','latex')
    xlabel('k','interpreter','latex')
    sgtitle('\chi 1,2 at boundary of a tentative third particle')
    axis tight
    %slow decay with coarse grid!

    %% Also, add evaluation of the coarse 1-body basis. Still nice and
    %smooth?
    u3_selfother = Nself3*tau_mapped;
    u4_selfother = Nself4*tau_mapped;

    uthird  = uthird + u3_selfother;
    uthird2 = uthird2 + u4_selfother; 

    figure(14);
    clf;
    x1 = uthird(1:end/2);
    y1 = uthird(end/2+1:end);
    x2 = uthird2(1:end/2);
    y2 = uthird2(end/2+1:end);
    plot(teval,x1);
    hold on
    plot(teval,y1);
    plot(teval,x2);
    plot(teval,y2);
  
    legend('touching particle, x coord','touching particle, y coord','slight offset, x coord','slight offset, y coord',...
        'interpreter','latex')

    title('\chi 1,2 + 1-body basis at boundary of a tentative third particle')

    x1 = u3_selfother(1:end/2);
    y1 = u3_selfother(end/2+1:end);
    x2 = u4_selfother(1:end/2);
    y2 = u4_selfother(end/2+1:end);
    figure(15);
    clf;
    subplot(1,2,1)
    plot(teval,x1);
    hold on
    plot(teval,y1);
    plot(teval,x2);
    plot(teval,y2);
    legend('touching particle, x coord','touching particle, y coord','slight offset, x coord','slight offset, y coord',...
        'interpreter','latex')
    axis tight
    subplot(1,2,2)
    f1 = fft(x1);
    g1 = fft(y1);
    f2 = fft(x2);
    g2 = fft(y2);
    semilogy(abs(f1(1:end/2)),'*');
    hold on
    semilogy(abs(g1(1:end/2)),'*');
    semilogy(abs(f2(1:end/2)),'*');
    semilogy(abs(g2(1:end/2)),'*');
    grid on
    sgtitle('1-body basis at boundary of a tentative third particle')
    legend('touching particle, x coord','touching particle, y coord','slight offset, x coord','slight offset, y coord',...
        'interpreter','latex')
    xlabel('k','interpreter','latex')
    axis tight
    %The decay here is limited by the magnitude of tau_mapped (large!)


    %% Same thing again... Evaluate flow field in all points 
    [ufmm,vfmm] = stokesSLPfmm([tau_mapped_pair(s_ind1_x); tau_mapped_pair(s_ind2_x)],...
        [tau_mapped_pair(s_ind1_y); tau_mapped_pair(s_ind2_y)],real(rvec_in),imag(rvec_in),real(rcheck),imag(rcheck),...
         0,5);
    u_stress = getStresslets(tau_mapped_pair(t_ind_x),tau_mapped_pair(t_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],...
        rcheck,real(nvec),imag(nvec));
    u_pot = getPotdip(tau_mapped_pair(p_ind_x),tau_mapped_pair(p_ind_y),[rimage_vec{i,p2};rimage_vec{p2,i}],rcheck);

    utot = [ufmm; vfmm]+ u_stress+u_pot; 

    %% Visualise flow field in the exterior domain
    figure(5)
    clf; 
    subplot(1,2,1)
    %surfir(real(rcheck),imag(rcheck),utot(1:end/2))
    surfir(real(rcheck),imag(rcheck),log10(abs(utot(1:end/2))))
    view(0,90)
    hold on
    removePatches([real(q(1:2)) imag(q(1:2))],0,rads(1:2),10)
    %caxis([-8,1])
    colorbar
    plot3(real(q(1))+0.9*cos(tstar),imag(q(1))+0.9*sin(tstar),10,'r*')
    subplot(1,2,2)
    %surfir(real(rcheck),imag(rcheck),utot(end/2+1:end))
    surfir(real(rcheck),imag(rcheck),log10(abs(utot(end/2+1:end))))
    hold on
    removePatches([real(q(1:2)) imag(q(1:2))],0,rads(1:2),10)
    plot3(real(q(1))+0.9*cos(tstar),imag(q(1))+0.9*sin(tstar),10,'r*')
    view(0,90)
   % caxis([-8,1])
    colorbar

    %pause(2);



end


end

function removePatches(q,h,rad,maxval)
t = linspace(0,2*pi,200);
for k = 1:size(q,1)
    r_range = linspace(0,rad(k)+10*h,2);
    [R,T] = meshgrid(r_range,t);
    [X,Y] = pol2cart(T,R);
    patch(X(:)+q(k,1),Y(:)+q(k,2),maxval*ones(size(X(:))),'w','EdgeColor', 'w');
end

end