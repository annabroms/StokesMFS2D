function [f1,g1,f2,g2] = testPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,L,Lf,Nc,Nf,a_c,a_f,rads)
%Take U,Y as input, or recompute for an even coarser grid?

%% Evaluation points
%create grid to check basis in, at the boundaries of the particles in the
%pair
%parameterisation to visualise in at boundary
teval = linspace(0,2*pi,a_c*Nc+1); 
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
rcheck4 = x+cos(teval)+(y+0.1)*1i+1i*sin(teval); 

%Also look in a bunch of points exterior to the particles
N = 100;
t = linspace(0,2*pi,N);
t = t(1:end-1)';
rcheck = [q(1)+cos(t)+1i*sin(t); q(2)+cos(t)+1i*sin(t)];
N = 100; 
x = linspace(-2,4,N);
y = linspace(-2,2,N);
[XX,YY] = meshgrid(x,y);
r = XX(:)+1i*YY(:);
ind1 = find(abs(r-q(1))<1);
ind2 = find(abs(r-q(2))<1);
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
mu = zeros(2*ceil(a_c*Nc),1);
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
Nother = stokSLPmat(rbase_in_c+q(i),rout_fine_other,1);

%Only used to check properties of the pair-basis
Nother_c = stokSLPmat(rbase_in_c+q(i),q(p2)+cos(teval)+1i*sin(teval),1);

tfiner = linspace(0,2*pi,1000+1);
tfiner = tfiner(1:end-1)';
Nother_finer = stokSLPmat(rbase_in_c+q(i),q(p2)+cos(tfiner)+1i*sin(tfiner),1);

%coarse self-interaction
Nself = stokSLPmat(rbase_in_c,cos(tcoarse)+1i*sin(tcoarse),1);
Nself3 = stokSLPmat(rbase_in_c,rcheck3,1);
Nself4 = stokSLPmat(rbase_in_c,rcheck4,1);

%% Check basis

mu(1) = 1;
step1 = U{1}*mu; %Map to data at proxy points
tau_mapped = Y{1}*step1;

%If mobility, need to project
if ~isempty(L)
    tau_mapped_non = tau_mapped(1:2*Nc);
    tau_i_x = tau_mapped(1:Nc);
    tau_i_y = tau_mapped(Nc+1:end);
    tau_mapped = [tau_i_x; tau_i_y]-L*[tau_i_x; tau_i_y];
end

%% Look first at the 1-body basis on the body itself. 
uself = Nself*tau_mapped; 



%A smooth function, but it has to be resolved on the neighbour...

uother = Nother*tau_mapped; 
uother_f = Nother_finer*tau_mapped; 



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

%     if ~isempty(L) %For mobility. Confusing... 
%         R1 = R1-ucomp; %want to cancel the one-body basis on the other particle 
% %          %R1 = -Nother*(tau_mapped+comp(2*N_coarse*(i-1)+1:2*i*N_coarse));
%     end



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

        tau_mapped_i = [tau_mapped_tot_xi;tau_mapped_tot_yi]-Lf*[tau_mapped_tot_xi;tau_mapped_tot_yi];
        tau_mapped_p2 = [tau_mapped_tot_xp2;tau_mapped_tot_yp2]-Lf*[tau_mapped_tot_xp2;tau_mapped_tot_yp2];
        tau_mapped_pair(1:4*Nf) = [tau_mapped_i(1:end/2); tau_mapped_p2(1:end/2); tau_mapped_i(end/2+1:end); tau_mapped_p2(end/2+1:end)];

end





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
    ucomp = B*K'*tau_mapped_pair_noni; 
    %ucomp = zeros(size(ucomp)); %debug; 
    %ucomp = ucomp/a/Nc;

    %u2_selfother = u2_selfother + ucomp;
    u1x = u1x+ucomp(1:end/2);
    u1y = u1y+ucomp(end/2+1:end); 

    %check identity


end

%evaluate 1-body basis in same points on body 2, I expect the sum with \chi 1,2
% to be close to zero 

u2_selfother = Nother_c*tau_mapped;  

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
x1 = uthird(1:end/2);
y1 = uthird(end/2+1:end);
x2 = uthird2(1:end/2);
y2 = uthird2(end/2+1:end);


f1 = abs(fft(x1));
g1 = abs(fft(y1));
f2 = abs(fft(x2));
g2 = abs(fft(y2));

%slow decay with coarse grid!










end

function removePatches(q,h,rad,maxval)
t = linspace(0,2*pi,200);
for k = 1:length(q)
    r_range = linspace(0,rad(k)+10*h,2);
    [R,T] = meshgrid(r_range,t);
    [X,Y] = pol2cart(T,R);
    patch(X(:)+q(k,1),Y(:)+q(k,2),maxval*ones(size(X(:))),'w','EdgeColor', 'w');
end

end