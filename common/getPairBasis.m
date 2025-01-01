function [Ubf,Ybf,Ucf,Ycf,nimage] = getPairBasis(q,a_f,N_f,rads,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,s,N_peanut,L,Kf,project,Lc,opt)
%Get SVD factorisations of the fine grid for all identified close pairs of circles.
%Center coordinates q, global upsampling factor a_f for fine collocation
%points, N_f fine source points (for Stokeslets), circle radii rads (only unit radius considered
%so far), rbase_in_f the fine base grid, rimage_vec a cell array of image
%points for every contact, refine a cell array of parameter values t for
%the refined collocation points corresponding to a specific contact, the
%vector s sets the source types at image points (now hard coded to
%stresslets and potential dipoles). Compression with peanut points is done
%if N_peanut>0. L,K, used for the mobility problem. 


if nargin<11
    N_peanut = 0; 
    L = [];
    project = 0; 
elseif nargin < 12 
    L = [];
    project = 0; 
    Lc = []; 
end

P = length(q); 
Ubf = cell(P);
Ybf = cell(P);

%If using peanut compression 
if N_peanut
    Ucf = cell(P);
    Ycf = cell(P);
else
    Ucf = [];
    Ycf = [];
end

nimage = cell(P);

for i = 1:P

    if ~isempty(pairs)
        neigh = find(pairs(:,1)==i);
        rin_1_f = q(i)+rbase_in_f;
        %rin_1_c = rbase_in_f;
        for k = 1:length(neigh)
            %send in two particles with the fine grid 
            p2 = pairs(neigh(k),2);  
            rin_2_f = q(p2)+rbase_in_f;

            %% get grid for a pair 
           % [rin_f,rout_f,nimage,rimage] = getPairGrid([q(i),q(p2)],rads([i,p2]),Rg_c,N_c,a);
            rimage = [rimage_vec{i,p2}; rimage_vec{p2,i}];
            ntemp =  randn(size(rimage,1),2);
            nimage{i,p2} = (ntemp(1:end/2,1)+1i*ntemp(1:end/2,2))./abs(ntemp(1:end/2,1)+1i*ntemp(1:end/2,2));
            nimage{p2,i} = (ntemp(end/2+1:end,1)+1i*ntemp(end/2+1:end,2))./abs(ntemp(end/2+1:end,1)+1i*ntemp(end/2+1:end,2));

            nout = ceil(a_f*N_f); 
          %  nout = ceil(a2*(N_c+size(rimage,1))); %Previous choice
            t = linspace(0,2*pi,nout+1);
            t = t(1:end-1)';
            rout_base = rads(k)*(cos(t)+1i*sin(t));
            t1 = refine{i,p2};
            fine_1 = q(i)+rads(k)*(cos(t1)+1i*sin(t1));
            t2 = refine{p2,i};
            fine_2 = q(p2)+rads(k)*(cos(t2)+1i*sin(t2));

            rout_f = [q(i)+rout_base; fine_1 ;q(p2)+rout_base; fine_2];
            rin_pair = [rin_1_f; rin_2_f];
            
            %% Projection (Mobility only)
            if project
                B1 = getKmat2D([q(i)+rout_base; fine_1],q(i));
                B2 = getKmat2D([q(p2)+rout_base; fine_2],q(p2));
                Lr = getLrPair(B1,B2,Kf);            
            else
                Lr = []; 

            end
            
            %% Compute fine basis pseudo inverse

            %Takes in fine grid with image enhancement and fine collocation
            %points. 
            [Up,Yp,A] = getPairBlock(rin_pair,rout_f,[nimage{i,p2}; nimage{p2,i}],rimage,s,L,Lr,project);

%             N_c = opt.N_c;
%             tau_mapped = rand(2*opt.N_c,1);
%             tau_mapped(end/2+1:end) = 0;
%             mapped = rand(2*opt.N_c,1);
%             mapped(end/2+1:end) = 0;

            if opt.precomp
%
                 Nother1 = singleLayer(rbase_in_c+q(i),rout_f(end/2+1:end),1); %points on p2
                 Nother2 = singleLayer(rbase_in_c+q(p2),rout_f(1:end/2),1); %points on i
    % 
                 Npair = getNPair(Nother2,Nother1);
    % 
    %             % DEBUG
    %             rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,i,p2); 
    %             %Nother = singleLayer(rbase_in_c+q(i),rout_fine_other,mu);
    %             %R2 = -Nother*tau_mapped; %read off on particle 2
    % 
    %             [u2,v2] = StokesletDirect(real(rbase_in_c+q(i)),imag(rbase_in_c+q(i)),...
    %                 real(rout_fine_other),imag(rout_fine_other),...
    %                 tau_mapped(1:N_c),tau_mapped(N_c+1:2*N_c),N_c);
    %             R2  = -[u2; v2];
    % 
    %             %% Do a similar thing for the other order of the particles in the pair
    %             rout_fine_other = getFineOther(opt.a_f,opt.N_f,opt.rads,refine,q,p2,i); 
    %             %Nother2 = singleLayer(rbase_in_c+q(p2),rout_fine_other,mu);
    %             %R1 = -Nother2*mapped; %read off on particle 1
    % 
    %             %To be replaced with C implementation?
    %             [u1,v1] = StokesletDirect(real(rbase_in_c+q(p2)),imag(rbase_in_c+q(p2)),...
    %                 real(rout_fine_other),imag(rout_fine_other),...
    %                 mapped(1:N_c),mapped(N_c+1:2*N_c),N_c);
    %             R1 = -[u1; v1];
    % 
    %         
    %             rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)]; 
    % 
    % 
    %             rhs2 = -Npair*[tau_mapped(1:end/2); mapped(1:end/2); tau_mapped(end/2+1:end); mapped(end/2+1:end)];
    
    
    
                %Npair = [Nother1(:,1:end/2) Nother2(:,1:end/2) Nother1(:,end/2+1:end) Nother2(:,end/2+1:end)];
    
               % Ubf{i,p2} = Up;
                Ubf{i,p2} = -Up*Npair;
            else
                Ubf{i,p2} = Up;
            end

            Ybf{i,p2} = Yp; 

            if N_peanut %If true: compress on peanut described by N_peanut points
                
                debug = 0; 
                rout_peanut = create_peanut(q(i),q(p2),N_peanut,debug);    
                rin_pair_c = [q(i)+rbase_in_c; q(p2)+rbase_in_c];
                [DC,YC] = getPeanutBlock(rin_pair_c,rin_pair,rout_peanut,[nimage{i,p2}; nimage{p2,i}],rimage,s,A,Lc); 
                          

                Ucf{i,p2} = DC;
                Ycf{i,p2} = YC;

            end
             
    
        end
        
    
    end
    

end

end





function [U,Y,A] = getPairBlock(rin_pair,rout_f,nimage,rimage,s,L,Lr,project)

mu = 1; 

% figure(8)
% clf; 
% plot(real(rin_pair),imag(rin_pair),'r+')
% hold on
% plot(real(rout_f),imag(rout_f),'bo')
% plot(real(rimage),imag(rimage),'k*')


if size(rimage,1)
    Nimage = getImageKernels2D(rimage,nimage,rout_f,mu,s);
else
    Nimage = [];
end
N = singleLayer(rin_pair,rout_f,mu); 

if project
    %There are four blocks in L, xx xy, yx, yy...
    %Need to build projection for the pair, A, the same way

    A = getLpair(L);
    N = [N*A+Lr Nimage];

else
    N = [N Nimage]; 
    A = []; 
end 

[UU,S,V] = svd(N); 
S = diag(S);

%tol = 1e-6;
tol = 1e-14; 
tol = 1e-11; %4 nov: seems to work well but results in peaks in close to touching points
%tol = 1e-9; 
%tol = 1e-8; 
tol = 1e-12; 
tol = 1e-10; 

%tol = 1e-14; %Good level of truncation it seems, for close contacts. Tried with a line of circles.
%tol = eps; % I don't think this is reasonable. Just to understand what happens now. 

%use relative tolerance 
ra = sum(S>max(S)*tol); 
S = S(1:ra);  
iS = 1./S; 
Y = V(:,1:ra)*diag(iS); 
U = UU(:,1:ra)'; 


end