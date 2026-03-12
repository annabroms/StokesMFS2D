%Visualise decay in fourier coeffs for the evaluation of the pair
%correction on a third neighbouring particle. 
clear; 
close all

delta = 1e-2;
q = [0; 2+delta];
N_c = 50; 
N_f = 200; 
a_f = 1.2;

tol_c = 1e-12; %I think this works reasonably
tol_c = 1e-8; %Curve moves closer to the surface -> smaller coeff 
tol_c = 1e-16; %Curve moves further from surface -> larger coeff.
%tol_c = 1e-19;% Maybe better as a lower a_c would be needed for the same
%decay? 
%tol_c = 1e-23; 
% But smoother coarse 1-body basis to evaluate on neighbour

s = [0 0 1 1]; %set type of singularities at image points
%s = [1 0 1 1]; %Other singularities? Currently not supported! But code can
%be changed!

sep_f = (1/N_f)*log(1/tol_c); %what to pick?

Rg_f = max([1-sep_f,0.01]);  % and fine grid
Rg_f = 0.4; 


P = length(q); 
NN = 20; 
a1vec = linspace(1.2,10,NN+1);
Rgvec = linspace(0.1,0.9,NN); 
rad = ones(2,1);  

for i = 1:length(Rgvec)
    Rg_c = Rgvec(i); 
    i
    for k = 1:length(a1vec)
        k
        a_c = a1vec(k);

        %accumulation point, given Rg and delta. Closed formula from fixed point of reflection formula
        accstop = (1-Rg_c)^2/Rg_c;  
        delta_pair = accstop; 

        %% CREATE GRID
        %Outer basic grid
        tout_c = linspace(0,2*pi,ceil(a_c*N_c)+1);
        tout_c = tout_c(1:end-1)';
        
        rbase_out_c = rad(1)*cos(tout_c)+1i*rad(1)*sin(tout_c);
        
        %Construct image grid
        basic = 1; %return only the basic outer grid, else refined outer grid 
        [rout,~,~,~,~,pairs,rimage_vec,refine,rbase_in_f] = get2DImageGrid(q,rad,Rg_c, a_c, N_c, 1, Rg_f,a_f,N_f, basic,delta_pair);

        tin = linspace(0,2*pi,N_c+1);
        tin = tin(1:end-1)';
        rbase_in_c = Rg_c*cos(tin)+Rg_c*1i*sin(tin);
        
        %Get pair basis
        opt_pair = struct();
        opt_pair.N_f = N_f;
        opt_pair.a_f = a_f;
        opt_pair.rad = rad;
        opt_pair.s = s;
        opt_pair.N_peanut = 0;
        opt_pair.precomp = 1;
        opt_pair.proj_all = 0;
        [Upf,Ypf,~,~,~,~,nimage] = getPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,refine,pairs,opt_pair,[],[],[]);
        
        
        %Get one-body pseduo inverse blocks -- enough to do this for single body.
        [U,Y] = getSelfPseudo(rbase_in_c,rbase_out_c);
        
        
        [f1,g1,f2,g2] = testPairBasis(q,rbase_in_c,rbase_in_f,rimage_vec,nimage,refine,Upf,Ypf,U,Y,[],[],N_c,N_f,a_c,a_f,rad);
        try
            m1(k,i) = f1(end/2);
            m2(k,i) = g1(end/2);
            m3(k,i) = f2(end/2);
            m4(k,i) = g2(end/2); 
        catch
            ind = length(f1); 
            ind = (ind-1)/2;
            m1(k,i) = f1(ind);
            m2(k,i) = g1(ind);
            m3(k,i) = f2(ind);
            m4(k,i) = g2(ind);
        end
    end
end

%% Visualise proper parameter selection
[A,R] = meshgrid(a1vec,Rgvec); 
A = a1vec*N_c;
figure()
subplot(2,2,1)
surf(A,R,log10(m1'))
colorbar
view(0,90)
caxis([-13,-5])
shading interp
colormap(parula(8))
axis tight
ylabel('$R_c$','interpreter','latex')
xlabel('$M_c$','interpreter','latex')


subplot(2,2,2)
surf(A,R,log10(m2'))
colorbar
view(0,90)
caxis([-13,-5])
shading interp
colormap(parula(8))
axis tight
ylabel('$R_c$','interpreter','latex')
xlabel('$M_c$','interpreter','latex')

subplot(2,2,3)
surf(A,R,log10(m3'))
colorbar
view(0,90)
caxis([-13,-5])
shading interp
colormap(parula(8))
axis tight
ylabel('$R_c$','interpreter','latex')
xlabel('$M_c$','interpreter','latex')

subplot(2,2,4)
surf(A,R,log10(m4'))
colorbar
view(0,90)
caxis([-13,-5])
shading interp
colormap(parula(8))
axis tight
ylabel('$R_c$','interpreter','latex')
xlabel('$M_c$','interpreter','latex')
