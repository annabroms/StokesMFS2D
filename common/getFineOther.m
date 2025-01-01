function [rout_fine_other,tout] = getFineOther(a,N_f,rads,refine,q,i,p2)
%Retrieve fine grid on circle p2 close to circle i
nout = ceil(a*N_f); % number of uniformly distributed collocation points on fine grid.
t = linspace(0,2*pi,nout+1);
t = t(1:end-1)';
rout_base = rads(p2)*(cos(t)+1i*sin(t));

t2 = refine{p2,i}; %get parameterisation of the clustered points close to contact
fine_2 = q(p2)+rads(p2)*(cos(t2)+1i*sin(t2));

rout_fine_other = [q(p2)+rout_base; fine_2];

tout = [t; t2];


end