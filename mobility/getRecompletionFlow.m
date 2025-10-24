function [u,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T)
%getRecompletionFlow computes completion flow for a 2D mobility problem due to
%given forces and torques
%
% Syntax: [u,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T)
%
% Input: 
% rin       - a vector of complex valued source points, with an equal number
%             of sources on each particle.
% rout      - a vector of complex valued target points, in which to
%             evaluated the completion flow
% q         - vector of length P with the complex valued center coordinates for
%             the particles
% F         - matrix of size P x 2 representing the x y coordinates of the of
%           the net forces acting on the particles 
% T         -  column vector of size P representing the torque acting on the 
%           particles  
%
% Output: 
% lambda_x  - x component of the computed completion sources (due to the given net
%               forces/torques)
% lambda_y  - y component of the computed completion sources 
% u         - the completion flow evaluated at targets rout due to
%             completion sources at rin of magnitude lambda_x,lambda_y
%
% Note: The grid of proxy sources is assumed to be the same on every
% particle
%
% Anna Broms April 4, 2025


if nargin==0, test_compflow; return; end


%% Checks
P = length(q);
N_small = length(rin)/P; % Number of sources per particle

assert(size(T,1)==P,'Wrong size of torque vector')
assert(size(F,1)==P,'Wrong size of force vector')
assert(size(F,2)==2,'Wrong size of force vector, should contain x y coordinates')


Ksmall = getKmat2D(rin(1:N_small),q(1)); 
A_s = Ksmall'*Ksmall; %assumed the same for everyone

lambda_x = zeros(N_small*P,1);
lambda_y = zeros(N_small*P,1);


for k = 1:P
    F_k = [F(k,:)'; T(k)];
        
    lambda0k = Ksmall*(A_s\F_k);
    lambda_x((k-1)*N_small+1:k*N_small) = lambda0k(1:N_small);
    lambda_y((k-1)*N_small+1:k*N_small) = lambda0k(N_small+1:2*N_small);

end

%Now, compute flow.
res = getVelocityField(rin,rout,lambda_x,lambda_y,[],[],[],[],[],[],[]);

u = -res; %minus sign as we move this term to rhs


end

function test_compflow
%Test script to visualise the completion flow on a particle, together with
%the completion sources

q = 0; %single body at the origin.
F = [1 0]; %force on particle
F = [0 0]; 
T = 1; %torque on the particle
%Generate proxy sources
tin = linspace(0,2*pi,60);
tin = tin(1:end-1)';
Rg = 0.5;
rin = Rg*(cos(tin)+1i*sin(tin));

%generate points exterior to the particle
x = linspace(-3,3,20); 
[XX,YY] = meshgrid(x,x); 
rout =  XX(:)+1i*YY(:); 
rout = rout(abs(rout-q)>1+1e-4); 

%get completion flow and visualise.
[u,lambda_x,lambda_y] = getRecompletionFlow(rin,rout,q,F,T);

figure()
subplot(1,2,1)
quiver(real(rin),imag(rin),lambda_x,lambda_y,'k'); %visualise force density
axis equal
title('Completion sources')
subplot(1,2,2)
quiver(real(rout),imag(rout),u(1:end/2),u(end/2+1:end)); 
axis equal
title('Competion flow')


end