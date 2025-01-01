function [u,lambda0] = getRecompletionFlow(rin,rout,q,f,t)

P = length(q); 
N_small = length(rin)/P; 

Ksmall = getKmat2D(rin(1:N_small),q(1));
A_s = Ksmall'*Ksmall;

%FT = zeros(3*P,1);
lambda0 = zeros(2*N_small*P,1); 
for k = 1:P
    F_k = [f(k,:)'; t(k)];
    lambda0k = Ksmall*(A_s\F_k);
    lambda0((k-1)*2*N_small+1:k*2*N_small) = lambda0k;
    lambda_x((k-1)*N_small+1:k*N_small) = lambda0k(1:end/2);
    lambda_y((k-1)*N_small+1:k*N_small) = lambda0k(end/2+1:end);

    %FT(3*(k-1)+1:3*k) = [f(k,:)'; t(k)]; 
end

% Kin = getKmat2D(rin,q); %large... Better to build for each body. Ok if all are the same.
% A  = Kin'*Kin;
% lambda0 = Kin*(A\FT);

%Kin sorts per particle... 
% for k = 1:P
%     lambda_x((k-1)*N_small+1:k*N_small) = lambda0((k-1)*2*N_small+1:(k-1)*2*N_small+N_small);
%     lambda_y((k-1)*N_small+1:k*N_small) = lambda0((k-1)*2*N_small+N_small+1:k*2*N_small);
% end


%Now, compute flow.

%debugging...
%Nii = singleLayer(rin,rout,1);
%u2 = Nii*[lambda_x';lambda_y'];
[ufmm,vfmm] = stokesSLPfmm(lambda_x,lambda_y,real(rin),imag(rin),real(rout),imag(rout),...
    0,5);

u = -[ufmm; vfmm]; %minus sign as we move this term to rhs




end