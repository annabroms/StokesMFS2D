clear;
close all;
d = 1e-6;
d = 0.05; % Example value for d
A0 = log(1 + d);              % Exact value

%Reacreates test Alex does on p.31, Section 4.4 (Nearly singular quadrature
%and reparameterization) of his BIE book. 

Nvec = 10:10:600;
for k = 1:length(Nvec) 
    N = Nvec(k);

    % Plain PTR
    %N = 1000;                     % Choose N large enough for decent accuracy
    t = linspace(0, 2*pi, N+1);   % N+1 points to include endpoint
    t(end) = [];                  % Remove the duplicated 2pi point for periodicity
    f = @(t) log(sqrt((1 + d - cos(t)).^2 + sin(t).^2)) / (2*pi);
    A(k) = 2*pi*1/N*(sum(f(t)));         % Basic trapezoidal integration 
    %fprintf('d=%g: plain PTR N=%d \t abs err %.3g (%.3g sec)\n', d, N, abs(A(k) - A0), tim);
    Terr(k) = abs(A(k) - A0);
    
    % Graded PTR
    %Ng = 80;  % 80 is sufficient for d = 1e-6, gives error < 1e-15
    beta = log(3 / d);
    [x, w] = gradedptr(N, beta);
    Ag(k) = sum(w .* f(x));
    
    %fprintf('\tgraded PTR      Ng=%d \t abs err %.3g\n', Ng, abs(Ag(k) - A0));
    Gerr(k) = abs(Ag(k) - A0);
end

figure()
semilogy(Nvec,Terr,'b.-')
hold on
semilogy(Nvec,Gerr,'r.-')

%%
% How many points are needed?
delta = 1e-1;

beta = log(3/sqrt(delta));
N = 20*round(beta);
[x, w] = gradedptr(N, beta);
x = real(x); 
figure();
plot(cos(x),sin(x),'.')

figure()
semilogy(abs(diff(x)))

% Taken from get2DImagegrid
% N_f = 150;
% sample_nbr = 30; 
% nout = N_f+2*2*(sample_nbr+1);
% a_loc = 6;
% h = 2*pi/nout/a_loc; 
% tstar = 0; 
% strip = logspace(-6,log10(4*h*sample_nbr),10*sample_nbr);
% t = [tstar-strip tstar+strip];
% figure()
% plot(cos(t),sin(t),'.')
% axis equal