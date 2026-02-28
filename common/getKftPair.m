function Kft_pair = getKftPair(K1,K2)
%GETKFPAIR Build pair map from pair-ordered source strengths to per-body RBM moments.
%
% Syntax:
%   Kft_pair = getKftPair(K1,K2)
%   getKftPair()   % run built-in self-test
%
% Description:
%   Assembles the 2-body map using two per-particle K matrices
%   (typically from getKmat2D). The returned matrix maps the pair-ordered
%   source vector
%
%       [f_x1; f_x2; f_y1; f_y2]
%
%   to
%
%       [F_x1; F_y1; T_1; F_x2; F_y2; T_2],
%
%   where [F_x; F_y; T] = K'*[f_x; f_y] for each particle.
%
% Inputs:
%   K1, K2 - Per-particle K matrices of size (2*N1)-by-3 and (2*N2)-by-3.
%
% Output:
%   Kf_pair - 6-by-(2*(N1+N2)) pair map in ordering [x1; x2; y1; y2].
%
% Anna Broms, Feb 2026

if nargin == 0
    run_self_test();
    return;
end

n1 = size(K1,1)/2;
n2 = size(K2,1)/2;

K1t = K1';
K2t = K2';

K1x = K1t(:,1:n1);
K1y = K1t(:,n1+1:end);
K2x = K2t(:,1:n2);
K2y = K2t(:,n2+1:end);

Kft_pair = [K1x, zeros(3,n2), K1y, zeros(3,n2); ...
           zeros(3,n1), K2x, zeros(3,n1), K2y];

end

function run_self_test()

q1 = 0.0 + 0.0i;
q2 = 2.2 - 0.3i;

N1 = 23;
N2 = 29;
t1 = linspace(0,2*pi,N1+1)'; t1(end) = [];
t2 = linspace(0,2*pi,N2+1)'; t2(end) = [];

r1 = q1 + (1.0 + 0.12*cos(2*t1)).*(cos(t1) + 1i*sin(t1));
r2 = q2 + (0.9 + 0.10*sin(3*t2)).*(cos(t2) + 1i*sin(t2));

K1 = getKmat2D(r1,q1);
K2 = getKmat2D(r2,q2);

Kft_pair = getKftPair(K1,K2);


f_pair = randn(2*(N1+N2),1); % [x1; x2; y1; y2]
ind_x1 = 1:N1;
ind_x2 = N1+1:N1+N2;
ind_y1 = N1+N2+1:N1+N2+N1;
ind_y2 = N1+N2+N1+1:2*(N1+N2);

f1 = [f_pair(ind_x1); f_pair(ind_y1)];
f2 = [f_pair(ind_x2); f_pair(ind_y2)];

mom_pair = Kft_pair*f_pair;
mom_ref = [K1'*f1; K2'*f2];

rel_err = norm(mom_pair - mom_ref,inf)/max(1,norm(mom_ref,inf));

fprintf('getKfPair self-test:\n');
fprintf('  map consistency rel. err : %.3e\n',rel_err);


end
