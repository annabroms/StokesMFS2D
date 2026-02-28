function Npair = evaluateCoarseOnPair(q,rbase_in_c,rout_f)
%EVALUATECOARSEONPAIR Pair map evaluating coarse proxy sources on fine grid
% for a pair of particles. 
%
% Syntax:
%   Npair = evaluateCoarseOnPair(q,rbase_in_c,rout_f)
%   evaluateCoarseOnPair()   % run built-in self-test
%
% Description:
%   Builds the matrix used when forming pair-correction right-hand sides.
%   The map includes only cross-interactions:
%   body 1 sources are evaluated on body 2 targets, and body 2 sources are
%   evaluated on body 1 targets.
%
% Inputs:
%   q          - [q1 q2], complex body centers
%   rbase_in_c - complex coarse proxy source locations relative to center
%   rout_f     - complex pair targets as [targets_body1; targets_body2]
%
% Output:
%   Npair      - map in pair ordering [x1; x2; y1; y2] -> [ux1; ux2; uy1; uy2]
%
% Notes:
%   rout_f is split in half into body 1 and body 2 targets.
%
% Anna Broms, Updated Feb 2026 

if nargin == 0
    run_self_test();
    if nargout > 0
        Npair = [];
    end
    return;
end

if nargin < 3
    error('evaluateCoarseOnPair requires q, rbase_in_c, and rout_f.');
end

q = q(:);
rbase_in_c = rbase_in_c(:);
rout_f = rout_f(:);

if numel(q) ~= 2
    error('q must contain exactly two complex centers [q1 q2].');
end
if mod(numel(rout_f),2) ~= 0
    error('rout_f must have even length so it can be split into two bodies.');
end

nhalf = numel(rout_f)/2;
rout_1 = rout_f(1:nhalf);
rout_2 = rout_f(nhalf+1:end);

% Cross interactions only.
N21 = singleLayer(rbase_in_c + q(1),rout_2,1); % body 1 -> body 2
N12 = singleLayer(rbase_in_c + q(2),rout_1,1); % body 2 -> body 1

% Reorder [x;y] block structure into pair ordering [x1;x2;y1;y2].
xx12 = N12(1:end/2,1:end/2);
xx21 = N21(1:end/2,1:end/2);
yy12 = N12(end/2+1:end,end/2+1:end);
yy21 = N21(end/2+1:end,end/2+1:end);
yx12 = N12(end/2+1:end,1:end/2);
yx21 = N21(end/2+1:end,1:end/2);
xy12 = N12(1:end/2,end/2+1:end);
xy21 = N21(1:end/2,end/2+1:end);

Nxx = [zeros(size(xx12)) xx12; xx21 zeros(size(xx12))];
Nyy = [zeros(size(yy12)) yy12; yy21 zeros(size(yy12))];
Nxy = [zeros(size(xy12)) xy12; xy21 zeros(size(xy12))];
Nyx = [zeros(size(yx12)) yx12; yx21 zeros(size(yx12))];

Npair = [Nxx Nyx; Nxy Nyy];

end

function run_self_test()

fprintf('evaluateCoarseOnPair self-test:\n');

q = [0; 2.01];
Nc = 60;
Nt = 1.2*Nc;

tc = linspace(0,2*pi,Nc+1)'; tc(end) = [];
t1 = linspace(0,2*pi,Nt+1)'; t1(end) = [];

rbase_in_c = 0.55*(cos(tc) + 1i*sin(tc));
rout_fine_1 = q(1) + cos(t1) + 1i*sin(t1);
rout_fine_2 = q(2) + cos(t1) + 1i*sin(t1);
rout_f = [rout_fine_1; rout_fine_2];

Npair = evaluateCoarseOnPair(q,rbase_in_c,rout_f);

tau1 = randn(2*Nc,1);
tau2 = randn(2*Nc,1);

Nother_2 = singleLayer(rbase_in_c + q(1),rout_fine_2,1);
R2 = -Nother_2*tau2; % read off on particle 2

Nother_1 = singleLayer(rbase_in_c + q(2),rout_fine_1,1);
R1 = -Nother_1*tau1; % read off on particle 1

rhs = [R1(1:end/2); R2(1:end/2); R1(end/2+1:end); R2(end/2+1:end)];
rhs2 = -Npair*[tau2(1:end/2); tau1(1:end/2); ...
               tau2(end/2+1:end); tau1(end/2+1:end)];

% rhs and rhs2 should be the same.
rel_err = norm(rhs-rhs2,inf)/max(1,norm(rhs,inf));
fprintf('  rhs consistency rel. err : %1.3e\n',rel_err);

% Optional quick sanity check on one own-body block remaining zero.
u_own = Npair*[tau2(1:Nc); zeros(Nc,1); tau2(Nc+1:end); zeros(Nc,1)];
own_err = norm(u_own([1:Nt,2*Nt+1:3*Nt]),inf);
fprintf('  own-body zero check      : %.3e\n',own_err);


end
