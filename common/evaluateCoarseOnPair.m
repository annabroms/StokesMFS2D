function Npair = evaluateCoarseOnPair(q, rbase_in_c, rout_f)
%EVALUATECOARSEONPAIR Pair map evaluating coarse proxy sources on fine grid
% for a pair of particles.
%
% Syntax:
%   Npair = evaluateCoarseOnPair(q, rbase_in_c, rout_f)
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
% Anna Broms, Updated Feb 2026 | Optimized May 2026

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

    % Validate inputs
    if numel(q) ~= 2
        error('q must contain exactly two complex centers [q1 q2].');
    end
    if mod(numel(rout_f), 2) ~= 0
        error('rout_f must have even length so it can be split into two bodies.');
    end

    % Precompute split indices once
    nhalf   = numel(rout_f) / 2;
    rout_1  = rout_f(1:nhalf);
    rout_2  = rout_f(nhalf+1:end);

    % Compute cross-interaction matrices
    N21 = stokSLPmat(rbase_in_c(:) + q(1), rout_2, 1); % body 1 -> body 2
    N12 = stokSLPmat(rbase_in_c(:) + q(2), rout_1, 1); % body 2 -> body 1

    % Precompute block index ranges once
    nc  = size(N12, 2) / 2;   % number of source DOFs per component
    nt  = size(N12, 1) / 2;   % number of target DOFs per component

    r1  = 1:nt;       % target row block 1
    r2  = nt+1:2*nt;  % target row block 2
    c1  = 1:nc;       % source col block 1
    c2  = nc+1:2*nc;  % source col block 2

    % Preallocate output directly — avoids all intermediate matrices
    Npair = zeros(4*nt, 4*nc);

    % Fill Npair in pair ordering [x1;x2;y1;y2] -> [ux1;ux2;uy1;uy2]
    % Row blocks: ux1=r1, ux2=r2, uy1=r1+2*nt, uy2=r2+2*nt
    % Col blocks: x1=c1, x2=c2, y1=c1+2*nc, y2=c2+2*nc

    % --- xx blocks ---
    Npair(r1,  c2)        = N12(r1, c1);   % xx12: x2 src -> ux1
    Npair(r2,  c1)        = N21(r1, c1);   % xx21: x1 src -> ux2

    % --- yx blocks (x source -> uy target) ---
    Npair(r1+2*nt, c2)    = N12(r2, c1);   % yx12
    Npair(r2+2*nt, c1)    = N21(r2, c1);   % yx21

    % --- xy blocks (y source -> ux target) ---
    Npair(r1,  c2+2*nc)   = N12(r1, c2);   % xy12
    Npair(r2,  c1+2*nc)   = N21(r1, c2);   % xy21

    % --- yy blocks ---
    Npair(r1+2*nt, c2+2*nc) = N12(r2, c2); % yy12
    Npair(r2+2*nt, c1+2*nc) = N21(r2, c2); % yy21

end